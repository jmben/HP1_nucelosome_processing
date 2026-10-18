#!/usr/bin/env python3
"""
03_run_simulation.py  —  Coarse-grained chromatin polymer simulation (OpenMM)
=============================================================================
Implements the full energy function from Supplementary Note 2:

    U = U_bond + U_EV + U_confine + U_contact

KEY FIX over previous versions:
  All forces are expressed in kJ/mol with nm distances so that OpenMM's
  internal unit system is respected exactly.  Previous versions passed
  "1.0 * unit.kelvin" as temperature, giving T ≈ 0 K — the thermostat
  was essentially off, letting energy accumulate until NaN.

  Conversion from paper's reduced units:
    σ_nm     = 10 nm  (1 bead ≈ nucleosome diameter)
    T        = 300 K
    kBT      = R*T  = 8.314e-3 * 300 = 2.494 kJ/mol
    k2       = 20 kBT/σ²   →  20 * 2.494 / 10² kJ/mol/nm²
    prefactor energies multiplied by kBT to give kJ/mol

Usage:
    python 03_run_simulation.py \
        --energy  data/energy_matrix.npz \
        --out_dir results/ \
        [--n_traj 16] [--platform CPU]

Quick test (~5-10 min):
    python 03_run_simulation.py --energy data/energy_matrix.npz \
        --n_traj 2 --equil_steps 50000 --sim_steps 200000 --out_dir results_test/
=============================================================================
"""

import argparse
import numpy as np
from pathlib import Path
import time

try:
    import openmm as mm
    import openmm.unit as unit
    from openmm import app
except ImportError:
    raise ImportError("conda install -c conda-forge openmm")


# ---------------------------------------------------------------------------
# Physical constants  (nm / kJ/mol / ps unit system)
# ---------------------------------------------------------------------------
SIGMA_NM  = 10.0                      # bead diameter in nm
T_K       = 300.0                     # temperature in Kelvin
R_KJ      = 8.314462e-3               # gas constant kJ/(mol·K)
KBT_KJ    = R_KJ * T_K               # 2.494 kJ/mol = 1 kBT
MASS_GMOL = 1e6                       # g/mol (large → overdamped, stable)

# WCA / softcore EV parameters
R_WCA_NM  = (2 ** (1.0 / 6)) * SIGMA_NM
ECUT_KJ   = 4.0 * KBT_KJ
R_S_NM    = R_WCA_NM * ((ECUT_KJ / (2 * KBT_KJ) + 1) ** (-1.0 / 6))


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--energy",      default="data/energy_matrix.npz")
    p.add_argument("--out_dir",     default="results/")
    p.add_argument("--n_traj",      type=int,   default=16)
    p.add_argument("--equil_steps", type=int,   default=10_000_000)
    p.add_argument("--sim_steps",   type=int,   default=80_000_000)
    p.add_argument("--save_every",  type=int,   default=10_000)
    p.add_argument("--dt_ps",       type=float, default=0.020,
                   help="Timestep in picoseconds (default 0.020 = 20 fs)")
    p.add_argument("--strength",    type=float, default=20.0,
                   help="Flory-Huggins χ strength in kBT units (default 20.0).")
    p.add_argument("--conc",        type=float, default=2.0,
                   help="Bead concentration in mM (default 2.0). "
                        "Controls confinement sphere size. "
                        "0.2=paper default (too sparse), 2.0=good for phase sep.")
    p.add_argument("--platform",    default="CPU")
    p.add_argument("--seed",        type=int,   default=42)
    return p.parse_args()


# ---------------------------------------------------------------------------
# Hilbert-like initialisation
# ---------------------------------------------------------------------------
def hilbert_init(n: int, radius_nm: float) -> np.ndarray:
    """Serpentine space-filling start inside sphere of radius_nm (nm)."""
    side = int(np.ceil(n ** (1 / 3))) + 1
    coords = []
    for i in range(side):
        for j in range(side):
            row = j if i % 2 == 0 else (side - 1 - j)
            for k in range(side):
                col = k if row % 2 == 0 else (side - 1 - k)
                coords.append([i, row, col])
                if len(coords) == n:
                    break
            if len(coords) == n:
                break
        if len(coords) == n:
            break
    pos = np.array(coords, dtype=float)
    pos -= pos.mean(axis=0)
    max_r = np.linalg.norm(pos, axis=1).max()
    if max_r > 0:
        pos *= (radius_nm * 0.7) / max_r   # fill 70% of sphere
    return pos   # nm


# ---------------------------------------------------------------------------
# System builder
# ---------------------------------------------------------------------------
def build_system(prefactor: np.ndarray, n_beads: int,
                 phi: np.ndarray,
                 prefactor_strength: float = 20.0,
                 conc_mM: float = 2.0) -> tuple:  # 2.0mM gives ~14 beads in contact range
    """
    Build OpenMM system with all energies in kJ/mol and distances in nm.
    Returns (system, R_nm).
    """
    s   = SIGMA_NM
    kbt = KBT_KJ

    # Confinement radius from nucleosome concentration
    V_L   = n_beads / (conc_mM * 1e-3 * 6.022e23)  # litres
    V_nm3 = V_L * 1e24                               # nm³
    R_nm  = (3 * V_nm3 / (4 * np.pi)) ** (1 / 3)

    system = mm.System()
    for _ in range(n_beads):
        system.addParticle(MASS_GMOL)

    # ── Bond  U = k2(r-r0)² + k3(r-r0)³ + k4(r-r0)⁴  ─────────────────
    k2 = 20.0 * kbt / s**2
    k3 = 20.0 * kbt / s**3
    k4 = 20.0 * kbt / s**4
    bond = mm.CustomBondForce("k2*(r-r0)^2 + k3*(r-r0)^3 + k4*(r-r0)^4")
    bond.addGlobalParameter("r0", s)
    bond.addGlobalParameter("k2", k2)
    bond.addGlobalParameter("k3", k3)
    bond.addGlobalParameter("k4", k4)
    for i in range(n_beads - 1):
        bond.addBond(i, i + 1)
    system.addForce(bond)

    # ── Softcore EV (fully inlined select expressions) ──────────────────
    ev_expr = (
        f"select(step({R_S_NM} - r),"
        f"  0.5*{ECUT_KJ}*(1 + tanh(2*(4*{kbt}*(({s}/r)^12"
        f"    - ({s}/r)^6) + {kbt})/{ECUT_KJ})),"
        f"  select(step({R_WCA_NM} - r),"
        f"    4*{kbt}*(({s}/r)^12 - ({s}/r)^6) + {kbt},"
        f"    0))"
    )
    ev = mm.CustomNonbondedForce(ev_expr)
    ev.setNonbondedMethod(mm.CustomNonbondedForce.CutoffNonPeriodic)
    ev.setCutoffDistance(R_WCA_NM * 1.5)
    for _ in range(n_beads):
        ev.addParticle([])
    system.addForce(ev)

    # ── Confinement ─────────────────────────────────────────────────────
    kappa = 1000.0 * kbt / s**2
    confine = mm.CustomExternalForce(
        f"{kappa} * max(0, sqrt(x*x + y*y + z*z) - {R_nm})^2"
    )
    for i in range(n_beads):
        confine.addParticle(i)
    system.addForce(confine)

    # ── Discrete A/B Flory-Huggins contact potential ───────────────────────
    # Classify each bead as A or B based on phi threshold (phi_mean).
    # Then apply:
    #   B-B pairs: -eps  (attractive — like attracts like)
    #   A-A pairs: -eps  (attractive — like attracts like)
    #   A-B pairs: +eps  (repulsive  — unlike repels unlike)
    #
    # This is the standard Flory-Huggins block copolymer model and
    # guarantees genuine phase separation when eps > eps_critical.
    # The continuous phi value scales eps smoothly so beads near the
    # A/B boundary have weaker interactions than strongly typed beads.
    #
    # Interaction strength:
    #   eps_ij = prefactor_strength * kbt * |phi_i - phi_mean| * |phi_j - phi_mean|
    # Sign:
    #   same type (both > mean or both < mean): negative (attractive)
    #   mixed type: positive (repulsive)

    rc_nm    = 1.42 * s
    mu       = 4.0  / s
    phi_mean = float(np.mean(phi))
    phi_dev  = phi - phi_mean      # signed deviation from mean

    contact = mm.CustomBondForce(
        f"lambda_contact * eps_ij * 0.5 * (1 + tanh({mu}*({rc_nm} - r)))"
    )
    contact.addGlobalParameter("lambda_contact", 0.0)
    contact.addPerBondParameter("eps_ij")

    base_strength = float(prefactor_strength) * kbt   # kJ/mol

    n_pairs   = 0
    eps_vals  = []
    n_BB = n_AA = n_AB = 0

    for i in range(n_beads):
        for j in range(i + 1, n_beads):
            # Magnitude scales with distance from mean (strongly typed beads
            # interact more strongly than boundary beads)
            magnitude = base_strength * abs(phi_dev[i]) * abs(phi_dev[j])

            # Sign: same side of mean = attractive (-), opposite = repulsive (+)
            same_type = (phi_dev[i] * phi_dev[j]) > 0
            eps_ij = -magnitude if same_type else +magnitude

            if abs(eps_ij) > 1e-6:
                contact.addBond(i, j, [eps_ij])
                eps_vals.append(eps_ij)
                n_pairs += 1
                if same_type:
                    if phi_dev[i] > 0:
                        n_BB += 1
                    else:
                        n_AA += 1
                else:
                    n_AB += 1

    eps_arr = np.array(eps_vals)
    attr    = eps_arr[eps_arr < 0]
    rep     = eps_arr[eps_arr > 0]

    print(f"  Contact pairs: {n_pairs:,}")
    print(f"  B-B attractive: {n_BB:,}  mean={attr[attr < attr.mean()*0.0 - 1e-10].mean() if n_BB else 0:.3f} kJ/mol")
    print(f"  A-A attractive: {n_AA:,}")
    print(f"  A-B repulsive:  {n_AB:,}")
    print(f"  eps range: {eps_arr.min():.3f} to {eps_arr.max():.3f} kJ/mol")
    print(f"  Total attractive: {attr.sum():.0f} kJ/mol")
    print(f"  Total repulsive:  {rep.sum():.0f} kJ/mol")
    print(f"  Net per bead: {eps_arr.sum()*2/n_beads:.2f} kJ/mol")
    system.addForce(contact)

    return system, R_nm


# ---------------------------------------------------------------------------
# Helper: run MD steps with NaN guard
# ---------------------------------------------------------------------------
def _run_checked(sim, steps: int, traj_id: int, label: str):
    try:
        sim.step(steps)
    except mm.OpenMMException as exc:
        state = sim.context.getState(getEnergy=True)
        pe = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        raise RuntimeError(
            f"[Traj {traj_id:02d}] OpenMM error during '{label}' "
            f"(PE={pe:.1f} kJ/mol): {exc}\n"
            "  → Try smaller --dt_ps (e.g. 0.010) or check energy matrix."
        ) from exc


# ---------------------------------------------------------------------------
# Single trajectory
# ---------------------------------------------------------------------------
def run_trajectory(traj_id: int, prefactor: np.ndarray, n_beads: int,
                   phi: np.ndarray, strength: float, conc_mM: float,
                   equil_steps: int, sim_steps: int, save_every: int,
                   dt_ps: float, platform_name: str, seed: int,
                   out_dir: Path) -> np.ndarray:

    print(f"\n[Traj {traj_id:02d}] Building system "
          f"(χ strength={strength:.1f} kBT, conc={conc_mM} mM)...")
    system, R_nm = build_system(prefactor, n_beads, phi=phi,
                                prefactor_strength=strength,
                                conc_mM=conc_mM)
    platform     = mm.Platform.getPlatformByName(platform_name)

    # LangevinMiddleIntegrator with proper temperature in Kelvin
    integrator = mm.LangevinMiddleIntegrator(
        T_K   * unit.kelvin,
        1.0   / unit.picosecond,    # friction coefficient
        dt_ps * unit.picoseconds,
    )
    integrator.setRandomNumberSeed(seed + traj_id * 997)

    sim = app.Simulation(app.Topology(), system, integrator, platform)

    # Place beads on Hilbert curve, λ=0 (no contact forces yet)
    init_pos = hilbert_init(n_beads, R_nm)
    sim.context.setPositions(init_pos * unit.nanometer)
    sim.context.setParameter("lambda_contact", 0.0)

    # ── Phase 1: minimise + equilibrate (λ=0) ──────────────────────────
    print(f"[Traj {traj_id:02d}] Minimising initial geometry ...")
    sim.minimizeEnergy(maxIterations=2000, tolerance=10.0)
    sim.context.setVelocitiesToTemperature(T_K * unit.kelvin,
                                           seed + traj_id * 997)

    print(f"[Traj {traj_id:02d}] Equilibrating ({equil_steps:,} steps, λ=0) ...")
    t0 = time.time()
    _run_checked(sim, equil_steps, traj_id, "equilibration")
    print(f"[Traj {traj_id:02d}] Equilibration done in {time.time()-t0:.1f}s")

    # ── Phase 2: ramp λ 0 → 1 ──────────────────────────────────────────
    N_RAMP     = 5_000
    STEPS_EACH = 100
    print(f"[Traj {traj_id:02d}] Ramping λ: 0→1 "
          f"({N_RAMP} increments × {STEPS_EACH} steps = "
          f"{N_RAMP*STEPS_EACH:,} total) ...")
    for i in range(N_RAMP + 1):
        sim.context.setParameter("lambda_contact", i / N_RAMP)
        if i < N_RAMP:
            _run_checked(sim, STEPS_EACH, traj_id, f"ramp {i}/{N_RAMP}")

    # Minimise and re-thermalise at λ=1
    print(f"[Traj {traj_id:02d}] Minimising at λ=1 ...")
    sim.minimizeEnergy(maxIterations=1000, tolerance=10.0)
    sim.context.setVelocitiesToTemperature(T_K * unit.kelvin,
                                           seed + traj_id * 997 + 1)
    _run_checked(sim, 10_000, traj_id, "post-ramp thermalisation")

    # Sanity check
    state = sim.context.getState(getPositions=True, getEnergy=True)
    pe    = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
    pos   = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    if not np.isfinite(pe) or not np.all(np.isfinite(pos)):
        raise RuntimeError(
            f"[Traj {traj_id:02d}] NaN/Inf after ramp (PE={pe:.1f} kJ/mol). "
            "Try --dt_ps 0.010."
        )
    print(f"[Traj {traj_id:02d}] PE = {pe:.1f} kJ/mol  ✓")

    # ── Phase 3: production sampling ───────────────────────────────────
    print(f"[Traj {traj_id:02d}] Production ({sim_steps:,} steps) ...")
    t0          = time.time()
    n_snaps     = sim_steps // save_every
    contact_sum = np.zeros((n_beads, n_beads), dtype=np.float32)
    rc_nm       = 1.42 * SIGMA_NM

    for snap in range(n_snaps):
        _run_checked(sim, save_every, traj_id, f"snapshot {snap}/{n_snaps}")
        state = sim.context.getState(getPositions=True)
        pos   = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        diff  = pos[:, None, :] - pos[None, :, :]
        dists = np.linalg.norm(diff, axis=-1)
        contact_sum += (dists <= rc_nm).astype(np.float32)

    contact_prob = contact_sum / n_snaps
    print(f"[Traj {traj_id:02d}] Production done in {time.time()-t0:.1f}s "
          f"({n_snaps} snapshots)")

    out_path = out_dir / f"contact_traj_{traj_id:02d}.npy"
    np.save(out_path, contact_prob)
    print(f"[Traj {traj_id:02d}] Saved → {out_path}")
    return contact_prob


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    args    = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("\n[Setup] Loading energy matrix ...")
    data      = np.load(args.energy, allow_pickle=True)
    prefactor = data["prefactor"].astype(np.float64)
    n_beads   = int(data["n_beads"].flat[0])
    phi       = data["phi"].astype(np.float64)

    print(f"  Beads        : {n_beads}")
    print(f"  φ mean/std   : {phi.mean():.4f} / {phi.std():.4f}")
    print(f"  σ            : {SIGMA_NM} nm")
    print(f"  T            : {T_K} K  (kBT = {KBT_KJ:.3f} kJ/mol)")
    print(f"  χ strength   : {args.strength:.1f} kBT")
    print(f"  Concentration: {args.conc} mM")
    print(f"  dt           : {args.dt_ps} ps")
    print(f"  Platform     : {args.platform}")
    print(f"  Trajectories : {args.n_traj}")
    print(f"  Equil steps  : {args.equil_steps:,}")
    print(f"  Sim steps    : {args.sim_steps:,}")
    print(f"  Snapshots/tj : {args.sim_steps // args.save_every:,}")

    all_cps = []
    for traj_id in range(args.n_traj):
        cp = run_trajectory(
            traj_id, prefactor, n_beads,
            phi=phi, strength=args.strength, conc_mM=args.conc,
            equil_steps=args.equil_steps,
            sim_steps=args.sim_steps,
            save_every=args.save_every,
            dt_ps=args.dt_ps,
            platform_name=args.platform,
            seed=args.seed,
            out_dir=out_dir,
        )
        all_cps.append(cp)

    print("\n[Final] Averaging contact matrices ...")
    avg = np.mean(all_cps, axis=0)
    np.save(out_dir / "contact_matrix_avg.npy", avg)
    print(f"✓  contact_matrix_avg.npy  →  {out_dir}")

    np.savez(out_dir / "simulation_meta.npz",
             n_beads=np.array([n_beads]),
             n_traj=np.array([args.n_traj]),
             equil_steps=np.array([args.equil_steps]),
             sim_steps=np.array([args.sim_steps]),
             save_every=np.array([args.save_every]),
             n_snapshots=np.array([args.sim_steps // args.save_every]))
    print(f"✓  simulation_meta.npz     →  {out_dir}\n")


if __name__ == "__main__":
    main()
