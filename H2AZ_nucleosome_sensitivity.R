#!/usr/bin/env Rscript
# Nucleosome sensitivity analysis across conditions
#
# Input: TSS_10A_*_sensitivity_overview.txt files
# - Deduplicates to unique nucleosomes (cols 1-4: chr, start, end, category)
# - Counts Sensitive / Resistant / Neutral per sample
# - Outputs summary table, per-sample stacked bar plots, and
#   combined condition plot with error bars
#
# Usage:
#   Rscript nucleosome_sensitivity.R

library(dplyr)
library(tidyr)
library(ggplot2)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
DATA_DIR  <- "/data2/jbenoit/reduced_condense_seq"
OUT_TABLE  <- file.path(DATA_DIR, "nucleosome_sensitivity_summary.txt")
OUT_PLOT1  <- file.path(DATA_DIR, "nucleosome_sensitivity_counts.pdf")
OUT_PLOT2  <- file.path(DATA_DIR, "nucleosome_sensitivity_proportions.pdf")
OUT_PLOT3  <- file.path(DATA_DIR, "nucleosome_sensitivity_combined.pdf")

# Category factor order and colors
CAT_LEVELS <- c("Sensitive", "Neutral", "Resistant")
CAT_COLORS <- c(Sensitive = "#D6604D", Neutral = "#BABABA", Resistant = "#4393C3")

# Sample order and condition grouping for plot x-axis
SAMPLE_META <- data.frame(
  file = c(
    "TSS_10A_Scr_clone1_sensitivity_overview.txt",
    "TSS_10A_Scr_clone2_sensitivity_overview.txt",
    "TSS_A_shH2AZ_rep1_sensitivity_overview.txt",
    "TSS_A_shH2AZ_rep2_sensitivity_overview.txt",
    "TSS_10A_HP1b_clone15_sensitivity_overview.txt",
    "TSS_10A_HP1b_cloneD_sensitivity_overview.txt"
  ),
  sample = c(
    "Scr_clone1",   "Scr_clone2",
    "H2AZ_clone1A",  "H2AZ_clone4",
    "HP1b_clone15", "HP1b_cloneD"
  ),
  condition = c("Scr", "Scr", "H2AZ_KD", "H2AZ_KD", "HP1b_KD", "HP1b_KD"),
  stringsAsFactors = FALSE
)

# ---------------------------------------------------------------------------
# Helper: read and deduplicate one file
# ---------------------------------------------------------------------------
read_sensitivity <- function(filepath, sample_name) {
  df <- read.table(filepath, header = FALSE, sep = "\t",
                   col.names = c("chr", "start", "end", "category",
                                 "gene_chr", "gene_start", "gene_end",
                                 "gene_name", "score", "strand"))
  df %>%
    distinct(chr, start, end, category) %>%
    mutate(sample   = sample_name,
           category = factor(category, levels = CAT_LEVELS))
}

# ---------------------------------------------------------------------------
# Step 1: Read all files
# ---------------------------------------------------------------------------
cat("Reading files...\n")
all_data <- mapply(
  read_sensitivity,
  filepath    = file.path(DATA_DIR, SAMPLE_META$file),
  sample_name = SAMPLE_META$sample,
  SIMPLIFY    = FALSE
) %>% bind_rows()

# ---------------------------------------------------------------------------
# Step 2: Per-sample counts and proportions
# ---------------------------------------------------------------------------
summary_df <- all_data %>%
  group_by(sample, category) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(sample) %>%
  mutate(
    total      = sum(count),
    proportion = count / total
  ) %>%
  ungroup() %>%
  left_join(SAMPLE_META[, c("sample", "condition")], by = "sample") %>%
  mutate(
    sample    = factor(sample,    levels = SAMPLE_META$sample),
    condition = factor(condition, levels = c("Scr", "H2AZ_KD", "HP1b_KD")),
    category  = factor(category,  levels = CAT_LEVELS)
  )

# ---------------------------------------------------------------------------
# Step 3: Combined condition summary (mean +/- SD of proportions)
# ---------------------------------------------------------------------------
combined_df <- summary_df %>%
  group_by(condition, category) %>%
  summarise(
    mean_prop = mean(proportion),
    sd_prop   = sd(proportion),
    n         = n(),
    .groups   = "drop"
  ) %>%
  mutate(
    se_prop   = sd_prop / sqrt(n),
    condition = factor(condition, levels = c("Scr", "H2AZ_KD", "HP1b_KD")),
    category  = factor(category,  levels = CAT_LEVELS)
  )

# ---------------------------------------------------------------------------
# Step 4: Write summary table
# ---------------------------------------------------------------------------
summary_wide <- summary_df %>%
  select(condition, sample, category, count, proportion) %>%
  arrange(condition, sample, category)

write.table(summary_wide, OUT_TABLE, sep = "\t", quote = FALSE,
            row.names = FALSE)
cat(sprintf("Saved summary table     -> %s\n", OUT_TABLE))

cat("\n--- Nucleosome counts per sample ---\n")
print(as.data.frame(summary_wide), row.names = FALSE)

cat("\n--- Combined condition summary (mean +/- SD) ---\n")
print(as.data.frame(combined_df), row.names = FALSE)

# ---------------------------------------------------------------------------
# Step 5: Plot 1 — raw counts stacked bar
# ---------------------------------------------------------------------------
p1 <- ggplot(summary_df,
             aes(x = sample, y = count, fill = category)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = count),
            position = position_stack(vjust = 0.5),
            size = 3, color = "white", fontface = "bold") +
  scale_fill_manual(values = CAT_COLORS, name = "Category") +
  facet_grid(~ condition, scales = "free_x", space = "free_x") +
  labs(
    title = "Nucleosome sensitivity at TSS — raw counts",
    x     = NULL,
    y     = "Number of nucleosomes"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "grey90"),
    strip.text       = element_text(face = "bold"),
    legend.position  = "right"
  )

pdf(OUT_PLOT1, width = 9, height = 6)
print(p1)
dev.off()
cat(sprintf("Saved counts plot       -> %s\n", OUT_PLOT1))

# ---------------------------------------------------------------------------
# Step 6: Plot 2 — per-sample proportions stacked bar
# ---------------------------------------------------------------------------
p2 <- ggplot(summary_df,
             aes(x = sample, y = proportion, fill = category)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", proportion * 100)),
            position = position_stack(vjust = 0.5),
            size = 3, color = "white", fontface = "bold") +
  scale_fill_manual(values = CAT_COLORS, name = "Category") +
  scale_y_continuous(labels = scales::percent_format()) +
  facet_grid(~ condition, scales = "free_x", space = "free_x") +
  labs(
    title = "Nucleosome sensitivity at TSS — proportions",
    x     = NULL,
    y     = "Proportion of nucleosomes"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "grey90"),
    strip.text       = element_text(face = "bold"),
    legend.position  = "right"
  )

pdf(OUT_PLOT2, width = 9, height = 6)
print(p2)
dev.off()
cat(sprintf("Saved proportions plot  -> %s\n", OUT_PLOT2))

# ---------------------------------------------------------------------------
# Step 7: Plot 3 — combined conditions, grouped bars with error bars
# ---------------------------------------------------------------------------

# Updated color palette
CAT_COLORS <- c(Sensitive = "#D4A017", Neutral = "#BABABA", Resistant = "#1B3A6B")

p3 <- ggplot(combined_df,
             aes(x = condition, y = mean_prop, fill = category)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8),
           width = 0.7, color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = mean_prop - sd_prop,
                    ymax = mean_prop + sd_prop),
                position = position_dodge(width = 0.8),
                width = 0.3, linewidth = 0.8, color = "black") +
  geom_point(data = summary_df,
             aes(x = condition, y = proportion,
                 fill = category, group = category),
             position = position_dodge(width = 0.8),
             shape = 21, size = 2.5, color = "black", stroke = 0.6,
             show.legend = FALSE) +
  scale_fill_manual(values = CAT_COLORS, name = NULL) +
  scale_x_discrete(labels = c("Scr" = "Scramble",
                               "H2AZ_KD" = "HP1\u03b1 KD",
                               "HP1b_KD" = "HP1\u03b2 KD")) +
  scale_y_continuous(labels = scales::percent_format(),
                     limits = c(0, 0.85),
                     expand = expansion(mult = c(0, 0))) +
  labs(
    x = NULL,
    y = "Proportion of nucleosomes"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text        = element_text(color = "black", size = 12),
    axis.title.y     = element_text(face = "bold", size = 13),
    axis.line        = element_line(color = "black", linewidth = 0.5),
    axis.ticks       = element_line(color = "black"),
    legend.position  = "right",
    legend.text      = element_text(size = 11),
    legend.key.size  = unit(0.5, "cm"),
    plot.margin      = margin(10, 10, 10, 10)
  )

pdf(OUT_PLOT3, width = 7, height = 5, useDingbats = FALSE)
print(p3)
dev.off()
cat(sprintf("Saved combined plot     -> %s\n", OUT_PLOT3))

cat("\nDone.\n")
