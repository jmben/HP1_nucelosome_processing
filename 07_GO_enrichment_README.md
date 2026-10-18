Markdown · 07_GO_enrichment_README.md

# 07 — GO and KEGG Enrichment Analysis

This section covers gene ontology (GO) and KEGG pathway enrichment analysis performed using clusterProfiler. Enrichments were run on three categories of gene sets: (1) genes associated with HP1 ChIP-seq and CUT&RUN consensus peaks (unique and overlapping, within and across conditions), and (2) genes assigned to k-means sensitivity clusters from the nucleosome heatmaps (Section 04).

All analyses used the human annotation database `org.Hs.eg.db` with all annotated genes as the background set.

---

## Gene set inputs

### ChIP-seq and CUT&RUN peak-associated genes

Gene names were assigned to consensus peaks via TSS ±1 kb bedtools intersect (as described in Section 05). Enrichments were run on:

- Genes unique to each condition or mark
- Genes shared (overlapping) between conditions or marks

This was done for all ChIP-seq targets (HP1α, HP1β, H2AZ) and CUT&RUN marks across all conditions, both within a single condition (e.g., HP1α-unique vs HP1β-unique vs HP1α∩HP1β in WT) and across conditions (e.g., WT-unique vs KD-unique vs shared).

### Sensitivity cluster-associated genes

Genes assigned to each k-means cluster (k=7) were extracted from the cluster BED files output by `plotHeatmap` (Section 04). Enrichments were run separately per cluster for each of the four k-means sorting schemes:

| k-means file | Sorted by |
| --- | --- |
| `*_sort_k7.bed` | Nucleosome occupancy |
| `*_sensitivity_k7.bed` | MNase sensitivity |
| `*_H2AZ_sensitivity_k7.bed` | H2AZ sensitivity |
| `*_shH2AZ_sensitivity_k7.bed` | shH2AZ sensitivity |

---

## GO and KEGG enrichment (clusterProfiler)

Enrichment was run for all three GO ontologies (BP, MF, CC) and KEGG pathways. The example below shows the pattern used for each gene set; replace the input gene list and output prefix as appropriate.

```r
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

# Load gene list (gene symbols)
genes <- read.delim("gene_list.txt", header = FALSE)$V1

# Convert gene symbols to Entrez IDs
entrez_ids <- bitr(genes,
                   fromType = "SYMBOL",
                   toType   = "ENTREZID",
                   OrgDb    = org.Hs.eg.db)

# GO enrichment — Biological Process
go_bp <- enrichGO(gene          = entrez_ids$ENTREZID,
                  OrgDb         = org.Hs.eg.db,
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  readable      = TRUE)

# GO enrichment — Molecular Function
go_mf <- enrichGO(gene          = entrez_ids$ENTREZID,
                  OrgDb         = org.Hs.eg.db,
                  ont           = "MF",
                  pAdjustMethod = "BH",
                  readable      = TRUE)

# GO enrichment — Cellular Component
go_cc <- enrichGO(gene          = entrez_ids$ENTREZID,
                  OrgDb         = org.Hs.eg.db,
                  ont           = "CC",
                  pAdjustMethod = "BH",
                  readable      = TRUE)

# KEGG enrichment
kegg <- enrichKEGG(gene         = entrez_ids$ENTREZID,
                   organism     = "hsa",
                   pAdjustMethod = "BH")

# Export results
write.csv(as.data.frame(go_bp), "output_prefix_GO_BP.csv", row.names = FALSE)
write.csv(as.data.frame(go_mf), "output_prefix_GO_MF.csv", row.names = FALSE)
write.csv(as.data.frame(go_cc), "output_prefix_GO_CC.csv", row.names = FALSE)
write.csv(as.data.frame(kegg),  "output_prefix_KEGG.csv",  row.names = FALSE)

# Dot plot visualization
dotplot(go_bp, showCategory = 20) +
    ggtitle("GO Biological Process") +
    theme_classic(base_size = 12)
ggsave("output_prefix_GO_BP_dotplot.png", width = 8, height = 10)
```

The same pattern was applied to GO MF, CC, and KEGG results. Background gene set: all annotated human genes (clusterProfiler default).

---

## Output files

| File | Description |
| --- | --- |
| `*_GO_BP.csv` | GO Biological Process enrichment results |
| `*_GO_MF.csv` | GO Molecular Function enrichment results |
| `*_GO_CC.csv` | GO Cellular Component enrichment results |
| `*_KEGG.csv` | KEGG pathway enrichment results |
| `*_GO_BP_dotplot.png` | Dot plot — top 20 enriched BP terms |
| `*_GO_MF_dotplot.png` | Dot plot — top 20 enriched MF terms |
| `*_GO_CC_dotplot.png` | Dot plot — top 20 enriched CC terms |
| `*_KEGG_dotplot.png` | Dot plot — top 20 enriched KEGG pathways |

Naming follows the pattern `<mark>_<condition>_<subset>_<ontology>.*`, e.g., `HP1a_WT_unique_GO_BP.csv` or `sensitivity_cluster3_GO_BP.csv`.
