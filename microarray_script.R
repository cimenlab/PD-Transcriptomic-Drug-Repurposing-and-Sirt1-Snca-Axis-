---
scriptname: "microarray_script.R"
title: "PD Transcriptomic Drug Repurposing and Sirt1 Snca Axis   "
author: "CoRD-OMICS - CimenLab - Fatma Hotaslier"
date: "21.07.2026"
license: "MIT License (see LICENSE file for details)"
description: "This script performs differential gene expression analysis on microarray
data using limma, comparing Parkinson's disease and control samples across
postmortem, iPSC, and cell-based models. It covers raw quality control,
outlier exclusion, quantile normalization, probe-to-gene annotation, PCA
and heatmap visualization, differential expression testing, volcano plots,
marker gene panels, GO/KEGG enrichment analysis, STRING protein-protein
interaction network construction, and Spearman correlation analysis
between SIRT1 and SNCA expression."
---

# Load required Bioconductor/CRAN packages for microarray analysis.
# ------------------------- 0) Packages -------------------------
library(GEOquery)
library(limma)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(dplyr)
library(tibble)
library(tidyr)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)
library(patchwork)
library(httr)
library(jsonlite)
library(tidygraph)
library(ggraph)
library(igraph)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

# DEFINE GLOBAL THEME (for visual consistency)
{
theme_supplement <- theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", size = 22, hjust = 0.5),
    axis.title = element_text(face = "bold", size = 20),
    axis.text = element_text(face = "bold", size = 16),
    legend.title = element_text(face = "bold", size = 18),
    legend.text = element_text(size = 14),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "grey90", size = 0.3),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 18),
    strip.background = element_blank()  
  )
}
# Load raw microarray data and define platform/sample metadata.
# ------------------------- 1) Veri & Platform -------------------------
gse <- getGEO("GSE20292", GSEMatrix = TRUE)
eset <- gse[[1]]
expr_raw <- exprs(eset)                    # PROBE x SAMPLE (raw)
pheno <- pData(eset)
platform_id <- annotation(eset)
gpl <- getGEO(platform_id, AnnotGPL = TRUE)
gpl_table <- Table(gpl)

# Group labels (e.g. from "Parkinson" / "Control" text)
char_col <- "characteristics_ch1"
group <- case_when(
  grepl("Parkinson", eset[[char_col]], ignore.case = TRUE) ~ "PD",
  grepl("Control",   eset[[char_col]], ignore.case = TRUE) ~ "CTRL",
  TRUE ~ NA_character_
)
stopifnot(!any(is.na(group)))             
# Standardize sample names (CTRL_1, PD_1 ...)
new_names <- paste0(group, "_", ave(seq_along(group), group, FUN = seq_along))
colnames(expr_raw) <- new_names

# Identify and exclude outlier samples before normalization.
# ------------------------- 2) Outlier selection (BEFORE normalization) -------------------------


# List of samples to exclude as outliers
keep_samples <- c()

# Take only these from the raw data
expr_raw_keep <- expr_raw[, keep_samples, drop = FALSE]

# Group information
group_keep <- factor(gsub("_\\d+$", "", keep_samples), levels = c("CTRL","PD"))
group<- group_keep


# Optional raw-data quality control (QC) checks prior to normalization.
# ------------------------- 3) Raw QC (optional, before normalization) -------------------------
expr_long <- reshape2::melt(expr_raw_keep)
colnames(expr_long) <- c("Probe","Sample","Intensity")
p_raw_hist <- ggplot(expr_long, aes(Intensity)) +
  geom_histogram(bins = 100) +
  facet_wrap(~ Sample, scales = "free") +
  theme_minimal() +
  ggtitle("Raw Intensities Histogram")
p_raw_hist 
ggsave("raw_histogram_GSE20292.png", p_raw_hist, width = 12, height = 8)

png("raw_boxplot_GSE20292.png", width = 1200, height = 600)
boxplot(expr_raw_keep, las = 2, main = "Raw Intensities Boxplot", outline = FALSE)
dev.off()

meds <- apply(expr_raw_keep, 2, median, na.rm = TRUE)
iqrs <- apply(expr_raw_keep, 2, IQR, na.rm = TRUE)
write.csv(data.frame(Sample = names(meds), Median = meds, IQR = iqrs),
          "raw_median_IQR_GSE20292.csv", row.names = FALSE)

# Quantile normalization at the probe level.
# ------------------------- 4) Quantile Normalization (PROBE level) -------------------------
expr_norm_probe <- normalizeBetweenArrays(expr_raw_keep, method = "quantile")

write.csv(expr_norm, "expr_quantile_normalized_raw_GSE20292.csv", row.names = TRUE)

# Annotate probes with gene symbols and collapse probe-level data to gene-level expression.
# ------------------------- 5) Gene annotation & reduction to GENE level -------------------------
# Flexibly find the symbol column on the platform
cand <- c("Gene Symbol","Gene symbol","GENE_SYMBOL","Symbol","gene_symbol")
symbol_col <- cand[cand %in% colnames(gpl_table)][1]
if (is.na(symbol_col)) stop("Gene symbol column not found in the GPL table.")

annot <- gpl_table[, c("ID", symbol_col)]
colnames(annot) <- c("PROBE_ID","SYMBOL")
annot <- annot %>% filter(!is.na(SYMBOL), SYMBOL != "")

# Map Probe -> Symbol
probe_df <- as.data.frame(expr_norm_probe) %>%
  rownames_to_column("PROBE_ID") %>%
  inner_join(annot, by = "PROBE_ID")

# Multiple probes map to the same gene -> average (or limma::avereps)
expr_gene <- probe_df %>%
  select(-PROBE_ID) %>%
  group_by(SYMBOL) %>%
  summarize(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop") %>%
  column_to_rownames("SYMBOL")

write.csv(expr_gene, "expr_quantile_normalized_GSE20292_supp.csv", row.names = TRUE)
write.xlsx(expr_gene, file = "expr_quantile_normalized_GSE20292_supp.xlsx", rowNames = TRUE)

df_vsd_GSE20292 <- expr_gene

# Medians before
apply(expr_raw_keep, 2, median, na.rm = TRUE)

# Medians after
apply(expr_norm_probe, 2, median, na.rm = TRUE)
summary(expr_gene)


# PCA: visualize sample clustering for the normalized microarray data.
#PCA FOR MICROARRAY 
{
# expr_gene: normalized expression matrix (genes x samples)
class(expr_gene)
str(expr_gene)
summary(expr_gene)               # General summary
any(is.na(expr_gene))            # TRUE means NA present
any(is.infinite(expr_gene))   
expr_gene <- expr_gene[!apply(expr_gene, 1, function(x) any(is.na(x) | is.infinite(x))), ]
expr_gene <- expr_gene[apply(expr_gene, 1, sd, na.rm = TRUE) > 0, ]

# TRUE means Inf present
pca <- prcomp(t(expr_gene), scale. = TRUE)
pcdf <- as.data.frame(pca$x[, 1:2])
pcdf$Sample <- rownames(pcdf)
pcdf$Group <- ifelse(grepl("^ROT50_4W", pcdf$Sample), "ROT50_4W", "Control")

vars <- pca$sdev^2 / sum(pca$sdev^2) * 100

pca_microarray <- ggplot(pcdf, aes(PC1, PC2, color = Group, shape = Group, label = Sample)) +
  geom_point(size = 6, alpha = 0.9) +
  geom_text(vjust = 1.6, hjust = 0.5, size = 4.5, color = "black", fontface = "bold") +
  scale_color_manual(values = c("Control" = "forestgreen", "ROT50_4W" = "deeppink")) +
  scale_shape_manual(values = c("Control" = 16, "ROT50_4W" = 17)) +
  labs(
    title = " ",
    x = paste0("PC1: ", round(vars[1], 1), "% variance"),
    y = paste0("PC2: ", round(vars[2], 1), "% variance"),
    color = "Group", shape = "Group"
  ) +
  theme_supplement
print(pca_microarray)
ggsave("GSE35642_ROT50_12.04/SUPP_PCA_GSE36321_MICROARRAY_LAST.tiff",
       pca_microarray, width = 10, height = 10, dpi = 600, bg = "white")

}
# Sample-distance heatmap for microarray samples.
#heatmap  FOR MICROARRay
{
# ============================================================

vars_gene <- apply(expr_gene, 1, var)
# preparing annotation_col (must match sample names)
annotation_col <- data.frame(Group = ifelse(grepl("^ROT", colnames(expr_gene)), "Rot", "Control"))
rownames(annotation_col) <- colnames(expr_gene)
col <- colorRampPalette(c("yellow", "orange", "red"))(50)

heatmap_plot <- pheatmap(
  expr_gene,                        # like vsd in RNA-seq
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  color = col,
  annotation_col = annotation_col,
  scale = "row",
  fontsize = 10,
  fontsize_row = 10,
  fontsize_col = 14,
  treeheight_row = 10,
  treeheight_col = 10,
  border_color = NA
)
print(heatmap)
ggsave("figures/SUPP_heatmap_GSE20292_MICROARRAY.tiff",
       heatmap_plot$gtable,
       device = "tiff", width = 10, height = 10, dpi = 600, bg = "white")

}
# Differential expression analysis with limma at the gene level.
#Limma (GENE level)
{
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

contrast_matrix <- makeContrasts(PDvsCTRL = PD - CTRL, levels = design)

fit <- lmFit(expr_gene, design) |>
  contrasts.fit(contrast_matrix) |>
  eBayes()

res_all <- topTable(fit, coef = "PDvsCTRL", number = Inf, adjust.method = "BH")
deseq_GSE20292 <- res_all
res_sig <- res_all %>% filter(adj.P.Val < 0.05 & abs(logFC) >= 0.6)

write.csv(res_all, "DE_all_genes_GSE20292.csv")
write.xlsx(res_all, file = "DE_all_genes_GSE20292_supp.xlsx", rowNames = TRUE)

write.xlsx(res_sig, file = "DE_significant_genes_GSE20292_supp.xlsx", rowNames = TRUE)
write.csv(res_sig, "DE_significant_genes_GSE20292.csv")
}
# HEATMAP microarray
{
  # 1) We take only the significantly changed genes
  expr_gene_sig <- expr_gene[rownames(res_sig), ]
  
  # 2) Variance calculation (kept as you had it, if you want)
  vars_gene <- apply(expr_gene_sig, 1, var)
  
  # 3) annotation_col (same as the one you wrote)
  annotation_col <- data.frame(
    Group = ifelse(grepl("^ROT50_4W", colnames(expr_gene_sig)), "ROT50_4W", "Control")
  )
  rownames(annotation_col) <- colnames(expr_gene_sig)
  
  # 4) color palette (your palette)
  col<- colorRampPalette(c("blue", "white", "red"))(50)
  
  # 5) Heatmap (your heatmap code)
  heatmap_plot <- pheatmap(
    expr_gene_sig,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = FALSE,
    color = col,
    annotation_col = annotation_col,
    scale = "row",
    fontsize = 10,
    fontsize_row = 10,
    fontsize_col = 14,
    treeheight_row = 10,
    treeheight_col = 10,
    border_color = NA
  )
  
  print(heatmap_plot)
  
  # 6) Save (your format)
  ggsave(
    "GSE35642_ROT50_12.04/SUPP_heatmap_GSE35642_ROT50_MICROARRAY_LAST.tiff",
    heatmap_plot$gtable,
    device = "tiff",
    width = 10,
    height = 10,
    dpi = 600,
    bg = "white"
  )
}
# Volcano plot: visualize significant differentially expressed genes.
#Volcano 
{
# ============================================================
  # VOLCANO PLOT (microarray)
# ============================================================
n_up   <- sum(res_all$Regulation == "Up")
n_down <- sum(res_all$Regulation == "Down")

volcano_standard <- EnhancedVolcano(
  res_all,
  lab = rownames(res_all),
  x = 'logFC',              # for RNA-seq write 'log2FoldChange'
  y = 'adj.P.Val',          # for RNA-seq write 'padj'
  pCutoff = 0.05,
  FCcutoff = 0.6,
  title = "PD vs Healthy",
  subtitle = "Differential Expression",
  caption = NULL,
  axisLabSize = 18,
  titleLabSize = 22,
  subtitleLabSize = 18,
  pointSize = 5,
  widthConnectors = 0.5,
  colConnectors = 'black',
  cutoffLineCol = "black",
  cutoffLineType = "dotdash",
  cutoffLineWidth = 1,
  hlineCol = "black",
  hlineType = "dotdash",
  hlineWidth = 0.3,
  labSize = 10,
  legendLabels = c('NS','|logFC|>0.6','p<0.05','p<0.05 & |logFC|>0.6'),
  legendPosition = 'bottom',
  legendIconSize = 6,
  drawConnectors = FALSE,
  legendLabSize = 14,
  gridlines.major = TRUE,
  gridlines.minor = FALSE,
  border = 'partial',
  borderWidth = 1,
  borderColour = 'black'
)

x_lim <- ggplot_build(volcano_standard)$layout$panel_params[[1]]$x.range
y_lim <- ggplot_build(volcano_standard)$layout$panel_params[[1]]$y.range

final_volcano <- volcano_standard +
  annotate("text",
           x = x_lim[2] * 0.85, y = y_lim[2] * 0.95,
           label = paste0("Up: ", n_up),
           size = 6, fontface = "bold", color = "black", hjust = 1) +
  annotate("text",
           x = x_lim[1] * 0.85, y = y_lim[2] * 0.95,
           label = paste0("Down: ", n_down),
           size = 6, fontface = "bold", color = "black", hjust = 0) +
  theme_supplement

print(final_volcano)

ggsave(
  filename = "GSE35642_ROT50_12.04/SUPP_VOLCANO_GSE35642_ROT50_microarray_LAST.tiff",
  plot = final_volcano,
  width = 10, height = 10, dpi = 600, bg = "white"
)
}
# GO and KEGG enrichment analysis on the set of differentially expressed genes.
#GO KEGG
{
  # ============================================================
  # Microarray GO Enrichment (BP, MF, CC)
  # ============================================================
  
  deg_genes <- as.character(rownames(res_sig))   # significant genes (adj.P.Val < 0.05, |logFC|>0.6)
  universe <- rownames(res_all)                  # background of all genes
  
  library(clusterProfiler)
  library(org.Hs.eg.db)
  
  
  
  gene_ids <- bitr(rownames(res_sig),
                   fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)
  
  universe_ids <- bitr(rownames(res_all),
                       fromType = "SYMBOL",
                       toType = "ENTREZID",
                       OrgDb = org.Hs.eg.db)
  
  
  egoBP <- enrichGO(
    gene          = gene_ids$ENTREZID,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = "BP",
    universe      = universe_ids$ENTREZID,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  
  egoMF <- enrichGO(
    gene          = gene_ids$ENTREZID,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = "MF",
    universe      = universe_ids$ENTREZID,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  egoCC <- enrichGO(
    gene          = gene_ids$ENTREZID,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = "CC",
    universe      = universe_ids$ENTREZID,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  
  df_bp <- egoBP@result |> mutate(ontology = "BP")
  df_mf <- egoMF@result |> mutate(ontology = "MF")
  df_cc <- egoCC@result |> mutate(ontology = "CC")
  
  combined_go <- bind_rows(df_bp, df_mf, df_cc) |>
    mutate(
      GeneRatio = sapply(strsplit(GeneRatio, "/"), \(x) as.numeric(x[1]) / as.numeric(x[2])),
      GeneRatio = as.numeric(GeneRatio)
    )
  
  top_terms <- combined_go |>
    group_by(ontology) |>
    arrange(p.adjust, desc(Count)) |> 
    slice_head(n = 10) |>
    ungroup()
  
  go_plot <- ggplot(top_terms, aes(
    x = GeneRatio,
    y = reorder(Description, GeneRatio),
    color = p.adjust,
    size = Count
  )) +
    geom_point() +
    scale_color_gradient(low = "red", high = "blue") +
    facet_wrap(~ ontology, scales = "free_y") +
    labs(
      title = "Top 10 Enriched GO Terms per Category ",
      x = "Gene Ratio",
      y = "GO Term Description",
      color = "Adjusted p-value",
      size = "Gene Count"
    ) +
    theme_supplement
  
  
  print(go_plot)
  
  ggsave("GSE35642_ROT50_12.04/SUPP_GO_GSE35642_ROT50_12.04_MICROARRAY_LAST.tiff",
         go_plot, width = 30, height = 15, dpi = 600, bg = "white")
  
  # ============================================================
  # 5. KEGG ENRICHMENT DOTPLOT
  # ============================================================
  gene_ids <- bitr(deg_genes,
                   fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)
  head(gene_ids)
  
  kegg_result <- enrichKEGG(
    gene         = gene_ids$ENTREZID,
    organism     = 'hsa',
    pvalueCutoff = 0.05
  )
  
  kegg_plot <- dotplot(kegg_result, showCategory = 10) +
    ggtitle("Top 10 Enriched KEGG Pathways") +
    theme_supplement
  
  ggsave("GSE35642_ROT50_12.04/SUPP_KEGG_GSE35642_ROT50_12.04_microarray_LAST.tiff",
         kegg_plot, width = 10, height = 10, dpi = 600, bg = "white")
  
  # ============================================================
  # PACKAGES
  # ============================================================
  library(ggplot2)
  library(patchwork)
  
  # ============================================================
  # GO + KEGG PLOT COMBINATION
  # ============================================================
  
  go_plot <- go_plot +
    scale_y_discrete(labels = function(x) str_wrap(x, 30))
  print(go_plot)
  kegg_plot <- kegg_plot +
    scale_y_discrete(labels = function(x) str_wrap(x, 10))
  
  combined_plot <- go_plot + kegg_plot + 
    plot_layout(ncol = 2) +              # side by side (nrow = 2 would stack them)
    plot_annotation(
      title = "GO and KEGG Enrichment Summary",
      theme = theme(
        plot.title = element_text(face = "bold", size = 24, hjust = 0.5)
      )
    )
  
  
  ggsave(
    filename = "GSE35642_ROT50_12.04/SUPP_GO_KEGG_COMBINED_GSE35642_ROT50_12.04.tiff",
    plot = combined_plot,
    device = "tiff",
    width = 32,        # two plots side by side, double width
    height = 10,
    dpi = 600,
    bg = "white"
  )
  # ---------------- KEGG ----------------
  write.csv(kegg_result,
            file = file.path(folder_path, "GSE35642_ROT50_12.04_SUPP_kegg_result_SUPP.csv"),
            row.names = TRUE)
  
  write.xlsx(kegg_result,
             file = file.path(folder_path, "GSE35642_ROT50_12.04_SUPP_kegg_result_SUPP.xlsx"),
             rowNames = TRUE)
  
  # ---------------- GO ----------------
  write.csv(combined_go,
            file = file.path(folder_path, "GSE35642_ROT50_12.04_GO_terms_SUPP.csv"),
            row.names = TRUE)
  
  write.xlsx(combined_go,
             file = file.path(folder_path, "GSE35642_ROT50_12.04_SUPP_GO_terms_SUPP.xlsx"),
             rowNames = TRUE)}
# Marker gene panel: plot expression of selected Parkinson's disease-related marker genes across groups.
#marker gene
{
  genes_pd_related <- c("SNCA", "PRKN", "PINK1", "LRRK2", "VPS35", "ATP13A2", "PLA2G6", "FBXO7", "GBA2")
  genes_dopa_marker <- c("TH", "DDC", "SLC18A2", "NR4A2", "ALDH1A1", "KCNJ6", "CALB1", "DRD2")
  genes_sirtuins <- c("SIRT1", "SIRT2", "SIRT3", "SIRT4", "SIRT5", "SIRT6", "SIRT7")
  genes_tet <- c("TET1", "TET2", "TET3")
  genes_lewy <- c("UBB", "HSP90AA1", "TGM2", "STUB1")
  
  
  get_group_data <- function(gene_set, group_name) {
    res_all %>%
      tibble::rownames_to_column("symbol") %>%       # move rownames into the symbol column
      dplyr::filter(symbol %in% gene_set) %>%
      dplyr::select(symbol, logFC, adj.P.Val) %>%
      dplyr::mutate(Group = group_name)
  }
  
  
  # Merge all groups
  df_markers <- bind_rows(
    get_group_data(genes_pd_related, "Parkinson"),
    get_group_data(genes_dopa_marker, "Dopamine"),
    get_group_data(genes_sirtuins, "SIRT"),
    get_group_data(genes_tet, "TET"),
    get_group_data(genes_lewy, "LewyBody")
  )
  colnames(df_markers)
  
  
  marker_plot_small <- ggplot(df_markers, aes(
    x = reorder(symbol, logFC),
    y = logFC,
    fill = adj.P.Val < 0.05 & abs(logFC) > 0.6
  )) +
    geom_bar(stat = "identity") +
    geom_text(
      aes(label = ifelse(!is.na(adj.P.Val),
                         paste0("p=", signif(adj.P.Val, 2)),
                         "")),
      vjust = -0.4,
      size = 3,
      fontface = "bold",
      color = "black"
    ) +
    facet_wrap(~Group, scales = "free_x", nrow = 1) +
    scale_fill_manual(values = c("TRUE" = "firebrick3", "FALSE" = "grey70")) +
    labs(
      title = "Marker Gene Groups: logFC and Significance (adj.P.Val < 0.05 & |log2FC|> 0.6)",
      x = "Gene",
      y = "log2 Fold Change",
      fill = "Significant"
    ) +
    theme_supplement +
    theme(
      legend.position = "top",
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1,
        size = 18,
        face = "bold"
      )
    )
  
  
  
  print(marker_plot_small)
  ggsave("GSE35642_ROT50_12.04/SUPP_MARKER_GROUPS_GSE35642_ROT50_12.04_microarray_LAST.tiff",
         marker_plot_small,
         width = 30, height = 10, dpi = 600, bg = "white")
  }
# STRING protein-protein interaction network: query STRING and build the target-gene subnetwork.
#string 
{
  # ---------------------- 1. REQUIRED LIBRARIES ----------------------
  library(httr)
  library(jsonlite)
  library(tidyverse)
  library(igraph)
  library(ggraph)
  library(tidygraph)
  library(ggplot2)
  
  # ---------------------- 2. FETCH PPI DATA FROM THE STRING API ----------------------
  
  # Significant genes (example filter: adj.P.Val < 0.05 & |logFC| > 0.6)
  sig_genes <- rownames(res_all %>% filter(adj.P.Val < 0.05 & abs(logFC) > 0.6))
  
  # API URL
  string_api_url <- "https://string-db.org/api/json/network"
  
  # Split genes into chunks of 50
  gene_batches <- split(sig_genes, ceiling(seq_along(sig_genes) / 50))
  all_string_data <- list()
  
  # Fetch data from the API for each batch
  for (i in seq_along(gene_batches)) {
    params <- list(
      identifiers = paste(gene_batches[[i]], collapse = "%0d"),
      species = 9606,
      required_score = 700
    )
    response <- GET(string_api_url, query = params)
    
    if (http_status(response)$category == "Success") {
      string_data <- fromJSON(content(response, as = "text", encoding = "UTF-8"))
      all_string_data[[i]] <- string_data
    } else {
      warning(paste("STRING API error! Batch", i))
    }
  }
  
  # Combine the data
  ppi_network <- bind_rows(all_string_data)
  
  # Save as CSV (optional)
  write.csv(ppi_network, "STRING_PPI_network_GSE20292.csv", row.names = FALSE)
  
  # ---------------------- 3. BUILD GRAPH (igraph + ggraph) ----------------------
  
  # Edge table
  edges <- ppi_network[, c("preferredName_A", "preferredName_B", "score")]
  colnames(edges) <- c("from", "to", "weight")
  
  # Build network (undirected graph)
  ppi_graph <- graph_from_data_frame(edges, directed = FALSE)
  
  # ---------------------- 4. BUILD VISUAL ----------------------
  
  set.seed(42)  # make the layout reproducible
  
  p <- ggraph(ppi_graph, layout = "fr") +
    geom_edge_link(aes(width = weight), edge_colour = "grey80", alpha = 0.5) +
    geom_node_point(color = "steelblue", size = 4) +
    geom_node_text(aes(label = name), repel = TRUE, size = 3.2) +
    theme_void() +
    ggtitle("STRING PPI Network - GSE20292") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  print(p)
  # Save the figure as PNG
  ggsave("STRING_PPI_network_plot_GSE20292_Rapor.png", plot = p, width = 10, height = 8, dpi = 300)
  # Get all node names
  all_genes_in_ppi <- unique(c(ppi_network$preferredName_A, ppi_network$preferredName_B))
  
  # Is SIRT1 present?
  "SIRT1" %in% all_genes_in_ppi
  "SIRT1" %in% sig_genes  # if this is also FALSE -> it was never sent!
  
  # Get all node names
  all_genes_in_ppi <- unique(c(ppi_network$preferredName_A, ppi_network$preferredName_B))
  all_genes_in_ppi <- sort(unique(c(ppi_network$preferredName_A, ppi_network$preferredName_B)))
  all_genes_in_ppi <- all_genes_in_ppi[!is.na(all_genes_in_ppi)]  # remove NA
  
  # Print as a single comma-separated line in the console
  cat(paste(all_genes_in_ppi, collapse = ","))
  
  
  
  # ---------------------- 6) Subnetwork for genes of interest ----------------------
  target_genes <- c("NPC1","FOXO3")
  
  filtered_data <- ppi_network %>%
    filter(preferredName_A %in% target_genes | preferredName_B %in% target_genes)
  
  # Create the figures folder if it does not exist
  if (!dir.exists("figures")) dir.create("figures", recursive = TRUE)
  
  write.csv(filtered_data, "figures/cytoscape_network_target_genes.csv", row.names = FALSE)
  
  # Subnetwork
  subgraph_edges <- filtered_data %>%
    select(preferredName_A, preferredName_B, score) %>%
    rename(from = preferredName_A, to = preferredName_B, score = score) %>%
    filter(!is.na(from), !is.na(to), !is.na(score))
  
  if (nrow(subgraph_edges) > 0) {
    
    if (max(subgraph_edges$score, na.rm = TRUE) > 1) subgraph_edges$score <- subgraph_edges$score / 1000
    
    subgraph_data <- graph_from_data_frame(subgraph_edges, directed = FALSE)
    V(subgraph_data)$degree <- degree(subgraph_data)
    
    # Safe fallback in the color scale to avoid errors
    deg <- V(subgraph_data)$degree
    if (length(unique(deg[is.finite(deg)])) <= 1) {
      color_scale <- scale_color_gradientn(
        colours = c("blue", "cyan", "yellow", "orange", "red"),
        values  = c(0, 1)
      )
    } else {
      color_scale <- scale_color_gradientn(
        colours = c("blue", "cyan", "yellow", "orange", "red"),
        values = scales::rescale(c(
          min(deg, na.rm = TRUE),
          median(deg, na.rm = TRUE),
          max(deg, na.rm = TRUE)
        ))
      )
    }
    
    # Simple subnetwork visual
    p_sub_simple <- ggraph(subgraph_data, layout = "fr") +
      geom_edge_link(aes(edge_alpha = score), color = "grey40") +
      scale_edge_alpha_continuous(range = c(0.3, 1), name = "Interaction Score") +
      geom_node_point(size = 5, color = "steelblue") +
      geom_node_text(aes(label = name), repel = TRUE, size = 4) +
      theme_void() +
      labs(
        title = "Filtered STRING Network",
        subtitle = "Selected genes only"
      )
    print(p_sub_simple)
    ggsave("figures/GSE20291_RESVERATROL_RAPOR.tiff", plot = p_sub_simple,
           width = 12, height = 10, dpi = 600, device = "tiff", bg = "white")
    
    # Subnetwork colored/sized by degree
    p_sub_deg <- ggraph(subgraph_data, layout = "fr") +
      geom_edge_link(aes(edge_alpha = score), color = "grey40", show.legend = TRUE) +
      scale_edge_alpha_continuous(range = c(0.3, 1), name = "Interaction Score") +
      geom_node_point(aes(size = degree, color = degree), alpha = 1, show.legend = TRUE) +
      color_scale +
      scale_size_continuous(range = c(3, 8), name = "Node Degree") +
      geom_node_text(aes(label = name), size = 3, repel = TRUE) +
      guides(color = guide_colorbar(title = "Node Degree"),
             size  = guide_legend(title = "Node Degree")) +
      theme_void() +
      labs(
        title = "",
        subtitle = ""
      )
    print(p_sub_deg)
    
    ggsave("figures/GSE20292_String_resveratrol_renkli.tiff", plot = p_sub_deg,
           width = 12, height = 10, dpi = 600, device = "tiff", bg = "white")
    
  } else {
    message()
  }
}
# Spearman correlation between SIRT1 and SNCA expression across samples.
#correlation
{
  
  df <- df_vsd_GSE20292 
  
  gene1 <- "SIRT1"
  gene2 <- "SNCA"
  
  # 1) Select the two relevant genes as rows
  gene_pair <- c(gene1, gene2)
  
  if ("symbol" %in% colnames(df)) {
    df_sub <- df[df$symbol %in% gene_pair, ]
    rownames(df_sub) <- df_sub$symbol
    df_sub <- df_sub[, colnames(df_sub) != "symbol", drop = FALSE]
  } else {
    df_sub <- df[rownames(df) %in% gene_pair, ]
  }
  
  # 2) Make sure both genes are present
  if (!all(gene_pair %in% rownames(df_sub))) {
    stop()
  }
  
  # 3) Convert to numeric vectors (across samples)
  x <- as.numeric(df_sub[gene1, ])   # SIRT1 expression values
  y <- as.numeric(df_sub[gene2, ])   # SNCA expression values
  
  # 4) Check how many samples (n) were used
  n_used <- sum(!is.na(x) & !is.na(y))
  cat("Kullan??lan ??rnek say??s?? (n):", n_used, "\n")
  
  # 5) Spearman correlation test
  test <- cor.test(x, y, method = "spearman", exact = FALSE)
  
  # 6) Print the results properly
  rho  <- round(as.numeric(test$estimate), 3)
  pval <- round(test$p.value, 4)
  
  cat("Spearman rho :", rho, "\n")
  cat("p-de??eri     :", pval, "\n")
  cat("S istatisti??i:", test$statistic, "\n")
  
  # 7) Optionally, inspect the full test object as well
  print(test)
}



