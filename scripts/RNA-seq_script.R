---
scriptname: "RNA-seq_script.R"
title: "PD Transcriptomic Drug Repurposing and Sirt1 Snca Axis   "
author: "CoRD-OMICS - CimenLab - Fatma Hotaslier"
date: "21.07.2026"
license: "MIT License (see LICENSE file for details)"
description: "This script performs differential gene expression analysis on RNA-seq count
data using DESeq2, comparing Parkinson's disease and control samples across
postmortem, iPSC, and cell-based models. It covers pre-filtering, variance-
stabilizing transformation, PCA and heatmap visualization, differential
expression testing, volcano plots, marker gene panels, GO/KEGG enrichment
analysis, STRING protein-protein interaction network construction, and
Spearman correlation analysis between SIRT1 and SNCA expression."
---

# Load required Bioconductor/CRAN packages for RNA-seq differential expression analysis.
#packages 
{
  library(GenomeInfoDb)
  library(gridExtra)
  library(tidyr)
  library(BiocGenerics)
  library(magrittr)
  library(stringr)
  library(ggplot2)
  library(org.Hs.eg.db)
  library(DOSE)
  library(DESeq2)
  library(clusterProfiler)
  library(RSQLite)
  library(enrichplot)
  library(msigdbr)
  library(purrr)
  library(stringr)
  library(curl)
  library(vsn)
  library("RColorBrewer") 
  library("pheatmap") 
  library("PoiClaClu")
  library("glmpca")
  library("ggbeeswarm")
  library("apeglm")
  library("genefilter")
  library(reshape)
  library(EnhancedVolcano)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(GenomicFeatures)
  library("tximport")
  library("tximportData") 
  library(DESeq2)
  library(ggplot2)
  library(plotly)
  install.packages("BiocManager") 
  BiocManager::install("DOSE")
  install.packages("pheatmap")
  install.packages("BiocManager")
  BiocManager::install("clusterProfiler")
  library(clusterProfiler)
  library(org.Hs.eg.db) # For human (as an example)
  # Install biomaRt via Bioconductor
  BiocManager::install("biomaRt")
  library(biomaRt)
  library(dplyr)
  library(ggplot2)
  library(dplyr)  
  library(tidyr)            
  library(readxl)
  library(dplyr)
  
  library(ggplot2)
  library(org.Hs.eg.db)  
  library(enrichplot)
  library(pheatmap)
  library(ggraph)
}

setwd("C://Users//User//OneDrive//Belgeler")

cts <- read.csv("C://Users//User//OneDrive//Belgeler//SNc_count.csv", 
                sep = ";", 
                header = TRUE,
                row.names= 1)      

cts<-as.matrix(cts) #convert to matrix
coldata<- read.csv("coldata_GSE114918.csv", header=TRUE, sep=";") 

# 2. Assign the Sample column to row names
rownames(coldata) <- coldata$Sample

# Define sample names to select the relevant samples
keep_samples <- c()

# Filter coldata
coldata <- coldata[coldata$Sample %in% keep_samples, ]

# Make sure the column names in cts match coldata$sample
cts <- cts[, colnames(cts) %in% coldata$Sample]


# Convert the Group column to a factor
coldata$Group <- trimws(coldata$Group) 
coldata$Group <- as.factor(coldata$Group)
levels(coldata$Group)
coldata$Group <- relevel(coldata$Group, ref = "Healthy")

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~Group)
dds



# Pre-filtering: remove genes with very low counts across samples to reduce noise before analysis.
#Pre-filtering
{
  keep <- rowSums(counts(dds)) >= 10  #Calculates each gene's total expression across all samples.
  dds <- dds[keep,] 
  dds
  rm(keep) #remove keep 
  dim(dds) #count of gene 
}
# Variance-stabilizing transformation (VST): normalize counts for visualization, clustering, and PCA.
#The variance stabilizing transformation and the rlog
{
  vsd <- vst(dds, blind = FALSE)  # Takes group information into account
  head(assay(vsd), 3)
  colData(vsd) #size factor = the scaling factor DESeq2 uses to correct for differences between samples. Values close to 1.0 -> samples that need little normalization. Values below 1 -> samples with a high total read count (scaled down). Values above 1 -> samples with a low total read count (scaled up)
  dim(vsd)
  
  df_vsd<- as.data.frame(vsd@assays@data@listData)  #Extracts the VST-transformed count matrix (vsd@assays@data@listData) and converts it to a data.frame
  df_vsd$symbol <- mapIds(org.Hs.eg.db,
                          keys=as.character(rownames(df_vsd)),
                          column="SYMBOL",
                          keytype="ENTREZID",
                          multiVals="first")   #Looks up the SYMBOL for each gene and adds it as a new column in df_vsd.
  
  df_vsd_GSE114918_SNc <-  df_vsd
  write.csv(df_vsd, file = "normalizedcounts_GSE114918_SNC_supp.csv", append = FALSE, quote = TRUE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")
  
  # Save as Excel
  write.xlsx(df_vsd, file = "normalizedcounts_GSE114918_snc_supp.xlsx", rowNames = TRUE)
  
  
}
# PCA: visualize sample clustering and check for outliers or batch effects across conditions.
#PCA
{
  #pca
  pcaData1 <- plotPCA(vsd, intgroup=c("Group", "Sample"), returnData=TRUE)
  percentVar <- round(100 * attr(pcaData1, "percentVar"))
  
  ggplot(pcaData1, aes(PC1, PC2, color=Sample, shape=Group )) +
    geom_point(size = 10) +
    geom_text(aes(label = Sample), vjust = -1, size = 3, fontface = "bold") +  # Note this line
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, face="bold", size= 10)) +
    theme(axis.text.y = element_text(angle = 90, hjust = 1, face="bold", size= 10)) +
    theme(axis.title.x = element_text(face="bold",size=10)) +
    theme(axis.title.y = element_text(face="bold",size=10)) +
    theme(legend.text = element_text(size = 10, face="bold", color = "black"))+
    theme(legend.title = element_text(size = 10, face="bold", color = "black"))+
    theme(plot.background = element_rect(fill =  "white"),
          panel.background = element_rect(fill = "white"),
          panel.border = element_rect(color = "black", fill = NA, linetype = 'solid'),
          panel.grid = element_line(size = 0.5, linetype = 'solid', colour = "grey"))
  
  
  ggsave(paste0('figures/PCA_analiz_5__SNC_GSE114918.tiff'), 
         units="in", width=10, height=10, dpi=600, compression = 'lzw')
  

  
}
# Sample-distance heatmap: assess overall similarity/clustering between samples.
#heatmap
{
  vsd<-as.data.frame(vsd@assays@data@listData)
  head(vsd)
  
  my_sample_col<- read.csv("coldata_data2.csv", sep = ";")
  
  my_sample_col<-coldata[1]
  my_sample_col<- data.frame(my_sample_col)
  
  my_sample_col<- my_sample_col[,-2]
  
  row.names(my_sample_col) <- colnames(vsd)
  
  
  col<- colorRampPalette(c("yellow", "orange", "red"))(50)
  
  a<-pheatmap(vsd, cluster_rows = T, show_rownames=F, color = col,annotation_col = my_sample_col,
              fontsize =20, scale="row", treeheight_col= 20,
              treeheight_row= 20, fontsize_row=20,
              fontsize_col = 20)
  
  
  ggsave("figures/heatmap_analiz_5_sirt1_SNC_GSE114918.tiff",
         a$gtable,  # Save the grid object
         device = "tiff",
         width = 10,    # Width in inches
         height = 12,   # Height in inches
         
         dpi = 600) 
}
# Differential expression analysis with DESeq2: compare Parkinson's disease vs. control groups.
#Differential expression analysis
{
  dds <- DESeq(dds)
  colData(dds)
  resultsNames(dds)
  
  
  #Building the results table
  res <- results(dds, contrast = c("Group", "PD", "Healthy"))
  res
  res<- na.omit(res)
  
  summary(res)
  
  deseq<- as.data.frame(res@listData)
  rownames(deseq) = res@rownames
  deseq<- na.omit(deseq) 
  
  summary(res$pvalue)
  colData(dds)
  
  levels(dds$Group)
  dds$Group <- factor(trimws(dds$Group))
  levels(dds$Group)
  dds$Group <- relevel(dds$Group, ref = "Healthy")
  
  levels(dds$Group)
  levels(dds$Group)
  
  deseq$symbol <- mapIds(org.Hs.eg.db,
                         keys=as.character(rownames(deseq)),
                         column="SYMBOL",
                         keytype="ENTREZID",
                         multiVals="first")
  
  head(deseq)
  deseq_GSE114918_SNc <- deseq
  
  write.csv(deseq, file = "deseq_GSE114918_SNC_supp.csv", append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")
  library(openxlsx)  # Save as Excel
  write.xlsx(deseq, file = "deseq_GSE114918_SNC_supp.xlsx", rowNames = TRUE)
}
# Extract p-values and adjusted p-values (padj) from the DESeq2 results table.
#pval#padj
{  
  resSig <- deseq[deseq$padj < 0.05,]
  resSig_up <- subset(resSig, log2FoldChange > 0.6 & padj < 0.05)
  resSig_down <- subset(resSig, log2FoldChange < -0.6 & padj < 0.05)
  resSig<-rbind(resSig_up, resSig_down)
  
  write.csv(resSig, file = "deseq_significant_pval_fc_GSE114918_SNc_supp.csv", append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")  
  write.xlsx(resSig, file = "deseq_significant_pval_fc_GSE114918_SNc_supp.xlsx", rowNames = TRUE)
  
  
  # significant vsd counts
  vsd2 <- vst(dds, blind = FALSE)
  vsd2<-as.data.frame(vsd2@assays@data@listData)
  sig_vsd <- vsd2[rownames(resSig),]
  
  col<- colorRampPalette(c("blue", "white", "red"))(50)
  library(RColorBrewer)
  color <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)
  
  my_sample_col<- coldata
  my_sample_col<- data.frame(my_sample_col)
  my_sample_col<- my_sample_col[,-2]
  row.names(my_sample_col) <- colnames(sig_vsd)
  
  
  
  adfipsc <-pheatmap(sig_vsd, cluster_rows = T, show_rownames=F, color=col ,annotation_col =my_sample_col,
                     fontsize =5, scale="row", treeheight_col= 20,
                     treeheight_row= 20, fontsize_row=5,
                     fontsize_col = 10)
  
  
  ggsave("figures/heatmap_ggplot_significant_gse114918_outlier2_SNC.tiff",
         adfipsc$gtable,  # Save the grid object
         device = "tiff",
         width = 8,    # Width in inches
         height = 8,   # Height in inches
         dpi = 500) 
}
# Classify genes as up- or down-regulated based on fold-change direction and significance threshold.
#FC up and down
{ 
  
  resSig_up <- subset(resSig, log2FoldChange > 0.6 & padj < 0.05)
  resSig_up$symbol <- mapIds(org.Hs.eg.db,
                             keys=as.character(rownames(resSig_up)),
                             column="SYMBOL",
                             keytype="ENTREZID",
                             multiVals="first")
  write.csv(resSig_up, file = "deseq_significant_upregulated_data.csv", append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ",", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")
  
  resSig_down <- subset(resSig, log2FoldChange < -0.6 & padj < 0.05)
  resSig_down$symbol <- mapIds(org.Hs.eg.db,
                               keys=as.character(rownames(resSig_down)),
                               column="SYMBOL",
                               keytype="ENTREZID",
                               multiVals="first")
  write.csv(resSig_down, file = "deseq_significant_downregulated_data.csv", append = FALSE, quote = TRUE, sep = " ",
            eol = "\n", na = "NA", dec = ",", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")
} 
# Volcano plot: visualize significant differentially expressed genes (fold-change vs. significance).
#volcanoplot
{
  # ============================================================
  # VOLCANO PLOT (RNA-seq)
  # ============================================================
  
  # Determine the number of upregulated and downregulated genes
  n_up <- nrow(resSig_up)
  n_down <- nrow(resSig_down)
  
  # Create the volcano plot
  library(EnhancedVolcano)
  library(ggplot2)
  
  avolcanoipsc <- EnhancedVolcano(
    res,
    lab = deseq$symbol,
    x = 'log2FoldChange',
    y = 'padj',
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
    legendLabels = c('NS','log2FC','padj','padj & log2FC'),
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
  
  # Define the ggplot axis ranges
  x_lim <- ggplot_build(avolcanoipsc)$layout$panel_params[[1]]$x.range
  y_lim <- ggplot_build(avolcanoipsc)$layout$panel_params[[1]]$y.range
  
  # Add annotations symmetrically
  final_plot <- avolcanoipsc +
    annotate("text", 
             x = x_lim[2] * 0.85, y = y_lim[2] * 0.95,
             label = paste0("Up: ", n_up),
             size = 6, fontface = "bold", color = "black", hjust = 1) +
    annotate("text", 
             x = x_lim[1] * 0.85, y = y_lim[2] * 0.95,
             label = paste0("Down: ", n_down),
             size = 6, fontface = "bold", color = "black", hjust = 0) +
    theme_minimal(base_size = 16) +
    theme(
      plot.title = element_text(face = "bold", size = 22, hjust = 0.5),
      axis.title = element_text(face = "bold", size = 20),
      axis.text = element_text(face = "bold", size = 16),
      legend.title = element_text(face = "bold", size = 18),
      legend.text = element_text(size = 14),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "grey90", size = 0.3),
      panel.grid.minor = element_blank()
    )
  plot(final_plot)
  # Save the figure (600 dpi, supplementary standard)
  ggsave(
    filename = "figures/SUPP_VOLCANO_GSE285507 LUHMES_24hr_Pb_10_uM_supp_SUPP.tiff",
    plot = final_plot,
    width = 10, height = 10, dpi = 600, bg = "white"
  )
  
  
}
# Marker gene panel: plot expression of selected Parkinson's disease-related marker genes across groups.
#marker tidyplot 
{ 
  # SUPPLEMENT BARPLOT: Marker Gene Groups (log2FC & padj)
  # ============================================================
  
  library(dplyr)
  library(ggplot2)
  
  
  
  # ------------------------------------------------------------
  # Extract data for each marker group
  # ------------------------------------------------------------
  theme_supplement <- theme_minimal(base_size = 18) +
    theme(
      plot.title = element_text(face = "bold", size = 28, hjust = 0.5),
      axis.title = element_text(face = "bold", size = 24),
      axis.text = element_text(face = "bold", size = 20),
      legend.title = element_text(face = "bold", size = 22),
      legend.text = element_text(size = 18),
      
      panel.border = element_blank(),
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white", color = NA),
      
      panel.grid.major = element_line(color = "grey90", size = 0.3),
      panel.grid.minor = element_blank(),
      
      strip.text = element_text(face = "bold", size = 22),
      strip.background = element_blank()
    )
  
  
  parkinson_genes <- c('SNCA', 'LRRK2', 'PRKN','PINK1', 'GBA1',"GBA2","GBA3", 'VPS35', 'ATP13A2', 'PLA2G6', 'FBXO7')
  genes_of_lewy_body <- c("UBB", "HSPA1A", "HSP90AA1", "TGM2", "STUB1")
  sirt_genes <- c("SIRT1", "SIRT2", "SIRT3", "SIRT4", "SIRT5", "SIRT6", "SIRT7")
  dopamine_marker_genes <- c("TH",  "DDC", "SLC18A2", "NR4A2", "ALDH1A1", "KCNJ6", "CALB1", "FOXA1", "DRD2")
  tet_marker<- c("TET1","TET2","TET3")
  
  
  # ============================================================
  # SUPPLEMENT BARPLOT: Marker Gene Groups (optimized fonts)
  # ============================================================
  
  library(dplyr)
  library(ggplot2)
  
  # ------------------------------------------------------------
  # Extract data for each marker group
  # ------------------------------------------------------------
  get_group_data <- function(gene_set, group_name) {
    deseq %>%
      filter(symbol %in% gene_set) %>%
      select(symbol, log2FoldChange, padj) %>%
      
      mutate(Group = group_name)}
  get_markers <- function(deseq_Fc, gene_set, group_name) {
    deseq %>%
      filter(symbol %in% gene_set) %>%
      select(symbol, log2FoldChange, padj) %>%
      mutate(Group = group_name)
  }
  
  # ------------------------------------------------------------
  # Merge the groups
  # ------------------------------------------------------------
  df_markers <- bind_rows(
    get_group_data(parkinson_genes, "Parkinson"),
    get_group_data(dopamine_marker_genes, "Dopamine"),
    get_group_data(sirt_genes, "SIRT"),
    get_group_data(tet_marker, "TET"),
    get_group_data(genes_of_lewy_body, "LewyBody")
  )
  
  # ------------------------------------------------------------
  # Visualization (small-font version)
  # ------------------------------------------------------------
  marker_plot_small <- ggplot(df_markers, aes(
    x = reorder(symbol, log2FoldChange),
    y = log2FoldChange,
    fill = padj < 0.05 & abs(log2FoldChange) > 0.6
  )) +
    geom_bar(stat = "identity") +
    geom_text(
      aes(
        label = ifelse(!is.na(padj),
                       paste0("p=", signif(padj, 2)),
                       "")
      ),
      vjust = -0.4,
      size = 3,
      fontface = "bold",
      color = "black"
    ) +
    facet_wrap(~Group, scales = "free_x", nrow = 1) +
    scale_fill_manual(values = c("TRUE" = "firebrick3", "FALSE" = "grey70")) +
    labs(
      title = "Marker Gene Groups: log2FC and Significance (padj < 0.05 & |log2FC|> 0.6)",
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
        size = 18,     # large gene label font
        face = "bold"
      )
    )
  
  marker_plot_small
  
  # ------------------------------------------------------------
  # Save (600 DPI, supplementary standard)
  # ------------------------------------------------------------
  ggsave("figures/_SUPP_.tiff",
         marker_plot_small,
         width = 30, height = 10, dpi = 600, bg = "white")
}
# GO and KEGG enrichment analysis on the set of differentially expressed genes.
#GO KEGG
{
   
  # 4. GO ENRICHMENT DOTPLOT
  # ============================================================
  deg_genes <- as.character(rownames(resSig))
  head(rownames(resSig), 5)
  head(resSig)
  genes <- gsub("\\..*", "", rownames(resSig))
  head(genes)
  library(org.Hs.eg.db)
  
  resSig$ENTREZID <- mapIds(
    org.Hs.eg.db,
    keys = genes,
    column = "ENTREZID",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  resSig <- resSig[!is.na(resSig$ENTREZID), ]
  
  
  genes <- gsub("\\..*", "", rownames(deseq))
  head(genes)
  library(org.Hs.eg.db)
  
  deseq$ENTREZID <- mapIds(
    org.Hs.eg.db,
    keys = genes,
    column = "ENTREZID",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  deseq <- deseq[!is.na(deseq$ENTREZID), ]
  
  
  
  # 1) Remove entries where symbol is NA
  resSig <- resSig[!is.na(resSig$symbol) & resSig$symbol != "", ]
  
  # 2) (duplicate symbols may still remain) clean those up too
  resSig <- resSig[!duplicated(resSig$symbol), ]
  
  # 3) Now assign rownames
  rownames(resSig) <- resSig$ENTREZID
  library(clusterProfiler)
  library(org.Hs.eg.db)
  
  # 1) Drop NA entries
  deseq <- deseq[!is.na(deseq$ENTREZID), ]
  
  # 2) Remove duplicated ENTREZID entries (keeps the first one)
  deseq <- deseq[!duplicated(deseq$ENTREZID), ]
  
  # 3) Now rownames are assigned
  rownames(deseq) <- deseq$ENTREZID
  
  deseq$symbol <- mapIds(
    org.Hs.eg.db,
    keys      = rownames(deseq),
    column    = "ENTREZID",
    keytype   = "SYMBOL",
    multiVals = "first"
  )
  rownames(deseq) <- deseq$ENTREZID
  
  resSig$ENTREZID <- mapIds(
    org.Hs.eg.db,
    keys      = rownames(resSig),
    column    = "ENTREZID",
    keytype   = "SYMBOL",
    multiVals = "first"
  )
  resSig <- resSig[!is.na(resSig$ENTREZID), ]
  resSig <- resSig[!duplicated(resSig$ENTREZID), ]
  rownames(resSig) <- resSig$ENTREZID
  
  resSig <- resSig[!is.na(resSig$symbol), ]
  resSig <- resSig[!duplicated(resSig$symbol), ]
  rownames(resSig) <- resSig$symbol
  uni<-row.names(resSig)
  
  egoBP <- enrichGO(
    gene = as.character(row.names(resSig)),
    OrgDb = "org.Hs.eg.db",
    ont = "BP", universe = row.names(deseq),
    pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE
  )
  egoMF <- enrichGO(
    gene = as.character(row.names(resSig)),
    OrgDb = "org.Hs.eg.db",
    ont = "MF", universe = row.names(deseq),
    pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE
  )
  egoCC <- enrichGO(
    gene = as.character(row.names(resSig)),
    OrgDb = "org.Hs.eg.db",
    ont = "CC", universe = row.names(deseq),
    pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE
  )
  library(dplyr)
  # Combine
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
    arrange(p.adjust, desc(Count), .by_group = TRUE) |>
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
      title = "Top 10 Enriched GO Terms",
      x = "Gene Ratio",
      y = "GO Term Description",
      color = "Adjusted p-value",
      size = "Gene Count"
    ) +
    theme_supplement
  print(go_plot)
  ggsave("figures/SUPP_GO_GSE185009_IPSC_SUPP.tiff",
         go_plot, width = 30, height = 10, dpi = 600, bg = "white")
  
  
  write.csv(combined_go, file = "GSE285507 LUHMES 24hr_Pb_10_uM_supp_SUPP_GO_terms_SUPP.csv", append = FALSE, quote = TRUE, sep = "\t",
            eol = "\n", na = "NA", dec = ",", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")
  # Save as Excel
  write.xlsx(combined_go, file = "GSE196190_t74_200qm_SUPP_GO_terms_SUPP.xlsx", rowNames = TRUE)
  library(openxlsx)   
  
  # 1) Convert SYMBOL to ENTREZ ID
  gene_ids <- bitr(resSig$symbol,
                   fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)
  
  # 2) KEGG enrichment (Homo sapiens)
  kegg_result <- enrichKEGG(
    gene         = gene_ids$ENTREZID,
    organism     = 'hsa',       # human
    pvalueCutoff = 0.05
  )
  
  
  kegg_plot <- dotplot(kegg_result, showCategory = 10) +
    ggtitle("Top 10 Enriched KEGG Pathways") +
    theme_supplement
  print(kegg_plot)
  is.null(kegg_result)
  class(kegg_result)
  
  # As for the clusterProfiler result:
  nrow(as.data.frame(kegg_result))
  head(as.data.frame(kegg_result))
  
  ggsave("figures/GSE196190_t24_400qm_SUPP.tiff",
         kegg_plot, width = 10, height = 10, dpi = 600, bg = "white")
  write.csv(kegg_result, file = "GSE285507 LUHMES 24hr_Pb_10_uM_supp_SUPP_kegg_result_SUPP.csv", append = FALSE, quote = TRUE, sep = "\t",
            eol = "\n", na = "NA", dec = ",", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")
  # Save as Excel
  write.xlsx(kegg_result, file = "GSE196190_t24_400qm_SUPP_kegg_result_SUPP.xlsx", rowNames = TRUE)
  # ============================================================
  # GO + KEGG PLOT COMBINATION
  # ============================================================
  library(patchwork)
  library(stringr)
  
  go_plot <- go_plot +
    scale_y_discrete(labels = function(x) str_wrap(x, 25))
  print(go_plot)
  kegg_plot <- kegg_plot +
    scale_y_discrete(labels = function(x) str_wrap(x, 15))
  
  
  combined_plot <- go_plot + kegg_plot +
    plot_layout(widths = c(5, 1)) +
    plot_annotation(
      theme = theme(
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_blank()
      )
    )}
# STRING protein-protein interaction network: query STRING and compute node centrality metrics.
#string 
{
  # Package installations (may only be needed on the first run)
  install.packages(c("httr", "jsonlite", "igraph", "ggraph", "tidygraph", "ggplot2"))
  BiocManager::install(c("org.Hs.eg.db", "BioNet"))
  
  # Required libraries
  library(httr)
  library(jsonlite)
  library(igraph)
  library(ggraph)
  library(tidygraph)
  library(ggplot2)
  library(org.Hs.eg.db)
  library(BioNet)
  
  # Packages
  library(httr)
  library(jsonlite)
  library(igraph)
  library(ggraph)
  library(ggplot2)
  library(dplyr)
  # Fetch network data from the STRING API
  gene_list<- resSig$symbol
  
  library(httr)
  library(jsonlite)
  
  
  # API URL
  string_api_url <- "https://string-db.org/api/json/network"
  
  # Split the list into chunks of 50 genes
  batch_size <- 50
  gene_batches <- split(gene_list, ceiling(seq_along(gene_list) / batch_size))
  
  # Empty list to collect all data fetched from the API
  all_string_data <- list()
  for (i in seq_along(gene_batches)) {
    params <- list(
      identifiers = paste(gene_batches[[i]], collapse = "%0d"),
      species = 9606,
      required_score = 400
    )
    
    response <- GET(string_api_url, query = params)
    
    if (http_status(response)$category != "Success") {
      warning(paste("API Sorgusu Ba??ar??s??z! Hata Kodu:", http_status(response)$message))
      next
    }
    
    string_data <- fromJSON(content(response, as = "text", encoding = "UTF-8"))
    
    all_string_data[[i]] <- string_data
  }
  # Combine all the data
  string_data <- do.call(rbind, all_string_data)
  
  # Check the results
  head(string_data)
  library(igraph)
  library(httr)
  library(httr)
  library(ggraph)
  library(jsonlite)
  # Build the STRING network
  graph_data <- graph_from_data_frame(
    string_data[, c("preferredName_A", "preferredName_B", "score")], 
    directed = FALSE
  )
  
  # **Save data for Cytoscape**
  write.csv(string_data[, c("preferredName_A", "preferredName_B", "score")], 
            "figures/cytoscape_network_GSE185009", row.names = FALSE)
  
  # STRING network visualization
  ggraph(graph_data, layout = "fr") +
    geom_edge_link(aes(edge_alpha = score)) +
    geom_node_point(size = 3, color = "red") +
    geom_node_text(aes(label = name), repel = TRUE) +
    theme_void()
  
  ggraph(graph_data, layout = "fr") +
    geom_edge_link(aes(edge_alpha = score), color = "grey40") +  # assumption: 'score' is the edge weight
    geom_node_point(size = 3, color = "red") +
    geom_node_text(aes(label = name), repel = TRUE, size = 3) +
    theme_void() +
    ggtitle("STRING Network",
            subtitle = "Edge transparency reflects interaction score") +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5)
    )
  
  
  ggsave("figures/string_GSE114928_SNc_rapor_score400_s??rt1.tiff", width = 12, height = 10, dpi = 600, device = "tiff", bg="white")
  
  
  # Make sure the network has been created
  graph_data <- graph_from_data_frame(string_data[, c("preferredName_A", "preferredName_B")], directed = FALSE)
  
  # Degree (number of connections)
  degree_centrality <- igraph::degree(graph_data)
  head(degree_centrality) 
  
  # Betweenness Centrality (intermediate nodes)
  betweenness_centrality <- betweenness(graph_data, directed = FALSE, normalized = TRUE)
  
  # Closeness Centrality (distance)
  closeness_centrality <- closeness(graph_data, normalized = TRUE)
  
  # Eigenvector Centrality (connection strength)
  eigenvector_centrality <- eigen_centrality(graph_data)$vector
  
  # Let's save the results into a data frame
  network_metrics <- data.frame(
    Gene = V(graph_data)$name,
    Degree = degree_centrality,
    Betweenness = betweenness_centrality,
    Closeness = closeness_centrality,
    Eigenvector = eigenvector_centrality
  )
  
  # Let's rank them to identify the most important genes
  network_metrics <- network_metrics[order(-network_metrics$Betweenness), ]
  head(network_metrics)
  colnames(network_metrics)
  write.csv(network_metrics, "GSE114928_SNc_network_metrics.csv", row.names = TRUE)
  
  
  # Let's find drug-gene interactions using the DGIdb API
  drug_query <- paste(network_metrics$Gene, collapse = ",")
  drug_query  
  library(clusterProfiler)
  library(org.Hs.eg.db)
  
}  
# Spearman correlation between SIRT1 and SNCA expression across samples.
#correlation
{
  
  df <- df_vsd_GSE114918_SN  
  
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
    stop("SIRT1 veya SNCA bu datasette bulunamad??")
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

  
 



