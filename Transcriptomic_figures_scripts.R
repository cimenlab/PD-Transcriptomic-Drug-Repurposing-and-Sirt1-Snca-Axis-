---
scriptname: "Transcriptomic_figures_scripts.R"
title: "PD Transcriptomic Drug Repurposing and Sirt1 Snca Axis   "
author: "CoRD-OMICS - CimenLab - Fatma Hotaslier"
date: "21.07.2026"
license: "MIT License (see LICENSE file for details)"
description: "This script generates the three main summary figures of the manuscript by
combining differential expression results from all 33 RNA-seq and
microarray datasets. Figure 1 is a resveratrol target gene x dataset dot
plot; Figure 2 classifies each dataset's SIRT1-SNCA expression relationship
into distinct patterns; Figure 3 is a chord diagram linking resveratrol,
SIRT1, SNCA, and datasets by regulatory direction."
---

# FIGURE 1: resveratrol target gene x dataset dot plot. Combines DGIdb drug-gene interactions with DESeq2/limma results across all 33 datasets to show direction and significance of each target gene per dataset.
##FIGURE 1 SCRIPT
{
  # ---------------- DESeq results: Postmortem ----------------
  deseq_GSE156928        <- read.csv("deseq_GSE156928.csv", row.names = 1)
  deseq_GSE114918_SNc    <- read.csv("deseq_GSE114918_SNc.csv", row.names = 1)
  deseq_GSE114918_VTA    <- read.csv("deseq_GSE114918_VTA.csv", row.names = 1)
  deseq_GSE20292         <- read.csv("deseq_GSE20292.csv", row.names = 1)
  deseq_GSE20291         <- read.csv("deseq_GSE20291.csv", row.names = 1)
  
  # ---------------- DESeq results: iPSC ----------------
  deseq_GSE185009_ipsc   <- read.csv("deseq_GSE185009_ipsc.csv", row.names = 1)
  deseq_GSE185009_npc    <- read.csv("deseq_GSE185009_npc.csv", row.names = 1)
  
  deseq_GSE196190_t24_100 <- read.csv("deseq_GSE196190_t24_100qM.csv", row.names = 1)
  deseq_GSE196190_t24_200 <- read.csv("deseq_GSE196190_t24_200qM.csv", row.names = 1)
  deseq_GSE196190_t24_400 <- read.csv("deseq_GSE196190_t24_400qM.csv", row.names = 1)
  
  deseq_GSE196190_t74_100 <- read.csv("deseq_GSE196190_t74_100qM.csv", row.names = 1)
  deseq_GSE196190_t74_200 <- read.csv("deseq_GSE196190_t74_200qM.csv", row.names = 1)
  deseq_GSE196190_t74_400 <- read.csv("deseq_GSE196190_t74_400qM.csv", row.names = 1)
  
  deseq_GSE120746         <- read.csv("deseq_GSE120746.csv", row.names = 1)
  deseq_GSE282494         <- read.csv("deseq_GSE282494.csv", row.names = 1)
  deseq_GSE181029_NP      <- read.csv("deseq_GSE181029_NP.csv", row.names = 1)
  deseq_GSE181029_DN      <- read.csv("deseq_GSE181029_DN.csv", row.names = 1)
  deseq_GSE315738         <- read.csv("deseq_GSE315738.csv", row.names = 1)
  
  deseq_GSE36321          <- read.csv("deseq_GSE36321.csv", row.names = 1)
  deseq_GSE51922          <- read.csv("deseq_GSE51922.csv", row.names = 1)
  deseq_GSE51922_LRRK2    <- read.csv("deseq_GSE51922_LRRK2.csv", row.names = 1)
  
  # ---------------- DESeq results: Cell-based ----------------
  deseq_GSE285507         <- read.csv("deseq_GSE285507.csv", row.names = 1)
  deseq_GSE229460         <- read.csv("deseq_GSE229460.csv", row.names = 1)
  deseq_GSE203522_2D      <- read.csv("deseq_GSE203522_2D.csv", row.names = 1)
  deseq_GSE203522_3D      <- read.csv("deseq_GSE203522_3D.csv", row.names = 1)
  
  deseq_GSE315111_rot     <- read.csv("deseq_GSE315111_rot.csv", row.names = 1)
  deseq_GSE315111_MPP     <- read.csv("deseq_GSE315111_MPP.csv", row.names = 1)
  
  deseq_GSE21305          <- read.csv("deseq_GSE21305.csv", row.names = 1)
  deseq_GSE191302         <- read.csv("deseq_GSE191302.csv", row.names = 1)
  
  deseq_GSE35642_5mM      <- read.csv("deseq_GSE35642_5mM.csv", row.names = 1)
  deseq_GSE35642_50mM     <- read.csv("deseq_GSE35642_50mM.csv", row.names = 1)
  deseq_GSE287941_6OHDA     <- read.csv("deseq_GSE287941_6OHDA.csv", row.names = 1)
  deseq_GSE287941_6OHDA_asyn    <- read.csv("deseq_GSE287941_6OHDA+asyn.csv", row.names = 1)
  library(readxl)
  library(dplyr)
  
  # File paths
  file_post  <- "C:/Users/User/OneDrive/Belgeler/ana_figure_postmortem_03-06/All_Approved_Drug_Interactions_Combined_postmortem_ALL_DATA_supp_3_06_2026_LAST_supp.xlsx"
  file_ipsc <- "C:/Users/User/OneDrive/Belgeler/ana_figure_ipsc_03.06/All_Approved_Drug_Interactions_Combined__IPSC_ALL_DATA_supp_3_06_2026_LAST_supp.xlsx"
  file_cell <- "C:/Users/User/OneDrive/Belgeler/ana_figure_cell_based_03-06/All_Approved_Drug_Interactions_Combined__Cell_based_ALL_DATA_supp_3_06_2026_LAST_supp.xlsx"
  
  # Read the files
  cell_df <- read_excel(file_cell)
  ipsc_df <- read_excel(file_ipsc)
  post_df <- read_excel(file_post)
  
  # Filtering
  cell_resveratrol <- cell_df %>% filter(drug == "RESVERATROL")
  ipsc_resveratrol <- ipsc_df %>% filter(drug =="RESVERATROL")
  post_resveratrol <- post_df %>% filter(drug == "RESVERATROL")
  
  # Combine everything if you want
  
  all_resveratrol <- bind_rows(
    cell_based = cell_df %>% filter(drug == "RESVERATROL"),
    ipsc       = ipsc_df %>% filter(drug == "RESVERATROL"),
    postmortem = post_df %>% filter(drug == "RESVERATROL"),
    .id = "Group"
  ) %>%
    select(-any_of(c("regulatory approval", "indication", "regulatory.approval")))
  
  # Check the results
  cell_resveratrol
  ipsc_resveratrol
  post_resveratrol
  all_resveratrol
  
  library(dplyr)
  unique(all_resveratrol$dataset)
  head(all_resveratrol)
  ## =============================================================================
  ## STEP 0 --- DIAGNOSIS: show the ACTUAL column names of every deseq object
  ## Run this and paste ALL the output back. That way we can see which file uses
  ## log2FoldChange/logFC or another name, padj/adj.P.Val or another name,
  ## and add the correct overrides to the get_gene_stats() function
  ## accordingly --- we will not guess.
  ## =============================================================================
  
  deseq_object_names <- c(
    "deseq_GSE285507", "deseq_GSE229460", "deseq_GSE203522_2D", "deseq_GSE203522_3D",
    "deseq_GSE315111_rot", "deseq_GSE315111_MPP", "deseq_GSE21305", "deseq_GSE191302",
    "deseq_GSE35642_5mM", "deseq_GSE35642_50mM", "deseq_GSE287941_6OHDA",
    "deseq_GSE287941_6OHDA_asyn", "deseq_GSE156928", "deseq_GSE114918_SNc",
    "deseq_GSE114918_VTA", "deseq_GSE20292", "deseq_GSE20291",
    "deseq_GSE185009_ipsc", "deseq_GSE185009_npc",
    "deseq_GSE196190_t24_100", "deseq_GSE196190_t24_200", "deseq_GSE196190_t24_400",
    "deseq_GSE196190_t74_100", "deseq_GSE196190_t74_200", "deseq_GSE196190_t74_400",
    "deseq_GSE120746", "deseq_GSE282494", "deseq_GSE181029_NP", "deseq_GSE181029_DN",
    "deseq_GSE315738", "deseq_GSE36321", "deseq_GSE51922_LRRK2", "deseq_GSE51922"
  )
  
  column_report <- lapply(deseq_object_names, function(nm) {
    if (!exists(nm)) {
      return(data.frame(object = nm, exists = FALSE, columns = NA_character_,
                        n_rows = NA_integer_, row_names_example = NA_character_,
                        stringsAsFactors = FALSE))
    }
    df <- get(nm)
    data.frame(
      object            = nm,
      exists            = TRUE,
      columns           = paste(colnames(df), collapse = " | "),
      n_rows            = nrow(df),
      row_names_example = paste(head(rownames(df), 3), collapse = ", "),
      stringsAsFactors  = FALSE
    )
  })
  
  column_report_df <- do.call(rbind, column_report)
  print(column_report_df, right = FALSE)
  
  ## also write it to a file, to make it easy to copy
  write.csv(column_report_df, "deseq_column_report.csv", row.names = FALSE)
  

  print(table(column_report_df$columns))
  
  
  ## =============================================================================
  ## Figure 2 (REAL DATA) --- resveratrol target gene x dataset
  ## all_resveratrol (Group, gene, drug, interaction score, dataset) + deseq_*
  ## objects: log2FC/padj are pulled and direction (up/down) and significance are determined.
  ## =============================================================================
  ## ASSUMPTION: the following objects are already loaded in the environment:
  ##   - all_resveratrol (Group, gene, drug, `interaction score`, dataset)
  ##   - deseq_GSE285507, deseq_GSE229460, deseq_GSE203522_2D, deseq_GSE203522_3D,
  ##     deseq_GSE315111_rot, deseq_GSE315111_MPP, deseq_GSE21305, deseq_GSE191302,
  ##     deseq_GSE35642_5mM, deseq_GSE35642_50mM, deseq_GSE287941_6OHDA,
  ##     deseq_GSE287941_6OHDA_asyn, deseq_GSE156928, deseq_GSE114918_SNc,
  ##     deseq_GSE114918_VTA, deseq_GSE20292, deseq_GSE20291,
  ##     deseq_GSE185009_ipsc, deseq_GSE185009_npc,
  ##     deseq_GSE196190_t24_100, deseq_GSE196190_t24_200, deseq_GSE196190_t24_400,
  ##     deseq_GSE196190_t74_100, deseq_GSE196190_t74_200, deseq_GSE196190_t74_400,
  ##     deseq_GSE120746, deseq_GSE282494, deseq_GSE181029_NP, deseq_GSE181029_DN,
  ##     deseq_GSE315738, deseq_GSE36321, deseq_GSE51922_LRRK2, deseq_GSE51922
  ##
  ## PURPOSE: for EVERY ROW in all_resveratrol (i.e. rows where DGIdb already
  ## says "this gene is a resveratrol target in this dataset") --- pull that
  ## gene's padj and log2FC from the deseq file for that dataset and attach them.
  ## We are not interested in gene-dataset combinations that were never observed (no full grid).
  ## =============================================================================
  
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
  
  ## ---------------------------------------------------------------------------
  # Step 1: map each dataset label to its DESeq2/limma results object, a clean display label, and its model group (postmortem / iPSC / cell-based).
  ## 1) dataset string -> deseq object + clean label + model group
  ##    (built and manually verified from unique(all_resveratrol$dataset) --- 33/33)
  ## ---------------------------------------------------------------------------
  dataset_lookup <- tibble::tribble(
    ~dataset_raw,                                        ~deseq_name,                    ~dataset_clean,             ~group,
    "GSE285507\r\nLUHMES\r\nPb",                          "deseq_GSE285507",              "GSE285507",                "Cell-based",
    "GSE229460\r\nLUHMES\r\nMPP+",                        "deseq_GSE229460",              "GSE229460",                "Cell-based",
    "GSE203522\r\n2D\r\nSH-SY5Y\r\nMPP+",                 "deseq_GSE203522_2D",           "GSE203522 2D",             "Cell-based",
    "GSE203522\r\n3D\r\nSH-SY5Y\r\nMPP+",                 "deseq_GSE203522_3D",           "GSE203522 3D",             "Cell-based",
    "GSE315111\r\nSH-SY5Y\r\nrotenone",                   "deseq_GSE315111_rot",          "GSE315111 Rotenone",       "Cell-based",
    "GSE315111\r\nSH-SY5Y\r\nMPP+",                       "deseq_GSE315111_MPP",          "GSE315111 MPP+",           "Cell-based",
    "GSE21305\r\nSH-SY5Y\r\nparaquat",                    "deseq_GSE21305",               "GSE21305",                 "Cell-based",
    "GSE191302\r\nLUHMES\r\nasyn",                        "deseq_GSE191302",              "GSE191302",                "Cell-based",
    "GSE35642\r\nSK-N-MC\r\nrotenone\r\n5 mM",            "deseq_GSE35642_5mM",           "GSE35642 5nM",             "Cell-based",
    "GSE35642\r\nSK-N-MC\r\nrotenone\r\n50 mM",           "deseq_GSE35642_50mM",          "GSE35642 50nM",            "Cell-based",
    "GSE287941\r\nLUHMES\r\n6OHDA",                       "deseq_GSE287941_6OHDA",        "GSE287941 6-OHDA",         "Cell-based",
    "GSE287941\r\nLUHMES\r\n6OHDA+asyn",                  "deseq_GSE287941_6OHDA_asyn",   "GSE287941 6-OHDA+asyn",    "Cell-based",
    "GSE196190\r\nt74\r\n100 uM",                         "deseq_GSE196190_t74_100",      "GSE196190 t74 100uM",      "iPSC",
    "GSE196190\r\nt24\r\n400 uM",                         "deseq_GSE196190_t24_400",      "GSE196190 t24 400uM",      "iPSC",
    "GSE196190\r\nt24\r\n200 uM",                         "deseq_GSE196190_t24_200",      "GSE196190 t24 200uM",      "iPSC",
    "GSE196190\r\nt24\r\n100 uM",                         "deseq_GSE196190_t24_100",      "GSE196190 t24 100uM",      "iPSC",
    "GSE185009\r\nNeural\r\nprogenitor cells",            "deseq_GSE185009_npc",          "GSE185009 NPC",            "iPSC",
    "GSE185009\r\niPSC",                                  "deseq_GSE185009_ipsc",         "GSE185009 iPSC",          "iPSC",
    "GSE51922\r\nLRRK2\r\nmutation",                      "deseq_GSE51922_LRRK2",         "GSE51922 LRRK2",           "iPSC",
    "GSE51922\r\nidiopatic",                              "deseq_GSE51922",               "GSE51922 idiopathic",      "iPSC",
    "GSE36321\r\nNeural stem\r\ncell",                    "deseq_GSE36321",               "GSE36321",                 "iPSC",
    "GSE282494\r\nLRRK2\r\nmutation",                     "deseq_GSE282494",              "GSE282494",                "iPSC",
    "GSE181029\r\nNeural\r\nprogenitor",                  "deseq_GSE181029_NP",           "GSE181029 NP",             "iPSC",
    "GSE181029\r\nDopaminergic\r\nneuron",                "deseq_GSE181029_DN",           "GSE181029 DN",             "iPSC",
    "GSE315738\r\nGBA",                                   "deseq_GSE315738",              "GSE315738",                "iPSC",
    "GSE120746\r\nDopaminergic\r\nneuron",                "deseq_GSE120746",              "GSE120746",                "iPSC",
    "GSE196190\r\nt74\r\n400 uM",                         "deseq_GSE196190_t74_400",      "GSE196190 t74 400uM",      "iPSC",
    "GSE196190\r\nt74\r\n200 uM",                         "deseq_GSE196190_t74_200",      "GSE196190 t74 200uM",      "iPSC",
    "GSE156928\r\nFrontal Cortex",                        "deseq_GSE156928",              "GSE156928",                "Postmortem",
    "GSE114918\r\nVTA",                                   "deseq_GSE114918_VTA",          "GSE114918 VTA",            "Postmortem",
    "GSE114918\r\nSNc",                                   "deseq_GSE114918_SNc",          "GSE114918 SNc",            "Postmortem",
    "GSE20292\r\nSN",                                     "deseq_GSE20292",               "GSE20292",                 "Postmortem",
    "GSE20291\r\nPutamen",                                "deseq_GSE20291",               "GSE20291",                 "Postmortem"
  )
  
  stopifnot(nrow(dataset_lookup) == 33)
  
  ## check: does every dataset_raw value in all_resveratrol exist in the lookup?
  missing_datasets <- setdiff(unique(all_resveratrol$dataset), dataset_lookup$dataset_raw)
  if (length(missing_datasets) > 0) {
    warning("There are dataset strings not present in the lookup, please check:",
            paste(missing_datasets, collapse = "\n---\n"))
  }
  
  ## ---------------------------------------------------------------------------
  # Step 2: generic function to extract (gene, log2FC, padj) from a DESeq2 or limma results object, handling differing column names between platforms.
  ## 2) generic function to extract (gene, log2FC, padj) from each deseq object.
  ##    column names VERIFIED against the STEP0 diagnosis:
  ##      - RNA-seq (DESeq2): log2FoldChange, padj, symbol (sometimes SYMBOL/GeneSymbol)
  ##      - Microarray (limma): logFC, adj.P.Val --- if there is NO symbol column,
  ##        rownames already contain the gene symbol (e.g. EFHD1, CYTL1, DCC), fall back to that.
  ##    GSE287941_6OHDA / _asyn files have BOTH log2FoldChange/padj AND
  ##    logFC/adj.P.Val --- log2FoldChange/padj (DESeq2) is chosen with PRIORITY.
  ##    !! For these two, confirm what the logFC column actually represents !!
  ## ---------------------------------------------------------------------------
  get_gene_stats <- function(gene_symbol, deseq_obj_name) {
    if (!exists(deseq_obj_name)) {
      return(tibble::tibble(log2FC = NA_real_, padj = NA_real_, stat_source = NA_character_,
                            found = FALSE, note = paste("Object not found/not loaded:", deseq_obj_name)))
    }
    df <- get(deseq_obj_name)
    
    ## detect the column names --- list verified against STEP0
    lfc_col  <- intersect(c("log2FoldChange", "logFC"), colnames(df))
    padj_col <- intersect(c("padj", "adj.P.Val"), colnames(df))
    
    if (length(lfc_col) == 0 || length(padj_col) == 0) {
      return(tibble::tibble(log2FC = NA_real_, padj = NA_real_, stat_source = NA_character_,
                            found = FALSE,
                            note = paste0(deseq_obj_name, ": expected columns not found. Available columns: ",
                                          paste(colnames(df), collapse = ", "))))
    }
    lfc_col  <- lfc_col[1]   # takes log2FoldChange if present, otherwise logFC
    padj_col <- padj_col[1]  # takes padj if present, otherwise adj.P.Val
    
    ## --- find the gene name: FIRST the symbol/SYMBOL/GeneSymbol column, otherwise rownames ---
    symbol_col <- intersect(c("symbol", "SYMBOL", "GeneSymbol", "Symbol", "gene_symbol"), colnames(df))
    
    if (length(symbol_col) > 0) {
      symbol_col <- symbol_col[1]
      gene_vec <- trimws(as.character(df[[symbol_col]]))
      matches <- which(toupper(gene_vec) == toupper(trimws(gene_symbol)))
    } else {
      rn <- trimws(rownames(df))
      matches <- which(toupper(rn) == toupper(trimws(gene_symbol)))
    }
    
    if (length(matches) == 0) {
      return(tibble::tibble(log2FC = NA_real_, padj = NA_real_, stat_source = paste0(lfc_col, " / ", padj_col),
                            found = FALSE,
                            note = paste0(gene_symbol, " -> ", deseq_obj_name, " (may be outside the platform's coverage)")))
    }
    
    note_txt <- NA_character_
    if (length(matches) > 1) {
      ## if the same symbol appears in multiple rows (multiple transcripts/probes) --- pick
      ## the one with the lowest padj (the most significant hit) and note the status
      padj_vals <- suppressWarnings(as.numeric(df[[padj_col]][matches]))
      best <- matches[which.min(padj_vals)]
      note_txt <- paste0(gene_symbol, " icin ", length(matches), "N rows found for gene_symbol, lowest padj selected (row best)")
      idx <- best
    } else {
      idx <- matches
    }
    
    tibble::tibble(
      log2FC      = as.numeric(df[[lfc_col]][idx]),
      padj        = as.numeric(df[[padj_col]][idx]),
      stat_source = paste0(lfc_col, " / ", padj_col),
      found       = TRUE,
      note        = note_txt
    )
  }
  
  ## ---------------------------------------------------------------------------
  # Step 3: clean the resveratrol drug-gene interaction table, join it with the DE results lookup, and attach log2FC/padj to each gene-dataset pair.
  ## 3) clean all_resveratrol + merge with the lookup + for each row
  ##    pull log2FC/padj
  ## ---------------------------------------------------------------------------
  resveratrol_joined <- all_resveratrol %>%
    rename(interaction_score = `interaction score`) %>%
    left_join(dataset_lookup, by = c("dataset" = "dataset_raw"))
  
  ## check whether there are unmatched rows
  unmatched <- resveratrol_joined %>% filter(is.na(deseq_name))
  if (nrow(unmatched) > 0) {
    warning(nrow(unmatched),)
    print(unique(unmatched$dataset))
  }
  
  stats_list <- purrr::pmap(
    list(resveratrol_joined$gene, resveratrol_joined$deseq_name),
    function(g, d) get_gene_stats(g, d)
  )
  stats_df <- bind_rows(stats_list)
  
  resveratrol_full <- bind_cols(resveratrol_joined, stats_df) %>%
    mutate(
      direction = case_when(
        is.na(log2FC)       ~ NA_character_,
        log2FC > 0           ~ "Upregulated",
        log2FC < 0           ~ "Downregulated",
        TRUE                 ~ NA_character_
      ),
      significant = !is.na(padj) & padj < 0.05 & abs(log2FC) > 0.6
    ) %>%
    ## in the order you requested: gene, dataset (clean label), group, drug,
    ## interaction score, then log2FC/padj/direction/significance from deseq
    select(gene, dataset_clean, group, drug, interaction_score,
           log2FC, padj, direction, significant,
           stat_source, dataset_raw = dataset, deseq_name, found, note)
  
  ## ---------------------------------------------------------------------------
  # Step 4: data quality report ??? check for genes/datasets that failed to match between the interaction table and the DE results.
  ## 4) data quality report --- DO NOT SKIP THIS SILENTLY, read it first
  ## ---------------------------------------------------------------------------
 
  
  problem_rows <- resveratrol_full %>% filter(!found | (!is.na(padj) & !significant))
  if (nrow(problem_rows) > 0) {
    cat("\nAsagidaki satirlari elle kontrol edin (gene, dataset_clean, log2FC, padj, note):\n")
    print(problem_rows %>% select(gene, dataset_clean, log2FC, padj, found, note), n = Inf)
  }
  
  ## save this table so you don't need to re-read the deseq files repeatedly
  write.csv(resveratrol_full, "resveratrol_full_with_stats_all_datasets_22.07.2026.csv", row.names = FALSE)
  
  ## =============================================================================
  ## everything from here on is the same as the previous dot-plot script --- THE ONLY DIFFERENCE:
  ## `resveratrol_full` now contains REAL data, so instead of `gene_dataset_long`
  ## we use it directly. The presence_list / set.seed / sample() blocks
  ## have been COMPLETELY REMOVED.
  ## =============================================================================
  ## =============================================================================
  ## Figure 2 (FINAL, REAL DATA) --- resveratrol target gene x dataset
  ## Input: resveratrol_full (gene, dataset_clean, group, drug, interaction_score,
  ##        log2FC, padj, direction, significant, stat_source, dataset_raw,
  ##        deseq_name, found, note)
  ## Y axis = 33 datasets (grouped by model group)
  ## X axis = target genes (ordered by frequency of significant hits)
  ## Color   = direction (red=up / blue=down / grey=not significant or no data)
  ## Size    = interaction score
  ## Top bar = frequency, stacked by model group
  ## =============================================================================
  
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  
  ## ---------------------------------------------------------------------------
  ## 0) Safety: if there are multiple rows for the same (gene, dataset) pair
  ##    (e.g. the same gene came with multiple DGIdb records in the same dataset)
  ##    keep the one with the highest interaction_score.
  ## ---------------------------------------------------------------------------
  resveratrol_clean <- resveratrol_full %>%
    filter(!is.na(gene), !is.na(dataset_clean)) %>%
    group_by(gene, dataset_clean) %>%
    slice_max(order_by = interaction_score, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  ## ---------------------------------------------------------------------------
  ## 1) full list of the 33 datasets + model group (derived from resveratrol_full)
  ## ---------------------------------------------------------------------------
  datasets <- resveratrol_clean %>%
    distinct(dataset_clean, group) %>%
    mutate(group = factor(group, levels = c("Postmortem", "iPSC", "Cell-based"))) %>%
    arrange(group, dataset_clean)
  
  n_ds <- nrow(datasets)
  cat("Toplam veri seti sayisi (resveratrol_full'dan):", n_ds, "\n")
  if (n_ds != 33) warning("not 33 ", n_ds, " Please check")
  
  ## ---------------------------------------------------------------------------
  ## 2) determine the status for each (gene, dataset) cell
  ##    !! isTRUE() does not work vectorized, so in the previous version every row
  ##    was incorrectly set to "Not significant". Fixed: significant is already
  ##    a logical (TRUE/FALSE/NA) column, we use it directly.
  ## ---------------------------------------------------------------------------
  resveratrol_status <- resveratrol_clean %>%
    mutate(
      sig_clean = !is.na(significant) & significant,   # NA-safe
      status = case_when(
        sig_clean & log2FC > 0 ~ "Upregulated",
        sig_clean & log2FC < 0 ~ "Downregulated",
        TRUE                   ~ "Not significant"
      )
    )
  
  cat("Status dagilimi:\n"); print(table(resveratrol_status$status, useNA = "ifany"))
  
  ## ---------------------------------------------------------------------------
  ## 3) gene frequencies --- using ONLY the significant (Upregulated/Downregulated) counts
  ## ---------------------------------------------------------------------------
  gene_freq <- resveratrol_status %>%
    filter(status != "Not significant") %>%
    count(gene, name = "freq") %>%
    arrange(desc(freq))
  
  all_genes <- gene_freq$gene   # X-axis order (descending by frequency)
  

  
  ## ---------------------------------------------------------------------------
  ## 4) FULL GRID: every gene x every dataset --- missing combinations become "Not significant"
  ## ---------------------------------------------------------------------------
  full_grid <- expand_grid(gene = all_genes, dataset_clean = datasets$dataset_clean)
  
  plot_df <- full_grid %>%
    left_join(
      resveratrol_status %>% select(gene, dataset_clean, status, interaction_score, log2FC, padj),
      by = c("gene", "dataset_clean")
    ) %>%
    mutate(status = ifelse(is.na(status), "Not significant", status)) %>%
    left_join(datasets, by = "dataset_clean") %>%
    mutate(
      gene          = factor(gene, levels = all_genes),
      dataset_clean = factor(dataset_clean, levels = rev(datasets$dataset_clean)),  # Postmortem on top
      status        = factor(status, levels = c("Upregulated", "Downregulated", "Not significant")),
      is_sirt1      = gene == "SIRT1",
      is_sig        = status != "Not significant",
      log2FC_color  = ifelse(is_sig, log2FC, NA_real_),
      size_val      = ifelse(is_sig, interaction_score, 0.15)
    )
  
  raw_max_abs_lfc <- max(abs(plot_df$log2FC_color), na.rm = TRUE)
  cap_lfc <- max(2, quantile(abs(plot_df$log2FC_color), 0.90, na.rm = TRUE))
  cat("Ham max |log2FC|:", round(raw_max_abs_lfc, 2), " --- renk skalasi icin kullanilan sinir (90. persentil, min 2):",
      round(cap_lfc, 2), "\n")
  
  group_colors  <- c("Postmortem" = "#1E8449", "iPSC" = "#D4AC0D", "Cell-based" = "#6C3483")
  
  ## ---------------------------------------------------------------------------
  # Step 5: build the main dot plot (x = target gene, y = dataset), colored by regulation direction and sized by interaction score.
  ## 5) main dot plot --- X = gene, Y = dataset (updated to 12pt)
  ## ---------------------------------------------------------------------------
  p_dot <- ggplot(plot_df, aes(x = gene, y = dataset_clean)) +
    geom_vline(xintercept = seq_along(levels(plot_df$gene)), color = "grey90", linewidth = 0.3) +
    geom_hline(yintercept = seq_along(levels(plot_df$dataset_clean)), color = "grey90", linewidth = 0.3) +
    geom_point(data = filter(plot_df, is_sig),
               aes(color = log2FC_color, size = size_val)) +
    scale_color_gradient2(low = "#0033CC", mid = "white", high = "#E60000", midpoint = 0,
                          limits = c(-cap_lfc, cap_lfc), oob = scales::squish,
                          name = "log2FC") +
    scale_size_continuous(range = c(2, 9), name = "Interaction\nscore") +
    scale_x_discrete(expand = expansion(add = 1)) +
    scale_y_discrete(expand = expansion(add = 1)) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5,
                                 face = "bold", color = "black"),
      axis.text.y = element_text(size = 12, face = "bold", color = "black"),
      panel.grid = element_blank(),
      legend.position = "right",
      plot.margin = margin(t = 5, r = 10, b = 5, l = 5)
    )
  
  ## ---------------------------------------------------------------------------
  # Step 6: top bar chart ??? stacked count of significant hits per gene, grouped by model group.
  ## 6) top bar: each gene's SIGNIFICANT frequency, stacked by model group
  ## ---------------------------------------------------------------------------
  bar_df <- resveratrol_status %>%
    filter(status != "Not significant") %>%
    mutate(group = factor(group, levels = levels(datasets$group))) %>%
    count(gene, group, name = "n") %>%
    right_join(expand_grid(gene = all_genes, group = levels(datasets$group)), by = c("gene", "group")) %>%
    mutate(n = ifelse(is.na(n), 0, n),
           gene = factor(gene, levels = all_genes))
  
  totals_df <- gene_freq %>% mutate(gene = factor(gene, levels = all_genes))
  
  p_bar <- ggplot(bar_df, aes(x = gene, y = n, fill = group)) +
    geom_col(width = 0.65, color = "white", linewidth = 0.3) +
    scale_fill_manual(values = group_colors, name = "Model group") +
    geom_text(data = totals_df, aes(x = gene, y = freq, label = freq),
              inherit.aes = FALSE, vjust = -0.4, size = 4.2, fontface = "bold") +
    labs(
      title = "",
      subtitle = paste0(""),
      y = NULL, x = NULL
    ) +
    coord_cartesian(clip = "off") +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.y = element_text(size = 12),
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 11, color = "grey30"),
      legend.position = "right"
    )
  
  ## ---------------------------------------------------------------------------
  # Step 7: left-side strip indicating each dataset's model group, aligned with the dot plot rows.
  ## 7) left-hand strip: model group, aligned with the dataset rows
  ## ---------------------------------------------------------------------------
  strip_df <- datasets %>% mutate(dataset_clean = factor(dataset_clean, levels = rev(datasets$dataset_clean)))
  
  p_strip <- ggplot(strip_df, aes(x = 1, y = dataset_clean, fill = group)) +
    geom_tile(color = "white") +
    scale_fill_manual(values = group_colors, guide = "none") +
    theme_void()
  
  ## ---------------------------------------------------------------------------
  # Step 8: combine the dot plot, top bar, and side strip into one composite figure and save it to file.
  ## 8) combine and save
  ## ---------------------------------------------------------------------------
  design <- "
AB
CD
"
  
  layout <- plot_spacer() + p_bar + p_strip + p_dot +
    plot_layout(design = design, widths = c(0.6, 10), heights = c(2.4, 10), guides = "collect") &
    theme(legend.position = "right", legend.box = "vertical")
  
  print(layout)
  
  ## SAVING --- 33 datasets x 31 genes is cramped at 12pt, so we export wide and landscape
  output_dir  <- "figures"
  output_name <- "Figure1_resveratrol_target_dotplot_v5"
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  ggsave(file.path(output_dir, paste0(output_name, ".tiff")),
         plot = layout, width = 30, height = 20, units = "cm",
         dpi = 600, device = "tiff", compression = "lzw", bg = "white")
  ## ---------------------------------------------------------------------------
  ## Output folder / file name
  ## ---------------------------------------------------------------------------
  output_dir  <- "figures"
  output_name <- "Figure1_resveratrol_target_dotplot_v3"
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  ggsave(file.path(output_dir, paste0(output_name, ".tiff")),
         plot = layout, width = 18, height = 16, units = "cm",
         dpi = 600, device = "tiff", compression = "lzw", bg = "white")
  ggsave(file.path(output_dir, paste0(output_name, ".png")),
         plot = layout, width = 16, height = 12, dpi = 300, bg = "white")
  
  cat("\nKaydedildi:\n",
      " - ", file.path(output_dir, paste0(output_name, ".tiff")), " (600 dpi, yayin icin)\n",
      " - ", file.path(output_dir, paste0(output_name, ".png")),  " (300 dpi, on izleme icin)\n", sep = "")
  
}
# FIGURE 2: SIRT1-SNCA expression pattern classification across all 33 datasets, based on differential expression and correlation results.
##FIGURE 2 SCRIPT
{
  ## =============================================================================
  ## SIRT1-SNCA PATTERN ANALYSIS --- SINGLE BLOCK
  ## Section 1: reading ALL files (33 deseq + 33 df_vsd = 66 files)
  ## Section 2: pattern classification (producing pattern_data)
  ## Section 3: final figure (Figure 3, grouped by pattern)
  ## =============================================================================
  
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  
  ## =============================================================================
  # Section 1: load all 66 result files (33 DESeq2/limma tables + 33 normalized expression matrices) needed for the pattern analysis.
  ## SECTION 1 --- FILE READING
  ## Make sure your working directory (setwd) is set to the folder containing the files,
  ## or prepend the full path to the file names below.
  ## =============================================================================
  
  # ---------------- Postmortem (5) ----------------
  deseq_GSE156928        <- read.csv("deseq_GSE156928.csv", row.names = 1)
  deseq_GSE114918_SNc    <- read.csv("deseq_GSE114918_SNc.csv", row.names = 1)
  deseq_GSE114918_VTA    <- read.csv("deseq_GSE114918_VTA.csv", row.names = 1)
  deseq_GSE20292         <- read.csv("deseq_GSE20292.csv", row.names = 1)
  deseq_GSE20291         <- read.csv("deseq_GSE20291.csv", row.names = 1)
  
  df_vsd_GSE156928       <- read.csv("df_vsd_GSE156928.csv", row.names = 1)
  df_vsd_GSE114918_SNc   <- read.csv("df_vsd_GSE114918_SNc.csv", row.names = 1)
  df_vsd_GSE114918_VTA   <- read.csv("df_vsd_GSE114918_VTA.csv", row.names = 1)
  df_vsd_GSE20292        <- read.csv("df_vsd_GSE20292.csv", row.names = 1)
  df_vsd_GSE20291        <- read.csv("df_vsd_GSE20291.csv", row.names = 1)
  
  # ---------------- iPSC (16) ----------------
  deseq_GSE185009_ipsc   <- read.csv("deseq_GSE185009_ipsc.csv", row.names = 1)
  deseq_GSE185009_npc    <- read.csv("deseq_GSE185009_npc.csv", row.names = 1)
  deseq_GSE196190_t24_100 <- read.csv("deseq_GSE196190_t24_100qM.csv", row.names = 1)
  deseq_GSE196190_t24_200 <- read.csv("deseq_GSE196190_t24_200qM.csv", row.names = 1)
  deseq_GSE196190_t24_400 <- read.csv("deseq_GSE196190_t24_400qM.csv", row.names = 1)
  deseq_GSE196190_t74_100 <- read.csv("deseq_GSE196190_t74_100qM.csv", row.names = 1)
  deseq_GSE196190_t74_200 <- read.csv("deseq_GSE196190_t74_200qM.csv", row.names = 1)
  deseq_GSE196190_t74_400 <- read.csv("deseq_GSE196190_t74_400qM.csv", row.names = 1)
  deseq_GSE120746         <- read.csv("deseq_GSE120746.csv", row.names = 1)
  deseq_GSE282494         <- read.csv("deseq_GSE282494.csv", row.names = 1)
  deseq_GSE181029_NP      <- read.csv("deseq_GSE181029_NP.csv", row.names = 1)
  deseq_GSE181029_DN      <- read.csv("deseq_GSE181029_DN.csv", row.names = 1)
  deseq_GSE315738         <- read.csv("deseq_GSE315738.csv", row.names = 1)
  deseq_GSE36321          <- read.csv("deseq_GSE36321.csv", row.names = 1)
  deseq_GSE51922          <- read.csv("deseq_GSE51922.csv", row.names = 1)
  deseq_GSE51922_LRRK2    <- read.csv("deseq_GSE51922_LRRK2.csv", row.names = 1)
  
  df_vsd_GSE185009_ipsc      <- read.csv("df_vsd_GSE185009_ipsc.csv", row.names = 1)
  df_vsd_GSE185009_npc       <- read.csv("df_vsd_GSE185009_npc.csv", row.names = 1)
  df_vsd_GSE196190_t24_100qM <- read.csv("df_vsd_GSE196190_t24_100qM.csv", row.names = 1)
  df_vsd_GSE196190_t24_200qM <- read.csv("df_vsd_GSE196190_t24_200qM.csv", row.names = 1)
  df_vsd_GSE196190_t24_400qM <- read.csv("df_vsd_GSE196190_t24_400qM.csv", row.names = 1)
  df_vsd_GSE196190_t74_100qM <- read.csv("df_vsd_GSE196190_t74_100qM.csv", row.names = 1)
  df_vsd_GSE196190_t74_200qM <- read.csv("df_vsd_GSE196190_t74_200qM.csv", row.names = 1)
  df_vsd_GSE196190_t74_400qM <- read.csv("df_vsd_GSE196190_t74_400qM.csv", row.names = 1)
  df_vsd_GSE120746           <- read.csv("df_vsd_GSE120746.csv", row.names = 1)
  df_vsd_GSE282494           <- read.csv("df_vsd_GSE282494.csv", row.names = 1)
  df_vsd_GSE181029_NP        <- read.csv("df_vsd_GSE181029_NP.csv", row.names = 1)
  df_vsd_GSE181029_DN        <- read.csv("df_vsd_GSE181029_DN.csv", row.names = 1)
  df_vsd_GSE315738           <- read.csv("df_vsd_GSE315738.csv", row.names = 1)
  df_vsd_GSE36321            <- read.csv("df_vsd_GSE36321.csv", row.names = 1)
  df_vsd_GSE51922            <- read.csv("df_vsd_GSE51922.csv", row.names = 1)
  df_vsd_GSE51922_LRRK2      <- read.csv("df_vsd_GSE51922_LRRK2.csv", row.names = 1)
  
  # ---------------- Cell-based (12) ----------------
  deseq_GSE285507         <- read.csv("deseq_GSE285507.csv", row.names = 1)
  deseq_GSE229460         <- read.csv("deseq_GSE229460.csv", row.names = 1)
  deseq_GSE203522_2D      <- read.csv("deseq_GSE203522_2D.csv", row.names = 1)
  deseq_GSE203522_3D      <- read.csv("deseq_GSE203522_3D.csv", row.names = 1)
  deseq_GSE315111_rot     <- read.csv("deseq_GSE315111_rot.csv", row.names = 1)
  deseq_GSE315111_MPP     <- read.csv("deseq_GSE315111_MPP.csv", row.names = 1)
  deseq_GSE21305          <- read.csv("deseq_GSE21305.csv", row.names = 1)
  deseq_GSE191302         <- read.csv("deseq_GSE191302.csv", row.names = 1)
  deseq_GSE35642_5mM      <- read.csv("deseq_GSE35642_5mM.csv", row.names = 1)
  deseq_GSE35642_50mM     <- read.csv("deseq_GSE35642_50mM.csv", row.names = 1)
  deseq_GSE287941_6OHDA        <- read.csv("deseq_GSE287941_6OHDA.csv", row.names = 1)
  deseq_GSE287941_6OHDA_asyn   <- read.csv("deseq_GSE287941_6OHDA+asyn.csv", row.names = 1)
  
  df_vsd_GSE285507        <- read.csv("df_vsd_GSE285507.csv", row.names = 1)
  df_vsd_GSE229460        <- read.csv("df_vsd_GSE229460.csv", row.names = 1)
  df_vsd_GSE203522_2D     <- read.csv("df_vsd_GSE203522_2D.csv", row.names = 1)
  df_vsd_GSE203522_3D     <- read.csv("df_vsd_GSE203522_3D.csv", row.names = 1)
  df_vsd_GSE315111_rot    <- read.csv("df_vsd_GSE315111_rot.csv", row.names = 1)
  df_vsd_GSE315111_MPP    <- read.csv("df_vsd_GSE315111_MPP.csv", row.names = 1)
  df_vsd_GSE21305         <- read.csv("df_vsd_GSE21305.csv", row.names = 1)
  df_vsd_GSE191302        <- read.csv("df_vsd_GSE191302.csv", row.names = 1)
  df_vsd_GSE35642_5mM     <- read.csv("df_vsd_GSE35642_5mM.csv", row.names = 1)
  df_vsd_GSE35642_50mM    <- read.csv("df_vsd_GSE35642_50mM.csv", row.names = 1)
  df_vsd_GSE287941_6OHDA       <- read.csv("df_vsd_GSE287941_6OHDA.csv", row.names = 1)
  df_vsd_GSE287941_6OHDA_asyn  <- read.csv("df_vsd_GSE287941_6OHDA+asyn.csv", row.names = 1)
  
  cat("Bolum 1 tamam: 33 deseq + 33 df_vsd objesi okundu.\n")
  
  ## =============================================================================
  # Section 2: classify each dataset's SIRT1-SNCA relationship into a named expression pattern based on significance and correlation sign.
  ## SECTION 2 --- PATTERN CLASSIFICATION
  ## =============================================================================
  
  SIG_PADJ <- 0.05
  SIG_LFC  <- 0.6
  RHO_MIN  <- 0   ## <-- lowered from 0.2 to 0 at the user's request; there is no longer a |rho|
  ##     threshold, only significance (sirt1_sig/snca_sig) and
  ##     the availability of rho are checked. Even very weak (noise-level)
  ##     correlations now fall into a pattern.
  ##     "Unresolved" now only occurs when no gene is significant
  ##     or rho could not be computed (NA).
  
  ds_registry <- tibble::tribble(
    ~dataset_clean,              ~model_group,  ~deseq_obj,                    ~vsd_obj,                        ~fc_col,           ~p_col,       ~symbol_col,
    "GSE156928",                 "Postmortem",  "deseq_GSE156928",             "df_vsd_GSE156928",              "log2FoldChange",  "padj",       "symbol",
    "GSE114918 SNc",             "Postmortem",  "deseq_GSE114918_SNc",         "df_vsd_GSE114918_SNc",          "log2FoldChange",  "padj",       "symbol",
    "GSE114918 VTA",             "Postmortem",  "deseq_GSE114918_VTA",         "df_vsd_GSE114918_VTA",          "log2FoldChange",  "padj",       "symbol",
    "GSE20291",                  "Postmortem",  "deseq_GSE20291",              "df_vsd_GSE20291",               "logFC",            "adj.P.Val",  "symbol",
    "GSE20292",                  "Postmortem",  "deseq_GSE20292",              "df_vsd_GSE20292",               "logFC",            "adj.P.Val",  "symbol",
    
    "GSE185009 iPSC",            "iPSC",        "deseq_GSE185009_ipsc",        "df_vsd_GSE185009_ipsc",         "log2FoldChange",  "padj",       "symbol",
    "GSE185009 NPC",             "iPSC",        "deseq_GSE185009_npc",         "df_vsd_GSE185009_npc",          "log2FoldChange",  "padj",       "symbol",
    "GSE196190 t24 100uM",       "iPSC",        "deseq_GSE196190_t24_100",     "df_vsd_GSE196190_t24_100qM",    "log2FoldChange",  "padj",       "symbol",
    "GSE196190 t24 200uM",       "iPSC",        "deseq_GSE196190_t24_200",     "df_vsd_GSE196190_t24_200qM",    "log2FoldChange",  "padj",       "symbol",
    "GSE196190 t24 400uM",       "iPSC",        "deseq_GSE196190_t24_400",     "df_vsd_GSE196190_t24_400qM",    "log2FoldChange",  "padj",       "symbol",
    "GSE196190 t74 100uM",       "iPSC",        "deseq_GSE196190_t74_100",     "df_vsd_GSE196190_t74_100qM",    "log2FoldChange",  "padj",       "symbol",
    "GSE196190 t74 200uM",       "iPSC",        "deseq_GSE196190_t74_200",     "df_vsd_GSE196190_t74_200qM",    "log2FoldChange",  "padj",       "symbol",
    "GSE196190 t74 400uM",       "iPSC",        "deseq_GSE196190_t74_400",     "df_vsd_GSE196190_t74_400qM",    "log2FoldChange",  "padj",       "symbol",
    "GSE120746",                 "iPSC",        "deseq_GSE120746",             "df_vsd_GSE120746",              "log2FoldChange",  "padj",       "symbol",
    "GSE282494",                 "iPSC",        "deseq_GSE282494",             "df_vsd_GSE282494",              "log2FoldChange",  "padj",       "symbol",
    "GSE181029 NP",              "iPSC",        "deseq_GSE181029_NP",          "df_vsd_GSE181029_NP",           "log2FoldChange",  "padj",       "symbol",
    "GSE181029 DN",              "iPSC",        "deseq_GSE181029_DN",          "df_vsd_GSE181029_DN",           "log2FoldChange",  "padj",       "symbol",
    "GSE315738",                 "iPSC",        "deseq_GSE315738",             "df_vsd_GSE315738",              "log2FoldChange",  "padj",       "symbol",
    "GSE36321",                  "iPSC",        "deseq_GSE36321",              "df_vsd_GSE36321",               "logFC",            "adj.P.Val",  "symbol",
    "GSE51922 idiopathic",       "iPSC",        "deseq_GSE51922",              "df_vsd_GSE51922",               "logFC",            "adj.P.Val",  "symbol",
    "GSE51922 LRRK2",            "iPSC",        "deseq_GSE51922_LRRK2",        "df_vsd_GSE51922_LRRK2",         "logFC",            "adj.P.Val",  "symbol",
    
    "GSE285507",                 "Cell-based",  "deseq_GSE285507",             "df_vsd_GSE285507",              "log2FoldChange",  "padj",       "symbol",
    "GSE229460",                 "Cell-based",  "deseq_GSE229460",             "df_vsd_GSE229460",              "log2FoldChange",  "padj",       "symbol",
    "GSE203522 2D",              "Cell-based",  "deseq_GSE203522_2D",          "df_vsd_GSE203522_2D",           "log2FoldChange",  "padj",       "symbol",
    "GSE203522 3D",              "Cell-based",  "deseq_GSE203522_3D",          "df_vsd_GSE203522_3D",           "log2FoldChange",  "padj",       "symbol",
    "GSE315111 Rotenone",        "Cell-based",  "deseq_GSE315111_rot",         "df_vsd_GSE315111_rot",          "log2FoldChange",  "padj",       "symbol",
    "GSE315111 MPP+",            "Cell-based",  "deseq_GSE315111_MPP",         "df_vsd_GSE315111_MPP",          "log2FoldChange",  "padj",       "symbol",
    "GSE21305",                  "Cell-based",  "deseq_GSE21305",              "df_vsd_GSE21305",               "logFC",            "adj.P.Val",  "GeneSymbol",
    "GSE191302",                 "Cell-based",  "deseq_GSE191302",             "df_vsd_GSE191302",              "logFC",            "adj.P.Val",  "symbol",
    "GSE35642 5nM",              "Cell-based",  "deseq_GSE35642_5mM",          "df_vsd_GSE35642_5mM",           "logFC",            "adj.P.Val",  "symbol",
    "GSE35642 50nM",             "Cell-based",  "deseq_GSE35642_50mM",         "df_vsd_GSE35642_50mM",          "logFC",            "adj.P.Val",  "symbol",
    "GSE287941 6-OHDA",          "Cell-based",  "deseq_GSE287941_6OHDA",       "df_vsd_GSE287941_6OHDA",        "log2FoldChange",  "padj",       "SYMBOL",
    "GSE287941 6-OHDA+asyn",     "Cell-based",  "deseq_GSE287941_6OHDA_asyn",  "df_vsd_GSE287941_6OHDA_asyn",   "log2FoldChange",  "padj",       "SYMBOL"
  )
  stopifnot(nrow(ds_registry) == 33)
  
  get_gene_row <- function(df, symbol_col, gene, fc_col, p_col) {
    if (symbol_col %in% colnames(df)) {
      hit <- df[toupper(trimws(as.character(df[[symbol_col]]))) == gene, ]
    } else {
      hit <- df[toupper(trimws(rownames(df))) == gene, ]
    }
    if (nrow(hit) == 0) return(list(log2FC = NA_real_, padj = NA_real_))
    if (nrow(hit) > 1) hit <- hit[which.min(hit[[p_col]]), ]
    list(log2FC = as.numeric(hit[[fc_col]]), padj = as.numeric(hit[[p_col]]))
  }
  
  get_spearman <- function(vsd_df, symbol_col, dataset_name = NA_character_) {
    if (symbol_col %in% colnames(vsd_df)) {
      sub <- vsd_df[toupper(trimws(as.character(vsd_df[[symbol_col]]))) %in% c("SIRT1", "SNCA"), ]
      rownames(sub) <- toupper(trimws(as.character(sub[[symbol_col]])))
      sub <- sub[, setdiff(colnames(sub), symbol_col), drop = FALSE]
    } else {
      sub <- vsd_df[toupper(trimws(rownames(vsd_df))) %in% c("SIRT1", "SNCA"), ]
      rownames(sub) <- toupper(trimws(rownames(sub)))
    }
    if (!all(c("SIRT1", "SNCA") %in% rownames(sub))) {
      return(list(rho = NA_real_, pval = NA_real_, n = NA_integer_, dropped_cols = character(0)))
    }
    
    is_num_col <- vapply(sub, function(col) {
      suppressWarnings(!any(is.na(as.numeric(as.character(col))) & !is.na(col)))
    }, logical(1))
    dropped <- colnames(sub)[!is_num_col]
    if (length(dropped) > 0) {
      cat("!! [", dataset_name, "] ",
          paste(dropped, collapse = ", "), "\n", sep = "")
    }
    sub <- sub[, is_num_col, drop = FALSE]
    
    x <- as.numeric(sub["SIRT1", ]); y <- as.numeric(sub["SNCA", ])
    ok <- !is.na(x) & !is.na(y)
    if (sum(ok) < 3) {
      return(list(rho = NA_real_, pval = NA_real_, n = sum(ok), dropped_cols = dropped))
    }
    test <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
    list(rho = round(as.numeric(test$estimate), 3), pval = round(test$p.value, 4),
         n = sum(ok), dropped_cols = dropped)
  }
  
  results <- lapply(seq_len(nrow(ds_registry)), function(i) {
    row <- ds_registry[i, ]
    note <- character(0)
    if (!exists(row$deseq_obj)) note <- c(note, paste("deseq not found", row$deseq_obj))
    if (!exists(row$vsd_obj))   note <- c(note, paste("vsd not found", row$vsd_obj))
    
    if (length(note) > 0) {
      return(tibble::tibble(dataset_clean = row$dataset_clean, model_group = row$model_group,
                            sirt1_log2FC = NA, sirt1_padj = NA,
                            snca_log2FC = NA, snca_padj = NA,
                            rho = NA, pval = NA, n = NA,
                            note = paste(note, collapse = "; ")))
    }
    
    deseq_df <- get(row$deseq_obj)
    vsd_df   <- get(row$vsd_obj)
    
    sirt1 <- get_gene_row(deseq_df, row$symbol_col, "SIRT1", row$fc_col, row$p_col)
    snca  <- get_gene_row(deseq_df, row$symbol_col, "SNCA",  row$fc_col, row$p_col)
    sp    <- get_spearman(vsd_df, row$symbol_col, dataset_name = row$dataset_clean)
    
    sp_note <- if (length(sp$dropped_cols) > 0) {
      paste0("Atlanan sutun(lar): ", paste(sp$dropped_cols, collapse = ", "))
    } else NA_character_
    
    tibble::tibble(
      dataset_clean = row$dataset_clean, model_group = row$model_group,
      sirt1_log2FC = sirt1$log2FC, sirt1_padj = sirt1$padj,
      snca_log2FC  = snca$log2FC,  snca_padj  = snca$padj,
      rho = sp$rho, pval = sp$pval, n = sp$n,
      note = sp_note
    )
  })
  
  pattern_data <- bind_rows(results) %>%
    mutate(
      sirt1_sig = !is.na(sirt1_padj) & sirt1_padj < SIG_PADJ & abs(sirt1_log2FC) > SIG_LFC,
      snca_sig  = !is.na(snca_padj)  & snca_padj  < SIG_PADJ & abs(snca_log2FC)  > SIG_LFC,
      sirt1_dir = ifelse(is.na(sirt1_log2FC), NA_character_, ifelse(sirt1_log2FC > 0, "up", "down")),
      snca_dir  = ifelse(is.na(snca_log2FC),  NA_character_, ifelse(snca_log2FC  > 0, "up", "down"))
    )
  
  classify_pattern <- function(sirt1_sig, sirt1_dir, snca_sig, snca_dir, rho) {
    if (is.na(rho) || (!isTRUE(sirt1_sig) && !isTRUE(snca_sig))) return(list(pattern = "Unresolved", conf = "low"))
    if (abs(rho) < RHO_MIN) return(list(pattern = "Unresolved", conf = "low"))
    
    if (isTRUE(sirt1_sig)) {
      if (rho < 0 && sirt1_dir == "up")   return(list(pattern = "Compensatory",           conf = "high"))
      if (rho < 0 && sirt1_dir == "down") return(list(pattern = "Exhausted Compensation", conf = "high"))
      if (rho > 0 && sirt1_dir == "up")   return(list(pattern = "Acute Co-induction",     conf = "high"))
      return(list(pattern = "Unresolved", conf = "low"))
    }
    
    if (isTRUE(snca_sig)) {
      if (rho < 0 && snca_dir == "down") return(list(pattern = "Compensatory",            conf = "medium"))
      if (rho < 0 && snca_dir == "up")   return(list(pattern = "Exhausted Compensation",  conf = "medium"))
      if (rho > 0 && snca_dir == "up")   return(list(pattern = "Acute Co-induction",      conf = "medium"))
      return(list(pattern = "Unresolved", conf = "low"))
    }
    
    list(pattern = "Unresolved", conf = "low")
  }
  
  pattern_results <- lapply(seq_len(nrow(pattern_data)), function(i) {
    r <- pattern_data[i, ]
    cl <- classify_pattern(r$sirt1_sig, r$sirt1_dir, r$snca_sig, r$snca_dir, r$rho)
    tibble::tibble(pattern = cl$pattern, confidence = cl$conf)
  })
  pattern_data <- bind_cols(pattern_data, bind_rows(pattern_results)) %>%
    mutate(
      ## star assigned based on rho's OWN p-value (cor.test pval) --- identical to the old
      ## threshold: p<0.001=***, p<0.01=**, p<0.05=*, otherwise blank
      rho_sig_stars = as.character(cut(
        pval,
        breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
        labels = c("***", "**", "*", "")
      )),
      rho_sig_stars = ifelse(is.na(rho_sig_stars), "", rho_sig_stars),
      rho_label = ifelse(is.na(rho), NA_character_, paste0(sprintf("%.2f", rho), rho_sig_stars)),
      n_label   = ifelse(is.na(n), "n=NA", paste0("n=", n)),
      small_n_flag = !is.na(n) & n < 6  ## small-n warning --- interpret rho with caution
    )

  print(pattern_data %>%
          select(dataset_clean, model_group, pattern, confidence, rho, pval, n,
                 sirt1_dir, sirt1_sig, snca_dir, snca_sig, note) %>%
          arrange(pattern, dataset_clean), n = Inf)
  
  write.csv(pattern_data, "SIRT1_SNCA_pattern_classification.csv", row.names = FALSE)

  
  ## =============================================================================
  # Section 3: build the final figure, grouping datasets by their assigned SIRT1-SNCA pattern.
  ## SECTION 3 --- FINAL FIGURE (Figure 3, grouped by pattern)
  ## =============================================================================
  
  pattern_levels <- c("Compensatory",
                      "Exhausted Compensation",
                      "Acute Co-induction",
                      "Unresolved")
  
  plot_data <- pattern_data %>%
    filter(!is.na(pattern)) %>%
    mutate(pattern = factor(pattern, levels = pattern_levels)) %>%
    arrange(pattern, dataset_clean) %>%
    mutate(dataset_clean = factor(dataset_clean, levels = unique(dataset_clean)),
           model_group = factor(model_group, levels = c("Postmortem", "iPSC", "Cell-based")))
  
  pattern_bounds <- plot_data %>%
    mutate(x = as.integer(dataset_clean)) %>%
    group_by(pattern) %>%
    summarise(xmin = min(x) - 0.5, xmax = max(x) + 0.5, xmid = mean(x), .groups = "drop")
  
  group_colors <- c("Postmortem" = "#1E8449", "iPSC" = "#D4AC0D", "Cell-based" = "#6C3483")
  pattern_bg   <- c("Compensatory" = "#EAF2E3",
                    "Exhausted Compensation" = "#FDEBD0",
                    "Acute Co-induction" = "#FADBD8",
                    "Unresolved" = "#EAEAEA")
  
  cap_fc <- max(2, quantile(abs(c(plot_data$sirt1_log2FC, plot_data$snca_log2FC)), 0.90, na.rm = TRUE))
  
  ## only SIGNIFICANT SIRT1/SNCA get colored (red=up, blue=down);
  ## if not significant, grey --- NOT a gradient, categorical color based on significance
  plot_data <- plot_data %>%
    mutate(
      sirt1_status = case_when(
        !is.na(sirt1_sig) & sirt1_sig & sirt1_log2FC > 0 ~ "Up (significant)",
        !is.na(sirt1_sig) & sirt1_sig & sirt1_log2FC < 0 ~ "Down (significant)",
        TRUE                                              ~ "Not significant"
      ),
      snca_status = case_when(
        !is.na(snca_sig) & snca_sig & snca_log2FC > 0 ~ "Up (significant)",
        !is.na(snca_sig) & snca_sig & snca_log2FC < 0 ~ "Down (significant)",
        TRUE                                            ~ "Not significant"
      )
    )
  
  cat("SIRT1 status dagilimi:\n"); print(table(plot_data$sirt1_status))
  cat("SNCA status dagilimi:\n");  print(table(plot_data$snca_status))
  
  FC_STATUS_COLORS <- c("Up (significant)" = "#E60000",
                        "Down (significant)" = "#0033CC",
                        "Not significant" = "#CFCFCF")
  
  ## compute an upper limit the pattern labels won't overflow (to fix the clipping issue)
  fc_range <- range(c(plot_data$sirt1_log2FC, plot_data$snca_log2FC), na.rm = TRUE)
  label_y_fc <- fc_range[2] + diff(fc_range) * 0.15
  ylim_top   <- label_y_fc + diff(fc_range) * 0.08
  ## ---------------------------------------------------------------------------
  ## Only "Exhausted Compensation" is split into two lines, the others stay as-is
  ## ---------------------------------------------------------------------------
  ## ---------------------------------------------------------------------------
  ## Only "Exhausted Compensation" is split into two lines, the others stay as-is
  ## ---------------------------------------------------------------------------
  pattern_bounds$pattern <- recode(pattern_bounds$pattern,
                                   "Exhausted Compensation" = "Exhausted\nCompensation")
  
  p_fc <- ggplot(plot_data, aes(x = dataset_clean)) +
    geom_rect(data = pattern_bounds, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = NULL),
              ymin = -Inf, ymax = Inf, fill = pattern_bg[pattern_bounds$pattern], alpha = 0.5, inherit.aes = FALSE) +
    geom_hline(yintercept = 0, color = "#666666", linetype = "dashed", linewidth = 0.5) +
    geom_col(aes(y = sirt1_log2FC, fill = sirt1_status), width = 0.45, color = "black", linewidth = 0.3) +
    geom_line(aes(y = snca_log2FC, group = 1), color = "#888888", linewidth = 0.9, alpha = 0.6) +
    geom_point(aes(y = snca_log2FC, fill = snca_status), shape = 23, size = 4.2, color = "black", stroke = 0.4) +
    scale_fill_manual(values = FC_STATUS_COLORS, name = "SIRT1 / SNCA\n(bar / diamond)") +
    geom_vline(data = pattern_bounds %>% filter(xmin > 0.5), aes(xintercept = xmin), color = "black", linewidth = 1) +
    ## pattern labels --- white-boxed geom_label, "Exhausted Compensation" on two lines
    geom_label(data = pattern_bounds, aes(x = xmid, y = label_y_fc, label = pattern),
               fontface = "bold", size = 4.2, color = "black", fill = "white",
               label.size = 0.3, label.padding = unit(0.25, "lines"), inherit.aes = FALSE,
               lineheight = 0.85) +
    coord_cartesian(ylim = c(fc_range[1] - diff(fc_range) * 0.05, ylim_top), clip = "off") +
    labs(title = "",
         subtitle = "",
         x = NULL, y = "log2 Fold Change") +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 12, color = "black"),
      plot.subtitle = element_text(size = 11, color = "grey30"),
      axis.text.x = element_blank(), axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 12, face = "bold", color = "black"),
      axis.title.y = element_text(size = 12, face = "bold"),
      panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "grey92", linewidth = 0.4),
      legend.position = "right",
      plot.margin = margin(30, 15, 5, 15)
    )
  
  ## ---------------------------------------------------------------------------
  ## RHO PANEL --- font size of the rho boxes reduced (4.2 -> 2.6), to reduce
  ## overlap. label.padding also reduced so the box shrinks proportionally.
  ## ---------------------------------------------------------------------------
  plot_data <- plot_data %>%
    mutate(rho_status = ifelse(rho >= 0, "Positive", "Negative"),
           arrow_shape = ifelse(rho >= 0, 24, 25))
  RHO_COLORS <- c("Positive" = "#E60000", "Negative" = "#0033CC")
  
  p_rho <- ggplot(plot_data, aes(x = dataset_clean)) +
    geom_rect(data = pattern_bounds, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
              ymin = -Inf, ymax = Inf, fill = pattern_bg[pattern_bounds$pattern], alpha = 0.5, inherit.aes = FALSE) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = 1, fill = "#E60000", alpha = 0.05) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -1, ymax = 0, fill = "#0033CC", alpha = 0.05) +
    geom_hline(yintercept = c(-1, -0.5, 0.5, 1), color = "grey88", linewidth = 0.4) +
    geom_hline(yintercept = 0, color = "#666666", linewidth = 0.7, linetype = "dashed") +
    geom_segment(aes(xend = dataset_clean, y = 0, yend = rho, color = rho_status),
                 linewidth = 5, lineend = "round", show.legend = FALSE) +
    geom_point(aes(y = rho, shape = arrow_shape, fill = rho_status), size = 7, stroke = 0, show.legend = FALSE) +
    geom_label(aes(y = 1.0, label = n_label), size = 3.9, color = "grey25", fill = "white",
               label.size = 0, label.padding = unit(0.12, "lines"), fontface = "bold") +
    ## rho box --- font size reduced
    geom_label(aes(label = rho_label, color = rho_status), y = 0, vjust = 0.5, size = 2.6,
               fontface = "bold", fill = "white", label.size = 0,
               label.padding = unit(0.05, "lines"), show.legend = FALSE) +
    scale_shape_identity() +
    scale_color_manual(values = RHO_COLORS) +
    scale_fill_manual(values = RHO_COLORS) +
    geom_vline(data = pattern_bounds %>% filter(xmin > 0.5), aes(xintercept = xmin), color = "black", linewidth = 1) +
    scale_y_continuous(limits = c(-1, 1), breaks = c(-1, -0.5, 0, 0.5, 1),
                       labels = c("-1", "-0.5", "0", "+0.5", "+1"),
                       expand = expansion(mult = c(0.08, 0.12))) +
    labs(x = NULL, y = "Spearman rho", caption = "* p<0.05  ** p<0.01  *** p<0.001 (cor.test)") +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5, face = "bold", color = "black"),
      axis.text.y = element_text(size = 12, face = "bold", color = "black"),
      axis.title.y = element_text(size = 12, face = "bold"),
      panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
      panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5),
      plot.caption = element_text(size = 10, color = "grey50")
    )
  
  p_strip <- ggplot(plot_data, aes(x = dataset_clean, y = 1, fill = model_group)) +
    geom_tile(color = "white") +
    scale_fill_manual(values = group_colors, name = "Model group") +
    theme_void() + theme(legend.position = "bottom")
  
  layout <- p_fc / p_rho / p_strip +
    plot_layout(heights = c(3.2, 2.2, 0.35), guides = "collect") &
    theme(legend.position = "right")
  
  print(layout)
  
  
  
  
  print(layout)
  ggsave("Figure2_mol_omic_SIRT1_SNCA_patterns.tiff", layout,
         width =40, height = 25, units = "cm", dpi = 600, compression = "lzw")
  ## SAVING --- width MUST be in cm (in the previous version units="in" 
  
  
  
  output_dir  <- "figures"
  output_name <- "Figure3_SIRT1_SNCA_pattern_grouped"
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  ggsave(file.path(output_dir, paste0(output_name, ".tiff")),
         plot = layout, width = 20, height = 15, dpi = 600, device = "tiff", bg = "white")
  ggsave(file.path(output_dir, paste0(output_name, ".png")),
         plot = layout, width = 20, height = 15, dpi = 300, bg = "white")
  
  cat("\nKaydedildi:\n",
      " - ", file.path(output_dir, paste0(output_name, ".tiff")), " (600 dpi)\n",
      " - ", file.path(output_dir, paste0(output_name, ".png")),  " (300 dpi, on izleme)\n", sep = "")
}
# FIGURE 3: chord diagram linking resveratrol, SIRT1, SNCA, and datasets by regulatory direction, built with the circlize package.
#FIGURE 3 SCRIPT
{
  # ==============================================================================
  #  RESVERATROL - SIRT1 - SNCA CHORD DIAGRAM
  #  Full code - data preparation + drawing
  # ==============================================================================
  
  library(readxl)
  library(dplyr)
  library(tibble)
  library(circlize)
  # Resolve all filter conflicts
  conflicts_prefer(dplyr::filter)
  conflicts_prefer(dplyr::select)
  conflicts_prefer(dplyr::rename)
  conflicts_prefer(base::intersect)
  
  conflicts_prefer(base::intersect)
  
  # ==============================================================================
  # Step 1: read all per-dataset DESeq2/limma differential expression result files.
  # 1. READ DESEQ DATA
  # ==============================================================================
  
  deseq_GSE156928          <- read.csv("deseq_GSE156928.csv",             row.names = 1)
  deseq_GSE114918_SNc      <- read.csv("deseq_GSE114918_SNc.csv",         row.names = 1)
  deseq_GSE114918_VTA      <- read.csv("deseq_GSE114918_VTA.csv",         row.names = 1)
  deseq_GSE20292           <- read.csv("deseq_GSE20292.csv",               row.names = 1)
  deseq_GSE20291           <- read.csv("deseq_GSE20291.csv",               row.names = 1)
  deseq_GSE185009_ipsc     <- read.csv("deseq_GSE185009_ipsc.csv",        row.names = 1)
  deseq_GSE185009_npc      <- read.csv("deseq_GSE185009_npc.csv",         row.names = 1)
  deseq_GSE196190_t24_100  <- read.csv("deseq_GSE196190_t24_100qM.csv",  row.names = 1)
  deseq_GSE196190_t24_200  <- read.csv("deseq_GSE196190_t24_200qM.csv",  row.names = 1)
  deseq_GSE196190_t24_400  <- read.csv("deseq_GSE196190_t24_400qM.csv",  row.names = 1)
  deseq_GSE196190_t74_100  <- read.csv("deseq_GSE196190_t74_100qM.csv",  row.names = 1)
  deseq_GSE196190_t74_200  <- read.csv("deseq_GSE196190_t74_200qM.csv",  row.names = 1)
  deseq_GSE196190_t74_400  <- read.csv("deseq_GSE196190_t74_400qM.csv",  row.names = 1)
  deseq_GSE120746          <- read.csv("deseq_GSE120746.csv",              row.names = 1)
  deseq_GSE282494          <- read.csv("deseq_GSE282494.csv",              row.names = 1)
  deseq_GSE181029_NP       <- read.csv("deseq_GSE181029_NP.csv",          row.names = 1)
  deseq_GSE181029_DN       <- read.csv("deseq_GSE181029_DN.csv",          row.names = 1)
  deseq_GSE315738          <- read.csv("deseq_GSE315738.csv",              row.names = 1)
  deseq_GSE36321           <- read.csv("deseq_GSE36321.csv",               row.names = 1)
  deseq_GSE51922           <- read.csv("deseq_GSE51922.csv",               row.names = 1)
  deseq_GSE51922_LRRK2    <- read.csv("deseq_GSE51922_LRRK2.csv",        row.names = 1)
  deseq_GSE285507          <- read.csv("deseq_GSE285507.csv",              row.names = 1)
  deseq_GSE229460          <- read.csv("deseq_GSE229460.csv",              row.names = 1)
  deseq_GSE203522_2D       <- read.csv("deseq_GSE203522_2D.csv",          row.names = 1)
  deseq_GSE203522_3D       <- read.csv("deseq_GSE203522_3D.csv",          row.names = 1)
  deseq_GSE315111_rot      <- read.csv("deseq_GSE315111_rot.csv",         row.names = 1)
  deseq_GSE315111_MPP      <- read.csv("deseq_GSE315111_MPP.csv",         row.names = 1)
  deseq_GSE21305           <- read.csv("deseq_GSE21305.csv",               row.names = 1)
  deseq_GSE191302          <- read.csv("deseq_GSE191302.csv",              row.names = 1)
  deseq_GSE35642_5mM       <- read.csv("deseq_GSE35642_5mM.csv",          row.names = 1)
  deseq_GSE35642_50mM      <- read.csv("deseq_GSE35642_50mM.csv",         row.names = 1)
  deseq_GSE287941_6OHDA    <- read.csv("deseq_GSE287941_6OHDA.csv",       row.names = 1)
  deseq_GSE287941_6OHDA_asyn    <- read.csv("deseq_GSE287941_6OHDA+asyn.csv", row.names = 1)
  
  # ==============================================================================
  # Step 2: read the resveratrol drug-gene interaction data and clean stray line-break characters from text fields.
  # 2. READ RESVERATROL DATA AND CLEAN \r\n
  # ==============================================================================
  file_post  <- "C:/Users/User/OneDrive/Belgeler/ana_figure_postmortem_03-06/All_Approved_Drug_Interactions_Combined_postmortem_ALL_DATA_supp_3_06_2026_LAST_supp.xlsx"
  file_ipsc <- "C:/Users/User/OneDrive/Belgeler/ana_figure_ipsc_03.06/All_Approved_Drug_Interactions_Combined__IPSC_ALL_DATA_supp_3_06_2026_LAST_supp.xlsx"
  file_cell <- "C:/Users/User/OneDrive/Belgeler/ana_figure_cell_based_03-06/All_Approved_Drug_Interactions_Combined__Cell_based_ALL_DATA_supp_3_06_2026_LAST_supp.xlsx"
  
  cell_df <- read_excel(file_cell)
  ipsc_df <- read_excel(file_ipsc)
  post_df <- read_excel(file_post)
  
  
  cell_resv <- cell_df %>% dplyr::filter(drug == "RESVERATROL")
  ipsc_resv <- ipsc_df %>% dplyr::filter(drug == "RESVERATROL")
  post_resv <- post_df %>% dplyr::filter(drug == "RESVERATROL")
  
  cat("Satir sayilari:", nrow(cell_resv), nrow(ipsc_resv), nrow(post_resv), "\n")
  
  all_resveratrol <- bind_rows(
    cell_based = cell_resv,
    ipsc       = ipsc_resv,
    postmortem = post_resv,
    .id = "Group"
  ) %>%
    select(-any_of(c("regulatory approval", "indication", "regulatory.approval"))) %>%
    mutate(
      dataset = gsub("\r\n", " ", dataset),
      dataset = gsub("\r",   " ", dataset),
      dataset = gsub("\n",   " ", dataset),
      dataset = trimws(dataset)
    )
  

  print(head(unique(all_resveratrol$dataset), 10))
  
  
  all_genes_clean <- unique(all_resveratrol$gene)
  all_genes_clean <- all_genes_clean[!all_genes_clean %in% c("SIRT1", "SNCA")]
  
  # ==============================================================================
  # Step 3: build the full list of 33 datasets and map each to its model group.
  # 3. DATASET LIST AND GROUP MAP
  # ==============================================================================
  
  datasets <- list(
    "GSE156928 Frontal Cortex"          = deseq_GSE156928,
    "GSE114918 SNc"                     = deseq_GSE114918_SNc,
    "GSE114918 VTA"                     = deseq_GSE114918_VTA,
    "GSE20292 SN"                       = deseq_GSE20292,
    "GSE20291 Putamen"                  = deseq_GSE20291,
    "GSE185009 iPSC"                    = deseq_GSE185009_ipsc,
    "GSE185009 Neural progenitor cells" = deseq_GSE185009_npc,
    "GSE196190 t24 100 uM"              = deseq_GSE196190_t24_100,
    "GSE196190 t24 200 uM"              = deseq_GSE196190_t24_200,
    "GSE196190 t24 400 uM"              = deseq_GSE196190_t24_400,
    "GSE196190 t74 100 uM"              = deseq_GSE196190_t74_100,
    "GSE196190 t74 200 uM"              = deseq_GSE196190_t74_200,
    "GSE196190 t74 400 uM"              = deseq_GSE196190_t74_400,
    "GSE120746 Dopaminergic neuron"     = deseq_GSE120746,
    "GSE282494 LRRK2 mutation"          = deseq_GSE282494,
    "GSE181029 Neural progenitor"       = deseq_GSE181029_NP,
    "GSE181029 Dopaminergic neuron"     = deseq_GSE181029_DN,
    "GSE315738 GBA"                     = deseq_GSE315738,
    "GSE36321 Neural stem cell"         = deseq_GSE36321,
    "GSE51922 idiopatic"                = deseq_GSE51922,
    "GSE51922 LRRK2 mutation"           = deseq_GSE51922_LRRK2,
    "GSE285507 LUHMES Pb"               = deseq_GSE285507,
    "GSE229460 LUHMES MPP+"             = deseq_GSE229460,
    "GSE203522 2D SH-SY5Y MPP+"         = deseq_GSE203522_2D,
    "GSE203522 3D SH-SY5Y MPP+"         = deseq_GSE203522_3D,
    "GSE315111 SH-SY5Y rotenone"        = deseq_GSE315111_rot,
    "GSE315111 SH-SY5Y MPP+"            = deseq_GSE315111_MPP,
    "GSE21305 SH-SY5Y paraquat"         = deseq_GSE21305,
    "GSE191302 LUHMES asyn"             = deseq_GSE191302,
    "GSE35642 SK-N-MC rotenone 5 mM"    = deseq_GSE35642_5mM,
    "GSE35642 SK-N-MC rotenone 50 mM"   = deseq_GSE35642_50mM,
    "GSE287941 LUHMES 6OHDA"            = deseq_GSE287941_6OHDA,
    "GSE287941 LUHMES 6OHDA+asyn"       = deseq_GSE287941_6OHDA_asyn
  )
  
  dataset_group_map <- c(
    "GSE156928 Frontal Cortex"          = "Postmortem",
    "GSE114918 SNc"                     = "Postmortem",
    "GSE114918 VTA"                     = "Postmortem",
    "GSE20292 SN"                       = "Postmortem",
    "GSE20291 Putamen"                  = "Postmortem",
    "GSE185009 iPSC"                    = "iPSC",
    "GSE185009 Neural progenitor cells" = "iPSC",
    "GSE196190 t24 100 uM"              = "iPSC",
    "GSE196190 t24 200 uM"              = "iPSC",
    "GSE196190 t24 400 uM"              = "iPSC",
    "GSE196190 t74 100 uM"              = "iPSC",
    "GSE196190 t74 200 uM"              = "iPSC",
    "GSE196190 t74 400 uM"              = "iPSC",
    "GSE120746 Dopaminergic neuron"     = "iPSC",
    "GSE282494 LRRK2 mutation"          = "iPSC",
    "GSE181029 Neural progenitor"       = "iPSC",
    "GSE181029 Dopaminergic neuron"     = "iPSC",
    "GSE315738 GBA"                     = "iPSC",
    "GSE36321 Neural stem cell"         = "iPSC",
    "GSE51922 idiopatic"                = "iPSC",
    "GSE51922 LRRK2 mutation"           = "iPSC",
    "GSE285507 LUHMES Pb"               = "Cell-based",
    "GSE229460 LUHMES MPP+"             = "Cell-based",
    "GSE203522 2D SH-SY5Y MPP+"         = "Cell-based",
    "GSE203522 3D SH-SY5Y MPP+"         = "Cell-based",
    "GSE315111 SH-SY5Y rotenone"        = "Cell-based",
    "GSE315111 SH-SY5Y MPP+"            = "Cell-based",
    "GSE21305 SH-SY5Y paraquat"         = "Cell-based",
    "GSE191302 LUHMES asyn"             = "Cell-based",
    "GSE35642 SK-N-MC rotenone 5 mM"    = "Cell-based",
    "GSE35642 SK-N-MC rotenone 50 mM"   = "Cell-based",
    "GSE287941 LUHMES 6OHDA"            = "Cell-based",
    "GSE287941 LUHMES 6OHDA+asyn"       = "Cell-based"
  )
  
  dataset_order <- names(dataset_group_map)
  
  group_colors <- c(
    "Cell-based" = "#6C3483",
    "iPSC"       = "#D4AC0D",
    "Postmortem" = "#1E8449"
  )
  # ==============================================================================
  # Step 4: extract SIRT1 and SNCA expression/statistics values for each dataset.
  # 4. EXTRACT SIRT1 / SNCA VALUES
  # ==============================================================================
  
  extract_genes <- function(df, gse_name, genes = c("SIRT1", "SNCA")) {
    if (!"gene" %in% colnames(df)) {
      df <- df %>% rownames_to_column("gene")
    }
    if ("symbol" %in% colnames(df)) {
      df <- df %>% mutate(gene = symbol)
    }
    logfc_col <- base::intersect(c("logFC", "log2FoldChange"), colnames(df))[1]
    padj_col  <- base::intersect(c("adj.P.Val", "padj"),       colnames(df))[1]
    if (is.na(logfc_col) || is.na(padj_col)) {
      warning(paste("Eksik sutun:", gse_name))
      return(NULL)
    }
    df %>%
      filter(gene %in% genes) %>%
      transmute(
        gene,
        logFC     = .data[[logfc_col]],
        adj.P.Val = .data[[padj_col]],
        GSE       = gse_name
      )
  }
  
  gene_expr_list <- mapply(
    FUN      = extract_genes,
    df       = datasets,
    gse_name = names(datasets),
    SIMPLIFY = FALSE
  )
  
  gene_expr_df <- bind_rows(gene_expr_list)
  
  print(table(gene_expr_df$gene))
  
  # ==============================================================================
  # Step 5: compute the regulation status (up/down/not significant) for each gene-dataset combination.
  # 5. COMPUTE STATUS
  # ==============================================================================
  
  gene_expr_df <- gene_expr_df %>%
    mutate(
      status = case_when(
        adj.P.Val < 0.05 & logFC >= 0.6  ~ "UP_sig",
        adj.P.Val < 0.05 & logFC <= -0.6 ~ "DOWN_sig",
        logFC > 0                          ~ "UP_ns",
        logFC < 0                          ~ "DOWN_ns",
        TRUE                               ~ "NC"
      )
    )
  
  sirt1_df <- gene_expr_df %>%
    filter(gene == "SIRT1") %>%
    rename(logFC_SIRT1 = logFC, adj.P.Val_SIRT1 = adj.P.Val, SIRT1_status = status) %>%
    select(GSE, logFC_SIRT1, adj.P.Val_SIRT1, SIRT1_status)
  
  snca_df <- gene_expr_df %>%
    filter(gene == "SNCA") %>%
    rename(logFC_SNCA = logFC, adj.P.Val_SNCA = adj.P.Val, SNCA_status = status) %>%
    select(GSE, logFC_SNCA, adj.P.Val_SNCA, SNCA_status)
  
  # ==============================================================================
  # Step 6: define the color-mapping function used for chord diagram links and sectors.
  # 6. COLOR FUNCTION
  # ==============================================================================
  
  status_color <- function(status) {
    dplyr::case_when(
      status == "UP_sig"   ~ "#C0392BDD",
      status == "DOWN_sig" ~ "#1A5276DD",
      status == "UP_ns"    ~ "#F1948A55",
      status == "DOWN_ns"  ~ "#85C1E955",
      TRUE                 ~ "#AAAAAA33"
    )
  }
  
  # ==============================================================================
  # Step 7: build the link tables (source-target pairs and weights) that define the chord diagram's ribbons.
  # 7. LINK TABLES
  # the dataset column is already cleaned, use it directly
  # ==============================================================================
  
  df <- all_resveratrol %>%
    mutate(gse_label = dataset)
  
  cat("Eslesen sayi:\n")
  print(sum(unique(df$gse_label) %in% dataset_order))
  
  cat("Eslesmeyen dataset'ler:\n")
  print(unique(df$gse_label)[!unique(df$gse_label) %in% dataset_order])
  
  link_drug_to_gene <- df %>%
    filter(!(gene %in% c("SIRT1", "SNCA"))) %>%
    group_by(from = "RESVERATROL", to = gene) %>%
    summarise(value = sum(abs(`interaction score`)) * 100, .groups = "drop") %>%
    mutate(color = "#E67E22A6")
  
  link_dataset_to_gene <- df %>%
    filter(!(gene %in% c("SIRT1", "SNCA"))) %>%
    filter(gse_label %in% dataset_order) %>%
    group_by(from = gse_label, to = gene) %>%
    summarise(value = sum(abs(`interaction score`)) * 100, .groups = "drop") %>%
    mutate(color = "#E8A0A055")
  
  link_sirt1 <- sirt1_df %>%
    filter(GSE %in% dataset_order) %>%
    group_by(from = "SIRT1", to = GSE) %>%
    summarise(
      value  = mean(abs(logFC_SIRT1)) * 50 + 1,
      status = case_when(
        any(adj.P.Val_SIRT1 < 0.05 & logFC_SIRT1 >= 0.6)  ~ "UP_sig",
        any(adj.P.Val_SIRT1 < 0.05 & logFC_SIRT1 <= -0.6) ~ "DOWN_sig",
        mean(logFC_SIRT1) >= 0                              ~ "UP_ns",
        TRUE                                                ~ "DOWN_ns"
      ),
      .groups = "drop"
    ) %>%
    mutate(color = status_color(status)) %>%
    select(from, to, value, color)
  
  link_snca <- snca_df %>%
    filter(GSE %in% dataset_order) %>%
    group_by(from = "SNCA", to = GSE) %>%
    summarise(
      value  = mean(abs(logFC_SNCA)) * 50 + 1,
      status = case_when(
        any(adj.P.Val_SNCA < 0.05 & logFC_SNCA >= 0.6)  ~ "UP_sig",
        any(adj.P.Val_SNCA < 0.05 & logFC_SNCA <= -0.6) ~ "DOWN_sig",
        mean(logFC_SNCA) >= 0                             ~ "UP_ns",
        TRUE                                              ~ "DOWN_ns"
      ),
      .groups = "drop"
    ) %>%
    mutate(color = status_color(status)) %>%
    select(from, to, value, color)
  
  link_resv_sirt1 <- df %>%
    filter(gene == "SIRT1") %>%
    filter(gse_label %in% dataset_order) %>%
    group_by(from = gse_label, to = "SIRT1") %>%
    summarise(value = sum(abs(`interaction score`)) * 100, .groups = "drop") %>%
    mutate(color = "#F39C12DD")
  
  all_links <- bind_rows(
    link_drug_to_gene,
    link_dataset_to_gene,
    link_sirt1,
    link_snca,
    link_resv_sirt1
  )
  

  
  # ==============================================================================
  # Step 8: determine sector order and build the adjacency matrix for the chord diagram.
  # 8. SECTOR ORDER AND MATRIX
  # ==============================================================================
  
  sector_order <- unique(c("RESVERATROL", all_genes_clean, "SIRT1", "SNCA", dataset_order))
  
  mat <- matrix(0,
                nrow = length(sector_order),
                ncol = length(sector_order),
                dimnames = list(sector_order, sector_order)
  )
  
  for (i in seq_len(nrow(all_links))) {
    f <- all_links$from[i]
    t <- all_links$to[i]
    v <- all_links$value[i]
    if (f %in% sector_order && t %in% sector_order) {
      mat[f, t] <- mat[f, t] + v
    }
  }
  
  link_colors_mat <- matrix("#AAAAAA22",
                            nrow = nrow(mat), ncol = ncol(mat),
                            dimnames = dimnames(mat)
  )
  
  for (i in seq_len(nrow(all_links))) {
    f <- all_links$from[i]
    t <- all_links$to[i]
    if (f %in% rownames(link_colors_mat) && t %in% colnames(link_colors_mat)) {
      link_colors_mat[f, t] <- all_links$color[i]
    }
  }
  
  cat("Matris toplam:", sum(mat), "\n")
  cat("Sifir olmayan hucre:", sum(mat > 0), "\n")
  
  # ==============================================================================
  # Step 9: assign sector colors and gap spacing between sectors before rendering.
  # 9. SECTOR COLORS AND GAP
  # ==============================================================================
  
  sector_colors <- c(
    "RESVERATROL" = "#E67E22",
    setNames(rep("#E8A0A0", length(all_genes_clean)), all_genes_clean),
    "SIRT1"       = "#1A5276",
    "SNCA"        = "#154360",
    sapply(dataset_order, function(d) group_colors[dataset_group_map[d]])
  )
  
  gap_vec <- setNames(rep(1, length(sector_order)), sector_order)
  gap_vec["RESVERATROL"]            <- 12
  gap_vec[tail(all_genes_clean, 1)] <- 10
  gap_vec["SNCA"]                   <- 10
  
  grp_seq <- dataset_group_map[dataset_order]
  transitions <- dataset_order[which(diff(
    as.integer(factor(grp_seq, levels = c("Cell-based", "iPSC", "Postmortem")))
  ) != 0)]
  gap_vec[transitions] <- 5
  
  while (!is.null(dev.list())) dev.off()
  circos.clear()
  
  ## ---------------------------------------------------------------------------
  ## SINGLE TIFF --- legend removed, using a single panel.
  ## Canvas size increased (9000x7000 -> 10500x9000) and since the whole image
  ## is a single panel, it can spread across a wider drawing area.
  ## ---------------------------------------------------------------------------
  tiff("Figure3_chord.tiff",
       width       = 10500,
       height      = 9000,
       res         = 600,
       compression = "lzw")
  
  par(mar = c(6, 6, 8, 6), bg = "white", ps = 12)
  
  ## ---------------------------------------------------------------------------
  ## GAP SETTING - spacing between ribbons/sectors SLIGHTLY increased.
  ## A fixed small amount is added to the original gap_vec for each sector.
  ## This is a controlled, non-aggressive increase for the request to open it up "a bit".
  ## ---------------------------------------------------------------------------
  gap_vec_wide <- gap_vec + 0.8
  
  circos.par(
    start.degree            = 90,
    gap.after               = gap_vec_wide,
    track.margin            = c(0.005, 0.005),
    cell.padding            = c(0.02, 0, 0.02, 0),
    points.overflow.warning = FALSE,
    clock.wise              = TRUE
  )
  
  chordDiagram(
    mat,
    order           = sector_order,
    grid.col        = sector_colors,
    col             = link_colors_mat,
    transparency    = 0,
    reduce          = 0,
    link.lwd        = 0.3,
    link.border     = NA,
    annotationTrack = "grid",
    preAllocateTracks = list(
      list(track.height = uh(8, "mm")),
      list(track.height = uh(4, "mm"))
    ),
    directional     = 1,
    direction.type  = c("diffHeight", "arrows"),
    link.arr.type   = "big.arrow",
    link.arr.length = 0.02,
    link.sort       = TRUE,
    link.decreasing = TRUE
  )
  
  ## Track 2 - model group background colors
  circos.trackPlotRegion(
    track.index = 2,
    panel.fun = function(x, y) {
      s <- get.cell.meta.data("sector.index")
      if (!(s %in% dataset_order)) return(invisible())
      grp  <- dataset_group_map[s]
      col  <- group_colors[grp]
      xlim <- get.cell.meta.data("xlim")
      ylim <- get.cell.meta.data("ylim")
      circos.rect(xlim[1], ylim[1], xlim[2], ylim[2],
                  col = col, border = NA)
    },
    bg.border = NA
  )
  
  ## Track 2 - RESVERATROL / SIRT1 / SNCA labels (full 12pt = cex 1.0)
  circos.trackPlotRegion(
    track.index = 2,
    panel.fun = function(x, y) {
      s    <- get.cell.meta.data("sector.index")
      xlim <- get.cell.meta.data("xlim")
      ylim <- get.cell.meta.data("ylim")
      
      if (s %in% c("RESVERATROL", "SIRT1", "SNCA")) {
        circos.text(
          mean(xlim), mean(ylim),
          s,
          facing     = "bending.inside",
          niceFacing = TRUE,
          cex        = 1.0,
          font       = 2,
          col        = "black"
        )
      }
    },
    bg.border = NA
  )
  
  ## Track 1 - combined gene/dataset labels in a SINGLE block
  circos.trackPlotRegion(
    track.index = 1,
    panel.fun = function(x, y) {
      s    <- get.cell.meta.data("sector.index")
      xlim <- get.cell.meta.data("xlim")
      ylim <- get.cell.meta.data("ylim")
      
      if (s %in% c("RESVERATROL", "SIRT1", "SNCA")) {
        circos.text(
          mean(xlim), mean(ylim),
          s,
          facing     = "bending.inside",
          niceFacing = TRUE,
          cex        = 1.0,
          font       = 2,
          col        = "black"
        )
      } else if (s %in% all_genes_clean) {
        circos.text(
          mean(xlim), ylim[2] + 0.0001,
          s,
          facing     = "clockwise",
          niceFacing = TRUE,
          cex        = 0.83,
          font       = 1,
          col        = "#7B241C"
        )
      } else {
        label_split <- sub(" ", "\n", s)
        circos.text(
          mean(xlim), ylim[2] + 0.0001,
          label_split,
          facing     = "clockwise",
          niceFacing = TRUE,
          cex        = 0.75,
          font       = 1,
          col        = "gray10"
        )
      }
    },
    bg.border = NA
  )
  
  circos.clear()
  dev.off()
  
  cat("Saved: Figure3_chord.tiff (no legend, wider gaps, larger canvas)\n")
}