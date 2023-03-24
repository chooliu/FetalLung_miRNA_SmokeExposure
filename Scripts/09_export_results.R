# ==============================================================================
# 09_export_results.R
# export more human-readable Excel format for collaborator/journal sharing
# also includes calculation of pairwise miRNA-mRNA correlations for tables
# ==============================================================================




# display top mRNA-miRNA correlations given a miRNA input ----------------------

calculate_top_miRNA_correlations <-
  function(mirna_name, npairs = 8, filter_subjects = T) {
    
    calculate_correlations <-
      cor(feature_miRNA_mainfx[filter_subjects , mirna_name],
          feature_microarray_mainfx[filter_subjects, ], method = "spearman") %>%
      t %>% as_tibble(rownames = "symbol")
    
    calculate_correlations %>%
      arrange(-abs(V1)) %>%
      filter(!grepl("/", symbol)) %>% # optional remove multisymbol Affy probes
      (function(x) { x[1:min(npairs, nrow(x)), ] } ) %>%
      # mutate(output = paste0(symbol, " (", formatC(V1, digits = 2), ")")) %>%
      mutate(output = paste0(symbol, if_else(V1 < 0, "-", "+"))) %>%
      .$output %>%
      paste(collapse = "; ")
  }




# main IUS effect --------------------------------------------------------------
# Table 1

FinalResults_IUS <-
  FinalResults_IUS %>%
  rowwise() %>%
  mutate(`Top Microarray Correlations` = calculate_top_miRNA_correlations(miRNA))



# IUS x sex interaction exists -------------------------------------------------
# supplementary tables

FinalResults_InteractExists <-
  DEresults_interact_Exists %>% 
  filter(qval < 0.1) %>%
  left_join(annotation_miRNAs %>% dplyr::select(miRNA, Chr), by = c("row" = "miRNA")) %>%
  arrange(-log2FoldChange) %>%
  transmute(`Most Positive IUS-Effect In` =
              if_else(log2FoldChange > 0, "Male", "Female"),
            miRNA = row,
            Chrom = Chr,
            `Interaction Effect (log2FC)` = log2FoldChange,
            `SE (log2FC)` = lfcSE,
            `p-value` = pvalue,
            `q-value` = qval)  %>%
  mutate_if(is.numeric, format, digits = 2)

# IUS effects by sex  ----------------------------------------------------------
# e.g., H0: IUS effect in females = 0

FinalResults_Female <-
  DEresults_interact_Females %>% 
  filter(qval < 0.1) %>%
  left_join(annotation_miRNAs %>% dplyr::select(miRNA, Chr), by = c("row" = "miRNA")) %>%
  arrange(-log2FoldChange) %>%
  transmute(`Higher Levels In` = if_else(log2FoldChange > 0, "IUS", "non-exposed"),
            miRNA = row,
            Chrom = Chr,
            `IUS Effect (log2FC)` = log2FoldChange,
            `SE (log2FC)` = lfcSE,
            `p-value` = pvalue,
            `q-value` = qval) %>%
  rowwise() %>%
  mutate(`Top Microarray Correlations` =
           calculate_top_miRNA_correlations(miRNA, filter_subjects = metadata_miRNA$Sex == 0)) %>%
  mutate_if(is.numeric, format, digits = 2)

FinalResults_Male <-
  DEresults_interact_Males %>% 
  filter(qval < 0.1) %>%
  left_join(annotation_miRNAs %>% dplyr::select(miRNA, Chr), by = c("row" = "miRNA")) %>%
  arrange(-log2FoldChange) %>%
  transmute(`Higher Levels In` = if_else(log2FoldChange > 0, "IUS", "non-exposed"),
            miRNA = row,
            Chrom = Chr,
            `IUS Effect (log2FC)` = log2FoldChange,
            `SE (log2FC)` = lfcSE,
            `p-value` = pvalue,
            `q-value` = qval) %>%
  rowwise() %>%
  mutate(`Top Microarray Correlations` = 
           calculate_top_miRNA_correlations(miRNA, filter_subjects = metadata_miRNA$Sex == 1)) %>%
  mutate_if(is.numeric, format, digits = 2)




# final export -----------------------------------------------------------------
# each results set in its own tab

write_xlsx(list(
  IUS_AllSamples = FinalResults_IUS,
  Interaction = FinalResults_InteractExists,
  MaleIUS = FinalResults_Male,
  FemaleIUS = FinalResults_Female),
  path = "20230130_IUS_Results.xlsx", format_headers = F)
system('open "20230130_IUS_Results.xlsx"')
