# ==============================================================================
# 03_primary_model.R
# perform RUVr, DESeq2 main effects model (IUS vs unexposed, no interactions)
# ==============================================================================




# RUVr calculations ------------------------------------------------------------

# initial deseq model to obtain vst residuals
preruv_formula_mainfx <-
  as.formula("~ AgeLin + AgeSq + Sex + SmokeBinary + run")
metadata_final_mainfx <-
  metadata_miRNA

vst_counts_prelim_mainfx <-
  DESeqDataSetFromMatrix(
    countData = feature_miRNA,
    design = preruv_formula_mainfx,
    colData = metadata_final_mainfx
  ) %>%
  estimateSizeFactors() %>%
  estimateDispersions(fitType = "local") %>%
  varianceStabilizingTransformation(fitType = "local")

vst_residuals_mainfx <-
  vst_counts_prelim_mainfx %>%
  assay %>%
  apply(., 1,
        function(y) {
          tmp <-
            metadata_final_mainfx %>%
            bind_cols(y = y)
          lm(update(preruv_formula_mainfx, y ~ .),
             tmp) %>% 
            resid() }
  )

# check RUVr factor explanatory power on vst -- results in k = 4
check_ruv_kvals <-
  function(k) {
    ruv_tmp <-
      RUVr(x = assay(vst_counts_prelim_mainfx),
           k = k,
           residuals = vst_residuals_mainfx %>% t,
           isLog = T)$W
    adonis(
      vst_residuals_mainfx ~ ruv_tmp,
      method = "manhattan") %>%
      .$aov.tab %>%
      .$R2 %>%
      .[1]
  }

# choose k = 4 based on elbow plot
r2_by_k <- sapply(1:20, check_ruv_kvals)
ggplot(data = NULL, aes(x = 1:20, y = r2_by_k)) +
  geom_point() +
  geom_line() +
  xlab("# of RUVr Components (k)") +
  ylab("% of Residual Variance Explained\n(Multivariate PERMANOVA)") +
  theme_classic()

ruv_soln_mainfx <-
  RUVr(x = assay(vst_counts_prelim_mainfx),
       k = 4,
       residuals = vst_residuals_mainfx %>% t,
       isLog = T)

metadata_final_mainfx <-
  metadata_miRNA %>%
  bind_cols(
    ruv_soln_mainfx$W %>%
      as_tibble() %>%
      set_names(paste0("RUV", 1:4)))




# DESeq2 model (#1) ------------------------------------------------------------
# main effects only

base_formula_mainfx <-
  update(preruv_formula_mainfx,  ~ . + RUV1 + RUV2 + RUV3 + RUV4)

deseqiterations <- 10000

deseq_obj_mainfx <-
  DESeqDataSetFromMatrix(
    countData = feature_miRNA,
    design = base_formula_mainfx,
    colData = metadata_final_mainfx
  ) %>%
  estimateSizeFactors() %>%
  estimateDispersions(fitType = "local") %>%
  nbinomWaldTest(maxit = deseqiterations)

DEresults_IUS <-
  deseq_obj_mainfx %>% 
  .[mcols(.)$betaConv, ] %>%
  results(., name = "SmokeBinary", tidy = T, independentFiltering = T) %>%
  arrange(pvalue) %>%
  filter(!is.na(padj))

DEresults_IUS$qval <-
  qvalue(DEresults_IUS$pvalue)$qvalue

DEresults_IUS$FC <- 2^DEresults_IUS$log2FoldChange
DEresults_IUS$HigherIn <- if_else(DEresults_IUS$log2FoldChange > 0, "IUS", "noExposure")

FinalResults_IUS <-
  DEresults_IUS %>% filter((qval < 0.1)) %>%
  left_join(annotation_miRNAs %>% dplyr::select(miRNA, Chr),
            by = c("row" = "miRNA")) %>%
  arrange(-log2FoldChange) %>%
  transmute(`Higher Levels In` = if_else(log2FoldChange > 0, "IUS", "non-exposed"),
            miRNA = row,
            Chrom = Chr,
            # `Mean Count (non-IUS)` = baseMean,
            # `Mean Count (IUS-exposed)` = baseMean * FC,
            `Effect (log2FC)` = log2FoldChange,
            `SE (log2FC)` = lfcSE,
            # `LRT Statistic` = stat,
            `p-value` = pvalue,
            `q-value` = qval) %>%
  mutate_if(is.numeric, format, digits = 2) 




# residualize for plotting -----------------------------------------------------

# miRNA
feature_miRNA_vst <-
  varianceStabilizingTransformation(deseq_obj_mainfx, fitType = "local") %>%
  assay

feature_miRNA_mainfx <- 
  apply(feature_miRNA_vst, 1, 
        function(y) {
          lm(y ~ AgeLin + AgeSq + Sex + run + RUV1 + RUV2 + RUV3 + RUV4,
             metadata_final_mainfx) %>% resid
        })

# microarray
# k = 4 based on previous work / comparability
# (similar multivariate percent explained like above)
# no run variable since diff technical effects
microarray_residuals_mainfx <-
  feature_microarray %>%
  apply(., 1,
        function(y) {
          tmp <-
            metadata_miRNA %>%
            bind_cols(y = y)
          lm(update(preruv_formula_mainfx, y ~ . - run),
             tmp) %>% 
            resid() }
  )

ruv_microarray_mainfx <-
  RUVr(x = feature_microarray,
       k = 4,
       residuals = t(microarray_residuals_mainfx),
       isLog = T)
metadata_microarray <-
  metadata_miRNA %>% dplyr::select(-contains("RUV")) %>%
  bind_cols(ruv_microarray_mainfx$W %>% as_tibble %>% set_names(paste0("RUV", 1:4)))

feature_microarray_mainfx <- 
  apply(feature_microarray, 1, 
        function(y) {
          lm(y ~ AgeLin + AgeSq + Sex + run + RUV1 + RUV2 + RUV3 + RUV4,
             metadata_microarray) %>% resid
        })


