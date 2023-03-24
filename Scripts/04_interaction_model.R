# ==============================================================================
# 04_interaction_model.R
# add IUS x sex interact
# ==============================================================================




# run RUV for interact model ------------------------------------------------

preruv_formula_interact <-
  as.formula("~ AgeLin + AgeSq + Sex + SmokeBinary + SmokeBinary*Sex + run")

preruv_resid_interact <- 
  apply(feature_miRNA_vst, 1, 
      function(y) {
        lm(y ~ AgeLin + AgeSq + Sex + run + SmokeBinary + SmokeBinary*Sex,
           metadata_miRNA) %>% resid
      })

ruv_plusinteract <-
  RUVr(x = feature_miRNA,
       k = 4,
       isLog = T,
       residuals = preruv_resid_interact %>% t)

metadata_miRNA_interact <-
  metadata_miRNA %>%
  bind_cols(., ruv_plusinteract$W %>%
              as_tibble %>% set_names(paste0("RUVinteract", 1:4)))




# run model including RUV ------------------------------------------------------

base_formula_interact <-
  update(preruv_formula_interact,  ~ . +
           RUVinteract1 + RUVinteract2 + RUVinteract3 + RUVinteract4)

deseq_obj_interact <-
  DESeqDataSetFromMatrix(
    countData = feature_miRNA,
    design = base_formula_interact,
    colData = metadata_miRNA_interact
  ) %>%
  estimateSizeFactors() %>%
  estimateDispersions(fitType = "local") %>%
  nbinomWaldTest(maxit = deseqiterations)




# significant interacts -----------------------------------------------------
# turned off independent filtering to have the same miRNAs tested
# as main effect model in script 03
# i.e., H0: IUS x sex effect = 0

DEresults_interact_Exists <-
  deseq_obj_interact %>% 
  results(., name = "Sex.SmokeBinary", tidy = T, independentFiltering = F) %>%
  arrange(pvalue) %>%
  filter(!is.na(padj)) %>%
  filter(row %in% DEresults_IUS$row)

DEresults_interact_Exists$qval <-
  qvalue(DEresults_interact_Exists$pvalue)$qvalue

DEresults_interact_Exists %>%
  filter(qval < 0.10)




# significantly non-zero IUS in each sex ---------------------------------------
# e.g., H0: within females, mean IUS effect = 0
# use contrasts to extract relevant terms

DEresults_interact_Males <-
  deseq_obj_interact %>% 
  results(., contrast = as.numeric(resultsNames(deseq_obj_interact) %in%
                                     c("SmokeBinary", "Sex.SmokeBinary")),
          tidy = T, independentFiltering = F) %>%
  arrange(pvalue) %>%
  filter(!is.na(padj)) %>%
  filter(row %in% DEresults_IUS$row)

DEresults_interact_Males$qval <-
  qvalue(DEresults_interact_Males$pvalue)$qvalue

DEresults_interact_Females <-
  deseq_obj_interact %>% 
  results(., name = "SmokeBinary", tidy = T, independentFiltering = F) %>%
  arrange(pvalue) %>%
  filter(!is.na(padj)) %>%
  filter(row %in% DEresults_IUS$row)

DEresults_interact_Females$qval <-
  qvalue(DEresults_interact_Females$pvalue)$qvalue




# checking correlation of IUS effects ------------------------------------------

DEresults_interact_Collected <-
  left_join(DEresults_interact_Males,
          DEresults_interact_Females,
          suffix = c(".m", ".f"),
          by = "row") %>%
  mutate(Significance =
           case_when(qval.m < 0.1 & qval.f > 0.1 ~ "Male",
                     qval.m > 0.1 & qval.f < 0.1 ~ "Female",
                     qval.m < 0.1 & qval.f < 0.1 ~ "Both",
                     qval.m > 0.1 & qval.f > 0.1 ~ "Neither"))

# spearman cor = -0.12
cor(DEresults_interact_Collected$log2FoldChange.m[
  DEresults_interact_Collected$Significance != "Neither"],
    DEresults_interact_Collected$log2FoldChange.f[
      DEresults_interact_Collected$Significance != "Neither"],
    method = "spearman")

# spearman cor = -0.011
cor(DEresults_interact_Collected$pvalue.f,
    DEresults_interact_Collected$pvalue.m,
    method = "spearman")




# compile IUS effect estimates for plotting ------------------------------------
# shows limited correlation of male IUS and female IUS estimates

DEresults_interact_Collected <-
  DEresults_interact_Collected %>%
  left_join(annotation_miRNAs, by = c("row"= "miRNA"))

interact_largest_fx <-
  DEresults_interact_Collected %>%
  arrange(log2FoldChange.m) %>%
  bind_cols(index = 1:nrow(.)) %>%
  mutate(row = if_else(index %in% c(1:6, 510:515) & (qval.m < 0.10 | qval.f < 0.10),
                       gsub("hsa-miR-", "", row), ""))

Fig1B <-
    ggplot(data = DEresults_interact_Collected, # %>% filter(Significance != "Neither"),
         aes(log2FoldChange.m, log2FoldChange.f,
             color = Significance, shape = Significance)) +
  geom_abline(slope = 1, intercept = 0, lty = 3) +
  geom_point(alpha = 0.5, size = 2) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_text_repel(data = interact_largest_fx,
                  aes(log2FoldChange.m, log2FoldChange.f, label = row),
                  size = 3.5, min.segment.length = 0.15, show.legend = F, box.padding = 0.15) +
  theme_few(base_size = 14) +
  scale_color_manual(values = palette_color_sex_signif) +
  scale_shape_manual(values = palette_color_shape_signif) +
  scale_x_continuous(expression("IUS"~log[2]~"FC (Male)"), limits = c(-1.1, 1.1)) +
  scale_y_continuous(name = expression("IUS"~log[2]~"FC (Female)"), limits = c(-1.1, 1.1))

ggsave(Fig1B, filename = "Fig1B.png", width = 5, height = 4, dpi = 300)
  
# double check that miRNAs on sex chromosomes (since high density on chrX)
# aren't driving the majority of sex-specific effects -- seems OK
ggplot(data = DEresults_interact_Collected %>% filter(grepl("X", Chr)),
       aes(log2FoldChange.m, log2FoldChange.f,
           color = Significance, shape = Significance)) +
  geom_abline(slope = 1, intercept = 0, lty = 3) +
  geom_point(alpha = 0.5, size = 2) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  theme_few() +
  scale_color_manual(values = palette_color_sex_signif) +
  scale_shape_manual(values = palette_color_shape_signif) +
  scale_x_continuous(expression("IUS"~log[2]~"FC (Male)"), limits = c(-1.1, 1.1)) +
  scale_y_continuous(name = expression("IUS"~log[2]~"FC (Female)"), limits = c(-1.1, 1.1))




# checking post-hoc for any heteroscedasticity in sex --------------------------
# no clear systematic higher variance in female miRNA expression

ggplot(metadata_miRNA,
       aes(x = as.factor(Sex),
           y = feature_miRNA_vst[DEresults_interact_Exists$row[2], ])) +
  ylab(DEresults_interact_Exists$row[2]) +
  stat_summary() +
  theme_few()

left_join(DEresults_interact_Males,
          DEresults_interact_Females,
          suffix = c(".m", ".f"),
          by = "row") %>% mutate(ratio = lfcSE.f/lfcSE.m) %>%
  .$ratio %>% range


bind_rows(DEresults_interact_Males %>% bind_cols(Sex = "Male"),
          DEresults_interact_Females %>% bind_cols(Sex = "Female")) %>%
  ggplot(data = ., aes(x = abs(log2FoldChange), y = Sex)) +
  stat_density_ridges(aes(fill = Sex), alpha = 0.25, quantile_lines = T, 
                      quantiles = c(0.5), alpha = 0.7, size = 0.2) +
  theme_ridges() + theme(legend.position = "none") +
  scale_fill_manual(values = palette_color_sex_signif)

bind_rows(DEresults_interact_Males %>% bind_cols(Sex = "Male"),
          DEresults_interact_Females %>% bind_cols(Sex = "Female")) %>%
  ggplot(data = ., aes(x = abs(log2FoldChange), color = Sex)) +
  geom_vline(xintercept = 0, lty = 3) +
  stat_ecdf() +
  theme_few() +
  scale_color_manual(values = palette_color_sex_signif[c("Male", "Female")])

bind_rows(DEresults_interact_Males %>% bind_cols(Sex = "Male"),
          DEresults_interact_Females %>% bind_cols(Sex = "Female")) %>%
  ggplot(data = ., aes(x = lfcSE, color = Sex)) +
  stat_ecdf() +
  theme_few() + 
  scale_color_manual(values = palette_color_sex_signif[c("Male", "Female")])




# venn diagram showing number of IUS-miRNAs by sex -----------------------------

png(filename = "Fig1A.png", width = 1200, height = 1200, res = 300)
venn(
  list(`All-Sample\nModel` = filter(DEresults_IUS, qval < 0.1)$row,
       `Female Samples\n(Interact Model)` = filter(DEresults_interact_Females, qval < 0.1)$row,
       `Male Samples\n(Interact Model)` = filter(DEresults_interact_Males, qval < 0.1)$row)
) %>% plot(quantities = T, main = "# IUS-miRNAs", cex = 0.1)
dev.off()

# miRNAs that pop up only in joint modeling
setdiff(filter(DEresults_IUS, qval < 0.1)$row,
        filter(DEresults_interact_Males, qval < 0.1)$row) %>%
  setdiff(filter(DEresults_interact_Females, qval < 0.1)$row)

Fig1 <-
  plot_grid(
  as.ggplot(image_read("Fig1A.png")),
  Fig1B, labels = c("A", "B"), rel_widths = c(4, 6))
ggsave("Fig1.png", plot = Fig1,
       width = 10, height = 5)



# example code signif residualized miRNA values, -------------------------------
# allowing for sex interaction (done for ATS abstract, Rosenberg, Liu, et al 2022)

feature_miRNA_interactions <- 
  apply(feature_miRNA_vst, 1, 
        function(y) {
          lm(y ~ AgeLin + AgeSq + Sex + run +
               RUVinteract1 + RUVinteract2 + RUVinteract3 + RUVinteract4,
             metadata_miRNA_interact) %>% resid
        })

index_to_plot <- 1

ggplot(data = metadata_miRNA_interact %>%
         mutate(IUS = if_else(SmokeBinary == 1, "yes", "no")),
       aes(x = IUS,
           y = feature_miRNA_interactions[ , DEresults_interact_Exists$row[index_to_plot]],
           color = IUS, shape = IUS)) +
  geom_hline(yintercept = 0, lty = 3) +
  geom_quasirandom(alpha = 0.5) +
  stat_summary(geom = "errorbar", color = "black", size = 1.5, fun.min = median) +
  facet_grid(~ if_else(Sex == T, "Male", "Female")) +
  ggthemes::theme_few(base_size = 12) +
  ylab(DEresults_interact_Exists$row[index_to_plot]) +
  xlab("IUS") +
  scale_color_manual(values = palette_color_IUS) +
  scale_shape_manual(values = palette_shape_IUS) +
  theme(legend.position = "none")


