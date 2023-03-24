# ==============================================================================
# 08_DIABLO_misc.R
# robustness analyses and misc figures (e.g., sex-specificity figures/checks)
# (optional checks that don't appear in manuscript main text/analyses)
# ==============================================================================



# combine Fig2A and 2B  --------------------------------------------------------

plot_grid(Fig2A + theme(legend.position = "none"),
          Fig2B, rel_widths = c(4.5, 5.5),
          labels = c("A", "B"))
ggsave("Fig2.png", width = 10, height = 4)




# heatmaps mixing DIABLO models x data -----------------------------------------
# show miRNA-mRNA cor patterns different by sex (Suppl Fig 1)
# e.g., features selected in female model have weaker cors in male samples

# female data, subset to male model/features
crosscor_mat_female <-
  cbind(feature_microarray_mainfx %>%
          .[filter_samples_female, diabfull_mRNAs_male$gene],
        feature_miRNA_mainfx %>%
          .[filter_samples_female, diabfull_miRNAs_male$gene, drop = F]) %>%
  cor(method = "spearman")

png(filename = "FigS1_FemDat_MaleModel.png", height = 800, width = 800, res = 300)
corrplot(crosscor_mat_female, tl.col = "black", tl.cex = 0.25, cl.cex	= 0.2,
         type = "upper", method = "square", addgrid.col = "white")
dev.off()


# male data, subset to female model/features
crosscor_mat_male <-
  cbind(feature_microarray_mainfx %>%
          .[filter_samples_male, diabfull_mRNAs_female$gene],
        feature_miRNA_mainfx %>%
          .[filter_samples_male, diabfull_miRNAs_female$miRNA, drop = F]) %>%
  cor(method = "spearman")

png(filename = "FigS1_MaleDat_FemaleModel.png", height = 800, width = 800, res = 300)
corrplot(crosscor_mat_male, tl.col = "black", tl.cex = 0.25, cl.cex	= 0.2,
         type = "upper", method = "square", addgrid.col = "white")
dev.off()




# DIABLO alt model (miR-1323) --------------------------------------------------
# force with 1323 (commonly came up during sensitivity analyses on ncomp)

diablo_maleonly_alt <-
  block.splsda(  list(mRNA = feature_microarray_mainfx %>%
                        .[filter_samples_male, filter_genes_tidy & filter_genes_multigene],
                      miRNA = feature_miRNA_mainfx[filter_samples_male,
                           DEresults_interact_Males %>% filter(row == "hsa-miR-1323") %>% .$row, drop = F
                      ]) %>%
                   purrr::map(~ apply(.x, 2, function(x) { x - mean(x) / sd(x) })),
                 Y_male,
                 keepX = tune_malefull,
                 ncomp = 1,
                 design = "full")




# DIABLO null models -----------------------------------------------------------

# set.seed(1234)
# tune_femalenull <-
#   tune.block.splsda(X_female, Y_female,
#                     design = "null",
#                     ncomp = 1,
#                     test.keepX = list(mRNA = seq(5, 150, 10),
#                                       miRNA = seq(1, 10, 1)))
tune_femalenull <- list(mRNA = 85, miRNA = 1)
diablo_femaleonly_null <-
  block.splsda(X_female, Y_female,
               keepX = tune_femalenull,
               ncomp = 1,
               design = "null")


# set.seed(1234)
# tune_malenull <-
#   tune.block.splsda(X_male, Y_male,
#                     design = "null",
#                     ncomp = 1,
#                     test.keepX = list(mRNA = seq(5, 150, 10),
#                                       miRNA = seq(1, 10, 1)))
tune_malenull <- list(mRNA = 25, miRNA = 1)
diablo_maleonly_null <-
  block.splsda(X_male, Y_male,
               keepX = tune_malenull,
               ncomp = 1,
               design = "null")

# auc
perform_maleonly_null <-
  perf(object = diablo_maleonly_null, folds = 5, nrepeat = 20,
       progressBar = T, auc = T)
perform_femaleonly_null <-
  perf(object = diablo_femaleonly_null, folds = 5, nrepeat = 20,
       progressBar = T, auc = T)




# all-samples model ------------------------------------------------------------

X_all <-
  list(mRNA = feature_microarray_mainfx %>%
         .[T, filter_genes_tidy & filter_genes_multigene],
       miRNA = feature_miRNA_mainfx[T,
                                    DEresults_IUS %>% filter(qval < 0.2) %>% .$row
       ]) %>%
  purrr::map(~ apply(.x, 2, function(x) { x - mean(x) / sd(x) }))
Y_all <- metadata_miRNA$SmokeBinary %>% as.factor %>% .[T]

# set.seed(1234)
# tune_all_full <-
#   tune.block.splsda(X_all, Y_all,
#                     design = "full",
#                     ncomp = 1,
#                     test.keepX = list(mRNA = seq(5, 150, 10),
#                                       miRNA = seq(1, 10, 1)))
tune_all_full <- list(mRNA = 85, miRNA = 4)
diablo_all_full <-
  block.splsda(X_all, Y_all,
               keepX = tune_all_full,
               ncomp = 1,
               design = "full")

tune_all_null <- list(mRNA = 145, miRNA = 1)
diablo_all_null <-
  block.splsda(X_all, Y_all,
               keepX = tune_all_null,
               ncomp = 1,
               design = "null")




# compare gene lists -----------------------------------------------------------

# euler diagram
venn(combinations =
       list(FemaleFull = selectVar(diablo_femaleonly_full, block = 'mRNA', comp = 1)$mRNA$name,
            MaleFull = selectVar(diablo_maleonly_full, block = 'mRNA', comp = 1)$mRNA$name,
            MaleNull = selectVar(diablo_maleonly_full, block = 'mRNA', comp = 1)$mRNA$name,
            FemaleNull = selectVar(diablo_maleonly_null, block = 'mRNA', comp = 1)$mRNA$name,
            # MaleAlt = selectVar(diablo_maleonly_alt, block = 'mRNA', comp = 1)$mRNA$name,
            # AllNull = selectVar(diablo_all_null, block = 'mRNA', comp = 1)$mRNA$name,
            AllFull = selectVar(diablo_all_full, block = 'mRNA', comp = 1)$mRNA$name)) %>%
  plot

# limited correlation of loadings
loadings_compare_male <-
  as_tibble(selectVar(diablo_maleonly_null, block = 'mRNA', comp = 1)$mRNA$value,
            rownames = "gene") %>%
  left_join(.,
            as_tibble(selectVar(diablo_maleonly_full, block = 'mRNA', comp = 1)$mRNA$value,
                      rownames = "gene"), by = "gene", suffix = c("null", "full")) %>%
  mutate_all(function(x) { x[is.na(x)] <- 0; x })

ggplot(loadings_compare_male,
       aes(value.varfull, value.varnull)) +
  geom_point() +
  theme_few() +
  xlab("Null Loading") +
  ylab("Full Loading")




# heatmap of loadings ----------------------------------------------------------
# variable importance / suppl fig 4

loadings_compare_df <-
  as_tibble(selectVar(diablo_femaleonly_null, block = 'mRNA', comp = 1)$mRNA$value,
            rownames = "gene") %>%
  full_join(.,
            as_tibble(selectVar(diablo_femaleonly_full, block = 'mRNA', comp = 1)$mRNA$value,
                      rownames = "gene"), by = "gene") %>%
  full_join(.,
            as_tibble(selectVar(diablo_maleonly_full, block = 'mRNA', comp = 1)$mRNA$value,
                      rownames = "gene"), by = "gene") %>%
  full_join(.,
            as_tibble(selectVar(diablo_maleonly_null, block = 'mRNA', comp = 1)$mRNA$value,
                      rownames = "gene"), by = "gene")  %>%
  full_join(.,
            as_tibble(selectVar(diablo_maleonly_alt, block = 'mRNA', comp = 1)$mRNA$value,
                      rownames = "gene"), by = "gene") %>%
  full_join(.,
            as_tibble(selectVar(diablo_all_full, block = 'mRNA', comp = 1)$mRNA$value,
                      rownames = "gene"), by = "gene") %>%
  full_join(.,
            as_tibble(selectVar(diablo_all_null, block = 'mRNA', comp = 1)$mRNA$value,
                      rownames = "gene"), by = "gene") %>%
  set_names(c("gene", "Female (Null)", "Female (Full)",
              "Male (Full)", "Male (Null)", "Male (Alt)",
              "All-Samples (Full)", "All-Samples (Null)"))

png("FigS4_loading_heatmap.png", width = 800, height = 1000, res = 300)
Heatmap(matrix = loadings_compare_df %>%
          dplyr::select(-1) %>% as.matrix() %>% replace_na(replace = 0),
        name = "Variable\nImportance",
        col = colorRamp2(c(0, 0.05, 0.15), c("white", "#4f3475", "#1c0340")),
        column_names_gp = gpar(fontsize = 8),
        heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface = "plain"),
                                    labels_gp = gpar(fontsize = 6)))
dev.off()




# plot loadings (not ultimately included in the paper) -------------------------

loadings_df_clean <- 
  loadings_df %>% pivot_wider %>% .[rowSums(.[ , 2:ncol(.)]) != 0 , ] %>%
  pivot_longer(cols = 2:ncol(.))

ggplot(loadings_df_clean %>% filter(grepl('^female \\(', name)) %>% filter(value != 0) %>% 
         arrange(abs(value)) %>% mutate(mRNA = as.factor(mRNA) %>% fct_inorder()),
       aes(x = mRNA, y = value, fill = value > 0)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 0) +
  scale_fill_manual(values = palette_color_IUS) +
  coord_flip() +
  theme_few(base_size = 8) +
  facet_grid( ~ name) +
  theme(legend.position = "none") +
  ylab("Loading (Positive = IUS-Associated)") +
  xlab("DIABLO-Selected mRNAs")

