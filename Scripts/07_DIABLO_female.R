# ==============================================================================
# 07_DIABLO_female.R
# miRNA-mRNA integration in fetal lung samples, female-only
# ==============================================================================




# prep miRNA and mRNA for female samples ---------------------------------------

# female-only
filter_samples_female <-
  metadata_miRNA$Sex == 0

# little extra filtering of mRNA probesets here
# one probeset / gene only, no probesets missing gene symbols
filter_genes_tidy <- !duplicated(colnames(feature_microarray_mainfx))
filter_genes_multigene <- map_lgl(colnames(feature_microarray_mainfx),
                                  ~ grepl("/", .x) %>% `!`)

# also filtering of miRNAs for those with some nominal IUS-effect (q < 0.20)
# for interpretability purposes
X_female <-
  list(mRNA = feature_microarray_mainfx %>%
         .[filter_samples_female, filter_genes_tidy & filter_genes_multigene],
       miRNA = feature_miRNA_mainfx[
         filter_samples_female,
         DEresults_interact_Females %>% filter(qval < 0.20) %>% .$row]) %>%
  purrr::map(~ apply(.x, 2, function(x) { x - mean(x) / sd(x) }))
Y_female <- metadata_miRNA$SmokeBinary %>% as.factor %>% .[filter_samples_female]




# grid search: final values pasted below for reproducibility -------------------

# set.seed(1234)
# tune_femalefull <-
#   tune.block.splsda(X_female, Y_female,
#                     design = "full",
#                     ncomp = 1,
#                     test.keepX = list(mRNA = seq(5, 150, 10),
#                                       miRNA = seq(1, 10, 1)))
# check tuning
# tune_femalefull$choice.keepX
# enframe(tune_femalefull$error.rate) %>% arrange(value[, 'comp1'])

tune_femalefull <- list(mRNA = 65, miRNA = 2) # female only




# run DIABLO full --------------------------------------------------------------

diablo_femaleonly_full <-
  block.splsda(X_female, Y_female,
               keepX = tune_femalefull,
               ncomp = 1,
               design = "full")

# print variables
selectVar(diablo_femaleonly_full, block = 'miRNA', comp = 1)$miRNA$name 
selectVar(diablo_femaleonly_full, block = 'mRNA', comp = 1)$mRNA$name %>%
  sort %>% paste(collapse = "; ") 

# check AUC, stable features
perform_femaleonly_full <-
  perf(object = diablo_femaleonly_full, folds = 5, nrepeat = 50,
       progressBar = T, auc = T) 

perform_femaleonly_full$auc
perform_femaleonly_full$features$stable

diabfull_mRNAs_female <- diablo_femaleonly_full$loadings$mRNA %>%
  enframe %>% set_names(c("gene", "loading")) %>% filter(loading != 0)

diabfull_miRNAs_female <- diablo_femaleonly_full$loadings$miRNA %>%
  enframe %>% set_names(c("miRNA", "loading")) %>% filter(loading != 0)




# miRNA-mRNA correlations (heatmap display) ------------------------------------

cor_matrix_fem <-
  cbind(feature_microarray_mainfx %>%
          .[filter_samples_female, diabfull_mRNAs_female$gene],
        feature_miRNA_mainfx %>%
          .[filter_samples_female, diabfull_miRNAs_female$miRNA]) %>%
  cor(method = "spearman")

png(filename = "FigS1_FemDat_FemModel.png", height = 800, width = 800, res = 300)
corrplot(cor_matrix_fem, tl.col = "black", tl.cex = 0.25, cl.cex = 0.2,
         type = "upper", method = "square", addgrid.col = "white")
dev.off()




# miRNA-mRNA correlations (network display) ------------------------------------

# edge thresholds; note in that Supplementary Methods / sensitivity analyses
# use lowest 10th of negative corrs, highest 90th percentiles of positive corrs to prune
# (positive != 1 to avoid cormat diagonals)
# robustness of figures was checked by manually changing these thresholds (Suppl. Fig 3)
cor_matrix_fem[lower.tri(cor_matrix_fem, diag = T)] <- NA
cor_range_min <- -0.25 # quantile(c(cor_matrix_male[cor_matrix_male < 0], cor_matrix_fem[cor_matrix_fem < 0]), 0.10, na.rm = T)
cor_range_max <- 0.40 # quantile(c(cor_matrix_male[cor_matrix_male > 0 & cor_matrix_male != 1], cor_matrix_fem[cor_matrix_fem > 0 & cor_matrix_fem != 1]), 0.90, na.rm = T) 

# ggraph / tidygraph
# * note: two mRNAs manually removed b/c correlated with one another,
# but no other members of the networks after edge filtering
tidy_graph_fem <-
  as_tbl_graph(cor_matrix_fem, directed = F) %>%
  activate(edges) %>%
  filter(weight < cor_range_min | weight > cor_range_max) %>%
  filter(!is.na(weight)) %>%
  activate(nodes)  %>%
  mutate(deg = centrality_degree(),
         feature = if_else(grepl("miR", name), "miRNA", "mRNA")) %>%
  filter(deg != 0 & !name %in% c("PKNOX2", "USP36")) # *

# final figure
set.seed(1234)
Fig2B <-
  ggraph(tidy_graph_fem, layout = 'auto') +
  geom_edge_fan(aes(width = abs(weight), color = weight)) +
  geom_node_point(size = 3.5, aes(shape = feature)) +
  geom_node_text(aes(label = name), repel = T, size = 3) +
  scale_edge_colour_gradient2(name = "r", limits = c(-0.7, 0.7),
                              breaks = c(-0.6, -0.3, 0, 0.3, 0.6)) +
  scale_shape_manual(values = c(miRNA = 19, mRNA = 1)) +
  theme_void() +
  scale_edge_width(name = "|r|", range = c(0.5, 2.5), limits = c(0.2, 0.6))

ggsave("Fig2B.png", plot = Fig2B, width = 6, height = 5)

