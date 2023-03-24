# ==============================================================================
# 05_run_miEAA_GSEA.R
# export results for miEAA 2.0 webportal (preranked GSEA)
# ==============================================================================




dir.create("miEAA_Input")

# main effects

DEresults_IUS$row %>%
  write_lines(file = "miEAA_Input/all_GSEA.txt")

# interactions

DEresults_interact_Exists$row %>%
  write_lines(file = "miEAA_Input/interact_GSEA.txt")

DEresults_interact_Females$row %>%
  write_lines(file = "miEAA_Input/female_GSEA.txt")

DEresults_interact_Males$row %>%
  write_lines(file = "miEAA_Input/male_GSEA.txt")