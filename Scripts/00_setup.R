# ==============================================================================
# 00_setup.R
# libraries, plotting
# ==============================================================================




# load libraries ---------------------------------------------------------------

# data i/o
library(tidyverse)
library(data.table)
library(writexl)

# analysis
library(RUVSeq)
library(DESeq2)
library(qvalue)
library(biomaRt) # *
library(vegan) # *
library(mixOmics) # *

# plotting
library(eulerr)
library(ggbeeswarm)
library(ggridges)
library(ggrepel)
library(corrplot)
library(ggthemes)
library(tidygraph)
library(ggraph)
library(ComplexHeatmap)
library(circlize)
library(magick)
library(ggplotify)
library(cowplot)


# * starred pkgs - caution: the following may result in masked functions
# namely, specifying dplyr::select() & purrr::map()
# might be needed for some tidyverse operations

# save version control history
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")




# color palettes ---------------------------------------------------------------

palette_color_IUS <-
  c(IUS = "#d7191c", yes = "#d7191c", `TRUE` = "#d7191c",
    no = "#000000", `FALSE` = "#000000") # 2c7bb6

palette_shape_IUS <-
  c(IUS = 19, yes = 19, `TRUE` = 19,
    no = 15, `FALSE` = 15)

palette_color_sex_signif <-
  c(Female = "#e41a1c", Male = "#377eb8", Neither = "#606060")

palette_color_shape_signif <-
  c(Female = 19, Male = 19, Neither = 1)
