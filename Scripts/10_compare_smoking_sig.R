# ==============================================================================
# 10_compare_smoking_sig.R
# comparison to existing studies of IUS or tobacco exposure
# ==============================================================================





# add literature lists ---------------------------------------------------------

# Vives-Usano, M.; Hernandez-Ferrer, C.; Maitre, L.; Ruiz-Arenas, C.; Andrusaityte,
# S.; Borràs, E.; Carracedo, Á.; Casas, M.; Chatzi, L.; Coen, M.; et al.
# In utero and childhood exposure to tobacco smoke and multi-layer molecular
# signatures in children. BMC Med. 2020, 18, 1–19.

miRNAs_VivesUsano2020 <- "none"


# Herberth, G.; Bauer, M.; Gasch, M.; Hinz, D.; Röder, S.; Olek, S.; Kohajda, T.;
# Rolle-Kampczyk, U.; Von Bergen, M.; Sack, U.; et al. Maternal and cord blood
# miR-223 expression associates with prenatal tobacco smoke exposure and low
# regulatory T-cell numbers. J. Allergy Clin. Immunol. 2014, 133.

miRNAs_Herberth2014 <- c("hsa-miR-155", "hsa-miR-223")


# Sekovanić, A.; Dorotić, A.; Pašalić, D.; Orct, T.; Kljaković-Gašpić, Z.;
# Grgec, A.S.; Stasenko, S.; Mioč, T.; Piasek, M.; Jurasović, J. The effects of
# maternal cigarette smoking on cadmium and lead levels, miRNA expression and
# biochemical parameters across the feto-placental unit. Heliyon 2022, 8.

miRNAs_Sekovanic2022 <-
  c("hsa-miR-1537", "hsa-miR-190b", "hsa-miR-16", "hsa-miR-21", "hsa-miR-146a")


# Juracek, J.; Piler, P.; Janku, P.; Radova, L.; Slaby, O. Identification of
# microRNA signatures in umbilical cord blood associated with maternal
# characteristics. PeerJ 2019, 2019, 1–11.

miRNAs_Juracek2019 <- c(
  "hsa-miR-129-5p", "hsa-miR-30b-3p", "hsa-miR-187-3p", "hsa-miR-507",
  "hsa-miR-520b", "hsa-miR-33b-3p", "hsa-miR-138-1-3p", "hsa-miR-760"
)


# Dehmel, S.; Nathan, P.; Bartel, S.; El-Merhie, N.; Scherb, H.; Milger, K.;
# John-Schuster, G.; Yildirim, A.O.; Hylkema, M.; Irmler, M.; et al.
# Intrauterine smoke exposure deregulates lung function, pulmonary
# transcriptomes, and in particular insulin-like growth factor (IGF)-1
# in a sex-specific manner. Sci. Rep. 2018, 8, 1–12.

miRNAs_Dehmel2018 <- c(
  "mmu-miR-5117",
  "mmu-miR-181a-2",
  "mmu-miR-374",
  "mmu-miR-7a-1",
  "mmu-miR-539-5p",
  "mmu-miR-3962",
  "mmu-miR-3057-5p",
  "mmu-miR-374c",
  "mmu-miR-449a",
  "mmu-miR-183",
  "mmu-miR-714",
  "mmu-miR-4671",
  "mmu-miR-542-3p",
  "mmu-miR-376c",
  "mmu-miR-433",
  "mmu-miR-5099",
  "mmu-miR-669e",
  "mmu-miR-412-5p",
  "mmu-miR-130b",
  "mmu-miR-98",
  "mmu-miR-466e-3p",
  "mmu-miR-1983",
  "mmu-miR-676",
  "mmu-miR-5102",
  "mmu-miR-466d-3p",
  "mmu-miR-466p-3p",
  "mmu-miR-543",
  "mmu-miR-466b-3p",
  "mmu-miR-187",
  "mmu-miR-5110",
  "mmu-mir-432",
  "mmu-mir-5117",
  "mmu-miR-495",
  "mmu-miR-1193-5p",
  "mmu-miR-339-3p",
  "mmu-miR-497",
  "mmu-let-7i",
  "mmu-miR-1843b-5p",
  "mmu-miR-34c",
  "mmu-miR-1186",
  "mmu-miR-3473b",
  "mmu-miR-141",
  "mmu-miR-486",
  "mmu-miR-3107",
  "mmu-miR-34c",
  "mmu-miR-124",
  "mmu-miR-24-1",
  "mmu-miR-208a-5p",
  "mmu-miR-201",
  "mmu-miR-206"
) %>%
  gsub("mmu", "hsa", .)

# Willinger, C.M.; Rong, J.; Tanriverdi, K.; Courchesne, P.L.; Huan, T.;
# Wasserman, G.A.; Lin, H.; Dupuis, J.; Joehanes, R.; Jones, M.R.; et al.
# MicroRNA Signature of Cigarette Smoking and Evidence for a Putative Causal
# Role of MicroRNAs in Smoking-Related Inflammation and Target Organ Damage.
# Circ. Cardiovasc. Genet. 2017, 10, 1–12.

miRNAs_Willinger2017 <- c(
  "hsa-miR-1180", "hsa-miR-181a-2-3p", "hsa-miR-423-5p",
  "hsa-miR-25-5p", "hsa-miR-1285-3p", "hsa-miR-342-5p"
)


# Wang, G.; Wang, R.; Strulovici-Barel, Y.; Salit, J.; Staudt, M.R.;
# Ahmed, J.; Tilley, A.E.; Yee-Levin, J.; Hollmann, C.; Harvey, B.G.; et al.
# Persistence of smoking-induced dysregulation of MiRNA expression in the small
# airway epithelium despite smoking cessation. PLoS One 2015. 10, e0120824.

miRNAs_Wang2015 <-
  c(
    "hsa-miR-143",
    "hsa-miR-145",
    "hsa-miR-133b",
    "hsa-miR-214",
    "hsa-miR-634",
    "hsa-miR-126",
    "hsa-miR-139",
    "hsa-miR-199a-3p",
    "hsa-miR-199a-5p",
    "hsa-miR-133a",
    "hsa-miR-195",
    "hsa-miR-181c",
    "hsa-miR-675",
    "hsa-miR-127",
    "hsa-let-7b",
    "hsa-miR-1260",
    "hsa-miR-1226",
    "hsa-miR-636",
    "hsa-miR-193b",
    "hsa-miR-138-1",
    "hsa-miR-181a",
    "hsa-miR-550",
    "hsa-miR-181b",
    "hsa-miR-449b",
    "hsa-miR-224",
    "hsa-miR-1975",
    "hsa-miR-1979",
    "hsa-miR-218",
    "hsa-miR-146a",
    "hsa-miR-203",
    "hsa-miR-3201",
    "hsa-miR-1246"
  )


# Schembri, F.; Sridhar, S.; Perdomo, C.; Gustafson, A.M.; Zhang, X.;
# Ergun, A.; Lu, J.; Liu, G.; Zhang, X.; Bowers, J.; et al.
# MicroRNAs as modulators of smoking-induced gene expression changes in
# human airway epithelium. Proc. Natl. Acad. Sci. U. S. A. 2009, 106, 2319–2324.

miRNAs_Schembri2009 <-
  c(
    "hsa-miR-337",
    "hsa-miR-18a",
    "hsa-miR-189",
    "hsa-miR-365",
    "hsa-miR-181d",
    "hsa-miR-10b",
    "hsa-miR-150",
    "hsa-miR-218",
    "hsa-miR-338",
    "hsa-miR-362",
    "hsa-miR-17-3p",
    "hsa-miR-15a",
    "hsa-miR-652",
    "hsa-miR-106b",
    "hsa-miR-19b",
    "hsa-miR-106a",
    "hsa-miR-128a",
    "hsa-miR-30a-3p",
    "hsa-miR-128b",
    "hsa-miR-130a",
    "hsa-miR-500",
    "hsa-miR-363",
    "hsa-miR-199b",
    "hsa-miR-223",
    "hsa-miR-625",
    "hsa-miR-99a",
    "hsa-miR-125b",
    "hsa-miR-146a"
  )


# Huang, J.; Jiang, W.; Tong, X.; Zhang, L.; Zhang, Y.; Fan, H.; Ding, J.
# Identification of gene and microRNA changes in response to smoking in human
# airway epithelium by bioinformatics analyses. Med. (United States) 2019, 98.

miRNAs_Huang2019 <- c(
  "hsa-miR-627-5p",
  "hsa-miR-500a-3p",
  "hsa-miR-218-5p",
  "hsa-miR-9-5p",
  "hsa-miR-212-5p",
  "hsa-miR-143-3p",
  "hsa-miR-106b-3p"
)


# Izzotti, A.; Calin, G.A.; Arrigo, P.; Steele, V.E.; Croce, C.M.; De Flora, S.
# Downregulation of microRNA expression in the lungs of rats exposed to
# cigarette smoke. FASEB J. 2009, 23, 806–812.

miRNAs_Izzotti2009 <- c(
  "let-7a", "let-7b", "let-7c", "let-7f",
  "miR-10a", "miR-26a", "miR-30a", "miR-30c", "miR-34b", "miR-34c",
  "miR-99b", "miR-122a", "miR-123", "miR-124a", "miR-125a", "miR-125b",
  "miR-140s", "miR-145", "miR-146", "miR-191", "miR-192", "miR-219",
  "miR-222", "miR-223", "miR-294"
) %>% paste0("hsa-", .)


# our study of fetal lung
`miRNAs_FetalLung (all)` <-
  FinalResults_IUS$miRNA
`miRNAs_FetalLung (male)` <-
  FinalResults_Male$miRNA
`miRNAs_FetalLung (female)` <-
  FinalResults_Female$miRNA




# fxn to count overlaps & look for family matches  -----------------------------

# looks through all-samples, male, female lists of IUS-miRNAs
check_list_overlap <-
  function(list_in) {
    
    print("# Tested in Fetal Lung")
    sum(list_in %in%
          (DEresults_IUS %>% filter(!is.na(pvalue)) %>% .$row)) %>% print()
  
    print("Signif in All-Samples")
    sum(list_in %in% (FinalResults_IUS$miRNA)) %>% print()
    FinalResults_IUS %>%
      filter(miRNA %in% list_in) %>%
      print()
    FinalResults_IUS %>%
      filter(miRNA %in% list_in) %>%
      .$miRNA %>%
      sort() %>%
      paste0(collapse = "; ") %>%
      print()
  
    print("Signif in Male Samples")
    sum(list_in %in% (FinalResults_Male$miRNA)) %>% print()
    FinalResults_Male %>%
      filter(miRNA %in% list_in) %>%
      print()
    FinalResults_Male %>%
      filter(miRNA %in% list_in) %>%
      .$miRNA %>%
      sort() %>%
      paste0(collapse = "; ") %>%
      print()
  
    print("Signif in Female Samples")
    list_in[list_in %in% (FinalResults_Female$miRNA)] %>% print()
    FinalResults_Female %>%
      filter(miRNA %in% list_in) %>%
      .$miRNA %>%
      sort() %>%
      paste0(collapse = "; ") %>%
      print()
    
    NULL
}

# returns just the number e.g., hsa-30-5p --> 30
cleanup_vector <- "hsa-miR-|hsa-let-|-5p|-3p|-1-|-2-|[a-z]"

make_flexfamily_list <-
  function(list_in) {
    list_in <- list_in %>%
      map_chr(~ gsub(cleanup_vector, "", .x)) %>%
      map_chr(~ str_split(.x, "-", simplify = T) %>% .[1]) %>%
      .[`!=`(., "")]

    list_out <-
      DEresults_IUS %>%
      mutate(numbers_only = gsub(cleanup_vector, "", row)) %>%
      filter(numbers_only %in% list_in) %>%
      arrange(qval) %>%
      .$row %>%
      unique()

    return(list_out)
  }


# loop through all studies  ----------------------------------------------------

article_list_names <- ls() %>% .[grepl("^miRNAs_", .)]

fetallung_list_names <-
  c(
    "miRNAs_FetalLung (all)",
    "miRNAs_FetalLung (male)",
    "miRNAs_FetalLung (female)"
  )

article_lists %>% .[setdiff(names(.), fetallung_list_names)] %>%
  purrr::map(.f = check_list_overlap)




# check #, fraction of overlaps pairwise for all studies -----------------------

article_lists <-
  article_list_names %>%
  purrr::map(~ get(x = .x)) %>%
  set_names(article_list_names)

article_lists_flex <-
  article_list_names %>%
  purrr::map(~ get(x = .x) %>%
    map_chr(~ gsub(cleanup_vector, "", .x)) %>%
    map_chr(~ str_split(.x, "-", simplify = T) %>% .[1]) %>%
    .[`!=`(., "")] %>%
    unique()) %>%
  set_names(article_list_names)

article_comparisons <-
  expand_grid(
    Paper1 = article_list_names,
    Paper2 = article_list_names
  )

article_comparisons$NumOverlap <-
  article_comparisons %>%
  pmap_dbl(.f = function(Paper1, Paper2, ...) {
    intersect(
      article_lists[Paper1] %>% unlist(),
      article_lists[Paper2] %>% unlist()
    ) %>% length()
  })

article_comparisons$NumOverlapFlex <-
  article_comparisons %>%
  pmap_dbl(.f = function(Paper1, Paper2, ...) {
    intersect(
      article_lists_flex[Paper1] %>% unlist(),
      article_lists_flex[Paper2] %>% unlist()
    ) %>% length()
  })

article_comparisons$FracOverlap <-
  article_comparisons %>%
  pmap_dbl(.f = function(Paper1, Paper2, ...) {
    intersect(
      article_lists[Paper1] %>% unlist(),
      article_lists[Paper2] %>% unlist()
    ) %>%
      length() %>%
      `/`(., article_lists[Paper1] %>% unlist() %>% length())
  })

article_comparisons$FracOverlapFlex <-
  expand_grid(
    Paper1 = article_list_names,
    Paper2 = article_list_names
  ) %>%
  pmap_dbl(.f = function(Paper1, Paper2, ...) {
    intersect(
      article_lists_flex[Paper1] %>% unlist(),
      article_lists_flex[Paper2] %>% unlist()
    ) %>%
      length() %>%
      `/`(., article_lists_flex[Paper1] %>% unlist() %>% length())
  })

article_comparisons <-
  article_comparisons %>%
  mutate_at(1:2, function(x) {
    gsub("miRNAs_", "", x)
  }) %>%
  arrange(grepl("Fetal", Paper1), grepl("Fetal", Paper2)) %>%
  mutate(
    Paper1 = fct_inorder(Paper1),
    Paper2 = fct_inorder(Paper2)
  ) %>%
  mutate_at(
    c("NumOverlap", "NumOverlapFlex"),
    function(x) {
      if_else(x == 0, NA_real_, x)
    }
  ) %>%
  mutate_at(c("FracOverlap", "FracOverlapFlex"), function(x) {
    replace_na(x, 0)
  })




# represent overlaps in a figure -----------------------------------------------

# note: put paper 2 as axis X so the fraction displayed on geom_tile
# is based on rows as the denominator (may be more intuitive?)

FigS1A <-
  ggplot(
    data = article_comparisons,
    aes(Paper2, Paper1, fill = FracOverlap)
  ) +
  geom_tile(color = "black") +
  geom_text(aes(label = NumOverlap)) +
  scale_fill_gradient2("Frac", low = "white", high = "maroon") +
  theme_few() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("") +
  ylab("") +
  ggtitle("Fraction Overlap",
          "# shared miRNAs / # miRNA reported in row")

FigS1B <-
  ggplot(
    data = article_comparisons,
    aes(Paper2, Paper1, fill = FracOverlapFlex)
  ) +
  geom_tile(color = "black") +
  geom_text(aes(label = NumOverlapFlex)) +
  scale_fill_gradient2("Frac", low = "white", high = "maroon") +
  theme_few() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("") +
  ylab("") +
  ggtitle("Fraction Overlap (Flexible/Family Match)",
          "# shared miRNAs / # miRNA reported in row")

# note that this is marked "New FigS1" for added Suppl. File
# old DIABLO supplemental figures still marked as "Fig1"
# (didn't want to change code / check still reproducible)

NewFigS1 <-
  plot_grid(FigS1A + theme(legend.position = "none"),
    FigS1B,
    rel_widths = c(4, 5),
    labels = c("A", "B")
  )

ggsave(NewFigS1, filename = "NewFigS1.png", width = 10, height = 5)
system('open "NewFigS1.png"')




# print names of miRNAs --------------------------------------------------------

# stringent
found_in_lit <-
  intersect(
    article_lists[fetallung_list_names] %>% unlist(),
    article_lists %>% .[setdiff(names(.), fetallung_list_names)] %>% unlist()
  )

found_in_lit[
  purrr::map_dbl(found_in_lit, ~ gsub(cleanup_vector, "", .x) %>%
    as.numeric()) %>% order()
] %>%
  paste0(collapse = ", ")

# flex
found_in_lit_flex <-
  intersect(
    article_lists_flex[fetallung_list_names] %>% unlist(),
    article_lists_flex %>% .[setdiff(names(.), fetallung_list_names)] %>% unlist()
  )

found_in_lit_flex %>%
  as.numeric() %>%
  sort() %>%
  paste0("miR-", .) %>%
  paste0(collapse = ", ")

# most common miRNAs
c(
  article_lists_flex[fetallung_list_names] %>% unlist() %>% unique(),
  article_lists_flex %>% .[setdiff(names(.), fetallung_list_names)] %>% unlist()
) %>%
  table() %>%
  .[`>=`(., 3)] %>%
  names() %>%
  as.numeric() %>%
  sort() %>%
  paste0("miR-", .) %>%
  paste(collapse = ", ")


