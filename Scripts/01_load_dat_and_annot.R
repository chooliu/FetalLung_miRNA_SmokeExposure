# ==============================================================================
# 01_load_dat_and_annot.R
# load 'omics data, 'omics annotations, sample metadata
# ==============================================================================




# load post-bioinformatics pipeline versions of data ---------------------------
# miRNA counts, mRNA log2intensities, & phenodata

# miRNA: mapping, QC metadata clean-up identical to that of our first 
# fetal lung miRNA publication; see scripts 01 through 17 @
# https://github.com/chooliu/miRNASexDimorphismFetalLung

# mRNA: data normalization performed by A.T.K. (microarray --> normalization)
# see Kho, et al. "Age, Sexual Dimorphism, and Disease Associations in the 
# Developing Human Fetal Lung Transcriptome." Am J Respir Cell Mol Biol. 2016

# can download contents of "Data" in the current repo @ 
# https://github.com/chooliu/FetalLung_miRNA_SmokeExposure
# and place in "Data" folder in working directory

metadata_miRNA <- read_tsv("Data/sample_metadata.tsv")
feature_miRNA <- read_tsv("Data/miRNA_counts.tsv") %>%
    column_to_rownames("miRNA") %>% as.matrix()
feature_microarray <- read_tsv("Data/microarray_norm.tsv") %>%
    column_to_rownames("probeset") %>% as.matrix()



# miRNA annotations ------------------------------------------------------------
# genomic position & sequence checking / filtering

# mature miRNA sequences & annotations from mirbase v22.1 ftp
# https://www.mirbase.org

mirna_sequences_hsa <-
  readLines("Data/mature.fa") %>%
  data.frame(
    orig_header = .[c(T, F)],
    sequence = .[c(F, T)]) %>%
  transmute(id = str_split(orig_header, " ") %>% map_chr(., ~ .[[1]]),
            accession = str_split(orig_header, " ") %>% map_chr(., ~ .[[2]]),
            sequence_dna = gsub("U", "T", sequence)) %>%
  filter(grepl("hsa", id)) %>%
  mutate(id = gsub(">", "", id))

# annotate miRNAs with chr location
# important because high density of miRNA on chrX
mirna_chrom_locations <-
  fread("Data/hsa.gff3") %>%
  separate(col = V9, sep = "Name\\=", into = c("discard", "tmp")) %>%
  separate(col = tmp, sep = ";", into = c("Name", "discard")) %>%
  transmute(miRNA = Name, Chr = gsub("chr", "", `#`) %>%
              str_sort(., numeric = T), Left = V4, Right = V5) %>%
  group_by(miRNA) %>% summarize(Chr = unique(Chr) %>% sort %>% paste(collapse = ","))

annotation_miRNAs <-
  tibble(miRNA = rownames(feature_miRNA)) %>%
  left_join(., mirna_chrom_locations, by = "miRNA")





# microarray annotations download (ENSEMBL/biomaRt) ----------------------------
# saved as .Rdata to avoid repeated queries; run commented code below first

# annotation_ensembl <-
#   useMart(biomart = "ENSEMBL_MART_ENSEMBL",
#           dataset = "hsapiens_gene_ensembl",
#           host = "sep2019.archive.ensembl.org") %>%
#   getBM(attributes =
#           c("ensembl_gene_id",
#             "ensembl_transcript_id",
#             "hgnc_symbol", "external_gene_name",
#             "description",
#             "chromosome_name", "band",
#             "transcription_start_site", "transcript_start", "transcript_end",
#             "strand",
#             "gene_biotype",
#             "affy_hugene_1_0_st_v1",
#             "entrezgene_id"),
#         mart = .)
# save(annotation_ensembl, file = "./Data/biomart_grch38_ensembl.Rdata")

load("./Data/biomart_grch38_ensembl.Rdata")

annotations_mRNA <-
  annotation_ensembl %>%
  transmute(name = as.character(affy_hugene_1_0_st_v1), external_gene_name) %>%
  group_by(name, external_gene_name) %>%
  group_by(name) %>%
  summarize(symbol = external_gene_name %>%
              unique %>% .[`!=`(., "")] %>%
              paste(collapse = "/"))

tidy_mRNA_names <-
  tibble(name = row.names(feature_microarray)) %>%
  left_join(., annotations_mRNA, by = "name") %>%
  .$symbol %>%
  toupper()

row.names(feature_microarray) <- tidy_mRNA_names




# double check names correct ---------------------------------------------------
# both should return TRUE

identical(colnames(feature_microarray),
          colnames(feature_miRNA))

identical(metadata_miRNA$Sample %>% as.character(),
          colnames(feature_miRNA))

