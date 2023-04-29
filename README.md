## Intrauterine Smoke Exposure, microRNA Expression during Human Lung Development, and Childhood Asthma

### Research Summary

Intrauterine smoke exposure (IUS) is associated with increased post-natal asthma risk and immune dysregulation, but the underlying pathophysiology is poorly understood. Given the role of miRNAs in orchestrating gene regulatory programs during lung development, we hypothesized that IUS alters miRNA profiles during prenatal lung development.

Using paired miRNA and mRNA profiles from n = 298 human fetal lung samples, we demonstrate that: **(i)** IUS is associated with differential expression of miRNAs in early lung development. These fetal lung IUS-miRNAs are plausibly associated with the developmental origins of allergy and airway disease, based on their **(ii)** correlation with mRNA genes relevant to lung development and asthma/immunity within fetal lung (e.g., ORDML3, oxidative stress/mGLuR genes) and **(iii)** correlation to post-natal outcomes such as IgE and FEV1/FVC within serum miRNA profiles of a childhood asthmatics cohort (GACRS). We also find evidence for **(iv)** sex-specific impacts of IUS on fetal lung, with more statistically significant IUS-associated miRNAs in male samples than female samples and distinct miRNA-mRNA networks by sex.

### Files in Repo

`Data`: This repo includes processed data (miRNA counts; microarray normalized intensities) and sample metadata (IUS status, age) to reproduce the fetal lung analyses in the manuscript. **See GEO links below for raw data.**

`Scripts`: We also include code for the IUS differential expression analysis, IUS-by-sex interaction analysis, and multi-'omic DIABLO integration, starting from processed miRNA & mRNA data. We also performed a review of previously reported miRNA differential expression due to IUS/smoking (supplementary files).

```
File Structure
├─ Data
│   └─ miRNA_counts.tsv            # filtered counts, 898 miRNAs x 298 samples
│   └─ mRNA_norm.tsv               # log2, RMA-normalized Affy intensities,  18178 probesets x 298 samples
│   └─ sample_metadata.tsv         # 298 samples x 6 metadata variables
└─ Scripts
│   └─ 00_setup.R                  # load data, basic summary/EDA
│   └─ 01_load_dat_and_annot.R 
│   └─ 02_table1.R
│   └─ 03_primary_model.R          # test IUS effects
│   └─ 04_interaction_model.R      # IUS-by-sex effects
│   └─ 05_run_miEAA_GSEA.R         # pathway analyses
│   └─ 06_DIABLO_male.R            # miRNA-mRNA networks
│   └─ 07_DIABLO_female.R 
│   └─ 08_DIABLO_misc.R
│   └─ 09_export_results.R         # tidy results + pairwise correlations
│   └─ 10_compare_smoking_sig.R    # smoke exposure literature review 
└─ README.md                       # this readme
└─ sessionInfo.txt                 # version control
```

**Not shown here:**

* The code for miRNA-seq mapping, QC, and quantification (i.e., .fastq --> counts) are identical to that of our first publication on the miRNA data on sex-effects in fetal lung, shown in scripts 01 through 17 at the following Github repo: [\@chooliu/miRNASexDimorphismFetalLung](https://github.com/chooliu/miRNASexDimorphismFetalLung). 
* Details of fetal lung gene expression profiling are described in the Kho, et al. manuscript and associated GEO link. The expression data provided in the Data folder were subsetted from a larger set of samples and probesets processed by Dr. Alvin T. Kho.
* miRNA analysis in childhood asthma dataset (under preparation by Channing/Harvard collaborators; Genetics of Asthma in Costa Rica Study [GACRS])


### Key Links

**Manuscript & Citation**

Rosenberg L, Liu C, Sharma R, Wood C, Vyhlidal CA, Gaedigk R, Kho AT, Ziniti JP, Celedón JC, Tantisira KG, Weiss ST, McGeachie MJ, Kechris K, Sharma S. Intrauterine Smoke Exposure, microRNA Expression during Human Lung Development, and Childhood Asthma. International Journal of Molecular Sciences. 2023; 24(9):7727. 

* This Github Repo: [chooliu/FetalLung_miRNA_SmokeExposure](http://www.github.com/chooliu/FetalLung_miRNA_SmokeExposure) (analysis code)
* Manuscript: [original article](https://www.mdpi.com/1422-0067/24/9/7727) (open access)
* doi: [10.3390/ijms24097727](doi.org/10.3390/ijms24097727)

**Raw Data, Hosted with NCBI GEO**

* Fetal Lung miRNA: [GSE200153](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE200153) (small RNA-seq; .fastq.gz)
* Fetal Lung mRNA: [GSE68896](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE68896) (Affy microarray from Kho, 2015; .CEL.gz)
* Childhood Asthmatics miRNA: TBD (GARCS miRNA under preparation)

**Related Fetal Lung Manuscripts**

This manuscript is part of a continued multi-institution effort profiling the genomics of early lung development, including the multi-'omic impacts of tobacco smoke exposure on human fetal lung tissue.

* "Sex-Specific Differences in MicroRNA Expression During Human Fetal Lung Development" (Lin, 2022; Frontiers in Genetics) - [link](https://www.frontiersin.org/articles/10.3389/fgene.2022.762834/full)
* "Sex-specific associations with DNA methylation in lung tissue demonstrate smoking interactions" (Koo, 2021; Epigenetics) - [link](https://www.tandfonline.com/doi/full/10.1080/15592294.2020.1819662)
* "Age, Sexual Dimorphism, and Disease Associations in the Developing Human Fetal Lung Transcriptome" (Kho, et al., 2015; AJRCMB) - [link](https://www.tandfonline.com/doi/full/10.1080/15592294.2020.1819662)

