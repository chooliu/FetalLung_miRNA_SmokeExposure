### FetalLung_IUS_miRNA 


### Summary

**"The Impact of Intrauterine Smoke Exposure on MicroRNA Expression During Human Lung Development and Asthma and Allergy Outcomes in Childhood Asthma"** (submitted)

Intrauterine smoke exposure (IUS) is associated with increased post-natal asthma risk and immune dysregulation, but the underlying pathophysiology is poorly understood. Given the role of miRNAs in orchestrating gene regulatory programs during lung development, we hypothesized that IUS alters miRNA profiles during prenatal lung development.

Using paired miRNA and mRNA profiles from n = 298 whole human fetal lung samples, we demonstrate that: **(i)** IUS is associated with differential expression of miRNAs in early lung development and that **(ii)** these miRNAs are correlated with mRNA genes relevant to lung development and asthma/immunity. We also find evidence for **(iii)** sex-by-IUS interactions in fetal lung, with more statistically significant IUS-associated miRNAs in male samples than female samples and distinct miRNA-mRNA networks by sex. Lastly, within serum miRNA profiles within a childhood asthmatics cohort (GACRS), **(iv)** the fetal lung IUS-miRNAs are correlated to post-natal outcomes such as IgE and FEV1/FVC in a sex-specific manner, supporting our findings that IUS-miRNA impacts may be sex-specific and that IUS-miRNAs may be implicated in the developmental origins of allergy and airway disease.

### Files in Repo

`Data`: This repo includes processed data (miRNA counts; microarray normalized intensities) and sample metadata (IUS status, age) to reproduce the fetal lung analyses in the manuscript. See GEO links for raw data below.

`Scripts`: We also include code for the IUS differential expression analysis, IUS-by-sex interaction analysis, and DIABLO integration, starting from processed miRNA & mRNA data.

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
└─ README.md                       # this readme
└─ sessionInfo.txt                 # version control
```

Not shown here:

* The code for miRNA-seq mapping, QC, and quantification (i.e., .fastq --> counts) are identical to that of our first publication on the miRNA data, shown in scripts 01 through 17 at the following Github repo: [\@chooliu/miRNASexDimorphismFetalLung](https://github.com/chooliu/miRNASexDimorphismFetalLung). 
* Details of fetal lung gene expression profiling are described in the Kho, et al. manuscript and associated GEO link. The expression data provided in the Data folder were subsetted from a larger set of samples and probesets processed by Dr. Alvin T. Kho.
* miRNA analysis in childhood asthma dataset (under preparation by Channing/Harvard collaborators; Genetics of Asthma in Costa Rica Study [GACRS])


### Key Links

**Manuscript Links**

* Github Repo: [chooliu/FetalLung_miRNA_SmokeExposure](http://www.github.com/chooliu/FetalLung_miRNA_SmokeExposure)
* manuscript currently under review, link to final published work, doi, etc will be added

**Raw Data, Hosted with NCBI GEO**

* Fetal Lung miRNA: [GSE200153](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE200153) (small RNA-seq; .fastq.gz)
* Fetal Lung mRNA: [GSE68896](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE68896) (Affy microarray from Kho, 2015; .CEL.gz)
* Childhood Asthmatics miRNA: TBD (GARCS miRNA under preparation)

**Related Fetal Lung Manuscripts**

This manuscript is part of a continued multi-institution effort profiling the genomics of early lung development, including the multi-'omic impacts of tobacco smoke exposure on human fetal lung tissue.

* "Sex-Specific Differences in MicroRNA Expression During Human Fetal Lung Development" (Lin, 2022; Frontiers in Genetics) - [link](https://www.frontiersin.org/articles/10.3389/fgene.2022.762834/full)
* "Sex-specific associations with DNA methylation in lung tissue demonstrate smoking interactions" (Koo, 2021; Epigenetics) - [link](https://www.tandfonline.com/doi/full/10.1080/15592294.2020.1819662)
* "Age, Sexual Dimorphism, and Disease Associations in the Developing Human Fetal Lung Transcriptome" (Kho, et al., 2015; AJRCMB) - [link](https://www.tandfonline.com/doi/full/10.1080/15592294.2020.1819662)

