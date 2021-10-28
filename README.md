# Prenatal metal exposure, cord blood DNA methylation and persistence in childhood: an epigenome-wide association study of twelve metals

This repository contains the necessary scripts used in the analysis of Bozack et al.'s "Prenatal metal exposure, cord blood DNA methylation and persistence in childhood: an epigenome-wide association study of twelve metals".

The organization of the repository is as follows:

- <b>EWAS functions.R:</b> functions including plotting and wrapper functions for analyses of DMPs and DMRs /b
Note: This scrpt also contains a modified version of the combp function to exclude DMRs consisting of one probe
- <b>viva_DNAm_metals_dataProcessing_descriptiveStats.rmd:</b> Processing and descriptive statistics for cord blood DNAm data
- <b>viva_DNAm_metals_dataProcessing_descriptiveStats.rmd:</b> Processing and descriptive statistics for mid-childhood DNAm data
- <b>viva_DNAm_metals_cordblood_DMPs_DMRs_smk.rmd:</b> EWAS and DMR analysis for cord blood DNAm data
- <b>viva_DNAm_metals_age7_DMPs_DMRs_smk.rmd:</b> EWAS and DMR analysis for mid-childhood DNAm data
- <b>EWAS_results folder:</b> Complete output of cord blood adjusted EWAS for all metals 
