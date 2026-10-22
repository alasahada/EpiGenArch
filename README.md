# Distinct polygenic architecture of neurodevelopmental disorders with and without epilepsy

## Status and scope

The scripts illustrate the principal models and analysis steps. They are not a complete, validated reproduction of the manuscript results or an automated pipeline. Cohort construction, data extraction, genotype and phenotype quality control, and several harmonization steps are not included. Some analyses are illustrated with examples rather than by every comparison reported in the paper.

Analyses of the 100,000 Genomes Project were performed inside the secure Genomics England Research Environment. This repository contains code only; no individual-level genotype, phenotype, family, or polygenic score records are included. Reproducing analyses involving those records requires approved access to the relevant data in an authorised environment. Data access information is available from [Genomics England](https://www.genomicsengland.co.uk/research).

## Using the scripts

Read each script's input comments and the manuscript Methods before adapting it. Numbering groups analyses by topic; it does not mean that running scripts 1-9 in order creates all required inputs. Several scripts require independently prepared analysis tables or external software scripts and reference panels. Results are variously kept in R objects, printed, or written by external tools; there is no unified output directory or results export step.

The R scripts use `dplyr`, `purrr`, `broom`, `rcompanion`, `linemodels`, `TwoSampleMR`, and `MRPRESSO`, as applicable. Shell scripts require Bash and the relevant installations of PLINK2, SAIGE, LDSC, and GCTA; script 7 invokes PLINK for LD clumping. These tools and their reference files are not bundled. 

## Polygenic scores and association with case status

`1_calculate_pgs.sh`: calculates the three epilepsy and three cognitive polygenic scores with PLINK2.

`2_pgs_association.R`: standardizes the scores and tests their association with NDD-only and NDD-Epi case status (overall, and by solved status, rare-variant burden and acquired epilepsy risk factors) by logistic regression. It includes an example score × solved-status interaction for PGS_EPI and a DDD association example; DDD harmonization and subgroup preparation are not implemented here.

`3_variance_explained.R`: illustrates incremental Nagelkerke R² and conversion to the liability scale, using PGS_EPI in unsolved NDD-Epi versus cancer controls as an example.

## Polygenic transmission disequilibrium test (pTDT)

`4_ptdt.R`: tests whether affected probands carry more risk alleles than expected from their unaffected parents.

## Genome-wide analyses

`5_gwas_saige.sh`: runs the GWAS of NDD-Epi and NDD-only with SAIGE.

`6_heritability_genetic_correlation.sh`: estimates SNP heritability (LDSC, GCTA-GREML) and genetic correlations between epilepsy and cognitive traits (LDSC).

`7_locus_concordance.R`: compares effect sizes at independent loci between NDD-Epi and NDD-only, and with external GWAS of epilepsy and cognitive traits.

`8_linemodels.R`: classifies lead SNPs of epilepsy and cognitive-trait GWAS into epilepsy-specific, cognition-specific and shared effects (linemodels).

`9_figures.R`: Main figures code
