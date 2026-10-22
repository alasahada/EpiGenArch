#!/bin/bash
# SNP heritability and genetic correlation. File names are placeholders.
# GWAS summary statistics have columns SNP, A1, A2, BETA, P.

# 1. SNP heritability from the case-control GWAS (LD score regression),
#    liability scale assuming a population prevalence of 0.007
ldsc_h2() {    # arguments: trait, number of cases, number of controls
    munge_sumstats.py --sumstats gwas/$1.txt --N-cas $2 --N-con $3 \
                      --merge-alleles w_hm3.snplist --out ldsc/$1
    ldsc.py --h2 ldsc/$1.sumstats.gz --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ \
            --samp-prev $(awk -v a=$2 -v b=$3 'BEGIN{print a/(a+b)}') --pop-prev 0.007 \
            --out ldsc/$1_h2
}
ldsc_h2 ndd_epi  2351 9303
ldsc_h2 ndd_only 4983 9303

# 2. SNP heritability from individual-level data (GCTA-GREML)
#    GRM from a high-confidence, LD-pruned panel of common SNPs (MAF >= 5%) in unrelated
#    European-ancestry individuals; covariates age, sex and PC1-PC10;
#    liability scale with prevalence 0.0015 (NDD-Epi) and 0.0085 (NDD-only)
gcta64 --bfile genotypes/hq_snps --keep samples/unrelated_european.txt --maf 0.05 \
       --make-grm --out grm/common_snps
gcta64 --reml --grm grm/common_snps --pheno phenotypes/ndd_epi.txt --prevalence 0.0015 \
       --covar covariates_discrete.txt --qcovar covariates_quantitative.txt --out greml/ndd_epi
gcta64 --reml --grm grm/common_snps --pheno phenotypes/ndd_only.txt --prevalence 0.0085 \
       --covar covariates_discrete.txt --qcovar covariates_quantitative.txt --out greml/ndd_only

# 3. Genetic correlation between epilepsy and cognitive traits (public GWAS, each file
#    prepared with munge_sumstats.py as above): three ILAE epilepsy GWAS (all epilepsy,
#    non-acquired focal, genetic generalised) against intelligence, educational attainment
#    and cognitive performance
for EPI in epilepsy_all epilepsy_nafe epilepsy_gge; do
    ldsc.py --rg ldsc/${EPI}.sumstats.gz,ldsc/intelligence.sumstats.gz,ldsc/educational_attainment.sumstats.gz,ldsc/cognitive_performance.sumstats.gz \
            --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out ldsc/rg_${EPI}
done
