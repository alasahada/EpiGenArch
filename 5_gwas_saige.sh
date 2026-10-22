#!/bin/bash
# Genome-wide association studies with SAIGE (logistic mixed model that accounts for
# relatedness and case-control imbalance).
#
# Two GWAS, each against its own half of the control pool (cancer and non-neurological
# rare-disease participants, split into two non-overlapping halves stratified by recruitment
# source, fixed random seed) so that the two scans are independent:
#   ndd_epi   NDD-Epi vs controls  (a sensitivity GWAS uses unsolved NDD-Epi only)
#   ndd_only  NDD-only vs controls
# European ancestry; covariates age, sex and PC1-PC10; variants with MAF >= 1%.
# File names are placeholders; options follow the SAIGE documentation.

COVARS=age,sex,$(seq -s, -f 'PC%g' 1 10)

for TRAIT in ndd_epi ndd_only; do

    # Step 1: fit the null model
    Rscript step1_fitNULLGLMM.R \
        --plinkFile=genotypes/grm_snps \
        --phenoFile=phenotypes/${TRAIT}.txt --phenoCol=case --sampleIDColinphenoFile=IID \
        --covarColList=${COVARS} --qCovarColList=sex \
        --traitType=binary --outputPrefix=saige/${TRAIT}_null

    # Step 2: test variants with MAF >= 1%, chromosome by chromosome
    for CHR in {1..22}; do
        Rscript step2_SPAtests.R \
            --bgenFile=genotypes/chr${CHR}.bgen --bgenFileIndex=genotypes/chr${CHR}.bgen.bgi \
            --sampleFile=phenotypes/${TRAIT}_samples.txt \
            --chrom=${CHR} --minMAF=0.01 --minMAC=1 \
            --GMMATmodelFile=saige/${TRAIT}_null.rda \
            --varianceRatioFile=saige/${TRAIT}_null.varianceRatio.txt \
            --SAIGEOutputFile=saige/${TRAIT}_chr${CHR}.txt
    done
done

# Genomic inflation (lambda GC), Manhattan and QQ plots are made from the step 2 output.
