#!/bin/bash
# Calculate polygenic scores with PLINK2.
#
# Six scores, each with one weight file (variant ID, effect allele, weight):
#   epilepsy:  PGS_EPI (all epilepsy), PGS_NAFE (non-acquired focal), PGS_GGE (genetic generalised);
#              ILAE GWAS, weights from Heyne et al. 2024
#   cognitive: PGS_IQ (PGS003724), PGS_CP (PGS002135), PGS_EA (PGS002231); PGS Catalog
#
# File names are placeholders (the genotypes can only be accessed inside the
# Genomics England Research Environment).

GENOTYPES=genotypes/aggV2_qc     # quality-controlled genotypes (PLINK2 pfile)

for SCORE in PGS_EPI PGS_NAFE PGS_GGE PGS_IQ PGS_CP PGS_EA; do
    plink2 --pfile ${GENOTYPES} \
           --score weights/${SCORE}.txt 1 2 3 header cols=+scoresums \
           --out scores/${SCORE}
done

# Overlap between the ILAE discovery GWAS and the GEL sample was assessed with EraSOR
# (Choi et al. 2023). Scores are standardised in 2_pgs_association.R.
