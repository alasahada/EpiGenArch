# Locus-level concordance of effect sizes
# Input (placeholder files): summary statistics with columns SNP, A1 (effect allele), A2, BETA, P, MAF

library(dplyr)

read_gwas <- function(f) read.delim(f)

# LD clumping with PLINK (1000 Genomes Phase 3 European reference): r2 < 0.1, 500 kb
clump <- function(gwas, p_threshold, r2 = 0.1, kb = 500) {
  gwas <- filter(gwas, P < p_threshold)
  write.table(gwas[, c("SNP", "P")], "clump_input.txt", quote = FALSE, row.names = FALSE)
  system(paste("plink --bfile 1000G_EUR --clump clump_input.txt",
               "--clump-p1", p_threshold, "--clump-p2", p_threshold,
               "--clump-r2", r2, "--clump-kb", kb, "--out clump_output"))
  read.table("clump_output.clumped", header = TRUE)$SNP
}

# Match alleles between two GWAS: flip the effect where alleles are swapped;
# drop palindromic (A/T, C/G) SNPs and allele mismatches
harmonise <- function(a, b) {
  m <- inner_join(a, b, by = "SNP", suffix = c(".a", ".b"))
  palindromic <- paste0(m$A1.a, m$A2.a) %in% c("AT", "TA", "CG", "GC")
  same <- m$A1.a == m$A1.b & m$A2.a == m$A2.b
  swap <- m$A1.a == m$A2.b & m$A2.a == m$A1.b
  m <- m[!palindromic & (same | swap), ]
  m$BETA.b <- ifelse(m$A1.a == m$A1.b, m$BETA.b, -m$BETA.b)
  m
}

# Pearson correlation of effect sizes and two-sided binomial test of sign concordance
concordance <- function(x, y) {
  k <- sum(sign(x) == sign(y))
  data.frame(n_loci = length(x), r = cor(x, y), p_r = cor.test(x, y)$p.value,
             sign_concordance = k / length(x), p_sign = binom.test(k, length(x), p = 0.5)$p.value)
}

epi  <- read_gwas("gwas/ndd_epi.txt")
only <- read_gwas("gwas/ndd_only.txt")

# (a) suggestive SNPs (P < 5e-6, MAF >= 0.05) clumped in each GWAS separately; union of lead SNPs
leads <- union(clump(filter(epi,  MAF >= 0.05), 5e-6),
               clump(filter(only, MAF >= 0.05), 5e-6))
h <- harmonise(filter(epi, SNP %in% leads), only)
concordance(h$BETA.a, h$BETA.b)

# (b) independent lead SNPs of each external GWAS (P < 5e-8; P < 1e-6 for non-acquired focal
#     epilepsy, which has no genome-wide significant SNP), compared with the NDD GWAS
external <- list(
  educational_attainment = list(file = "external/educational_attainment.txt", p = 5e-8),
  cognitive_function     = list(file = "external/cognitive_function.txt",     p = 5e-8),
  epilepsy_all           = list(file = "external/ilae_all.txt",               p = 5e-8),
  epilepsy_nafe          = list(file = "external/ilae_nafe.txt",              p = 1e-6),
  epilepsy_gge           = list(file = "external/ilae_gge.txt",               p = 5e-8)
)
results <- bind_rows(lapply(names(external), function(trait) {
  ext   <- read_gwas(external[[trait]]$file)
  leads <- clump(ext, external[[trait]]$p)
  bind_rows(lapply(list("NDD-Epi" = epi, "NDD-only" = only), function(ndd) {
    h <- harmonise(filter(ext, SNP %in% leads), ndd)
    concordance(h$BETA.a, h$BETA.b)
  }), .id = "ndd_gwas") %>% mutate(external_trait = trait)
}))

# Sensitivity analysis without LD clumping: all genome-wide significant ILAE all-epilepsy SNPs
ilae <- read_gwas("external/ilae_all.txt")
h <- harmonise(filter(ilae, P < 5e-8), epi)
concordance(h$BETA.a, h$BETA.b)
