# Variance explained by a polygenic score, on the liability scale
#
# Increase in Nagelkerke R2 between a covariates-only model and a model that also includes
# the score (rcompanion), converted to the liability scale (Lee et al. 2012).
#
# Input (placeholder file): one row per individual with case (0/1), the score, age, sex, PC1-PC10.

library(rcompanion)
library(dplyr)

covars <- c("age", "sex", paste0("PC", 1:10))

nagelkerke_r2 <- function(fit) {
  nagelkerke(fit)$Pseudo.R.squared.for.model.vs.null["Nagelkerke (Cragg and Uhler)", 1]
}

# K: population prevalence; P: proportion of cases in the sample
liability_r2 <- function(r2, K, P) {
  z <- dnorm(qnorm(1 - K))
  r2 * (K * (1 - K) / z^2) * (K * (1 - K) / (P * (1 - P)))
}

variance_explained <- function(dat, score, K) {
  base <- glm(reformulate(covars, "case"), data = dat, family = binomial)
  full <- glm(reformulate(c(score, covars), "case"), data = dat, family = binomial)
  delta <- nagelkerke_r2(full) - nagelkerke_r2(base)
  tibble(score = score, r2_observed = delta,
         r2_liability = liability_r2(delta, K, P = mean(dat$case)))
}

# Population prevalences assumed for each comparison
#   NDD-only vs controls           K = 0.0085
#   NDD-Epi vs controls            K = 0.0015   (unsolved NDD-Epi: 0.001; solved: 0.0005)
#   NDD-Epi vs NDD-only            K = 0.10
# Example: unsolved NDD-Epi versus cancer controls
dat <- read.delim("unsolved_ndd_epi_vs_cancer_controls.tsv")
variance_explained(dat, "PGS_EPI", K = 0.001)
