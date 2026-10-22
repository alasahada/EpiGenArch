# Line-based classification of epilepsy and cognitive-trait lead SNPs (linemodels; Pirinen 2023)
# Input (placeholder files): for each epilepsy x cognitive-trait pair, the independent lead SNPs
# (variants reaching a 5% FDR in either trait, LD-clumped with PLINK at r2 < 0.01 within 1 Mb in
# 1000 Genomes Europeans) with columns SNP, MAF, beta_epi, se_epi, beta_cog, se_cog (effects aligned
# to the epilepsy effect allele), and the LD score regression genetic correlation (rg) of each pair.

library(linemodels)
set.seed(1)    # placeholder seed

classify_pair <- function(snps, rg, n_sims = 200) {
  # put the binary and continuous traits on a common per-allele scale
  sc <- sqrt(2 * snps$MAF * (1 - snps$MAF))
  X  <- cbind(snps$beta_epi, snps$beta_cog) * sc
  SE <- cbind(snps$se_epi,   snps$se_cog)   * sc

  # three line models in the (epilepsy, cognition) plane: epilepsy-specific (slope 0),
  # cognition-specific (vertical) and shared (slope fixed from sign(rg) x effect scale ratio)
  model_names <- c("epilepsy", "cognition", "shared")
  vertical <- 1e6                       # vertical line (infinite slope) approximated by a very large slope
  slopes <- c(0, vertical, sign(rg) * sd(X[, 2]) / sd(X[, 1]))
  cors   <- rep(0.99, 3)

  # scale of the two trait-specific models by maximum likelihood; the shared model
  # takes the scale of the cognition model
  fit <- line.models.optimize(X, SE,
           par.include = rbind(c(TRUE, FALSE, FALSE), c(TRUE, FALSE, FALSE)),
           init.scales = c(0.05, 0.05), init.slopes = slopes[1:2], init.cors = cors[1:2],
           model.names = model_names[1:2], r.lkhood = 0, assume.constant.SE = TRUE,
           print.steps = 0)
  scales <- c(fit$scales, fit$scales[2])

  # posterior membership; a SNP is assigned to a model if its probability exceeds 0.95
  post <- line.models(X, SE, scales = scales, slopes = slopes, cors = cors,
                      model.names = model_names, model.priors = rep(1, 3), r.lkhood = 0)
  colnames(post) <- model_names
  class <- ifelse(apply(post, 1, max) > 0.95, model_names[apply(post, 1, which.max)], "undetermined")

  # evidence for a shared class: log-likelihood ratio of the three-model mixture against the
  # mixture without the shared line, null distribution by simulation
  sim <- simulate.loglr(X, SE, n.sims = n_sims,
           par.include.null = rbind(c(TRUE, FALSE, FALSE), c(TRUE, FALSE, FALSE)),
           init.scales.null = scales[1:2], init.slopes.null = slopes[1:2], init.cors.null = cors[1:2],
           par.include.alt  = rbind(c(TRUE, FALSE, FALSE), c(TRUE, FALSE, FALSE), c(FALSE, FALSE, FALSE)),
           init.scales.alt  = scales, init.slopes.alt = slopes, init.cors.alt = cors,
           r.lkhood = 0, assume.constant.SE = TRUE, print.steps = c(0, 0))
  list(classes = table(class),
       p_shared = (1 + sum(sim$null.loglr >= sim$obs.loglr)) / (n_sims + 1))
}

# six pairs: two epilepsy phenotypes x three cognitive traits
pairs <- expand.grid(epilepsy  = c("epilepsy_all", "epilepsy_nafe"),
                     cognition = c("cognitive_performance", "educational_attainment", "intelligence"),
                     stringsAsFactors = FALSE)
rg_table <- read.csv("ldsc_rg.csv")      # columns epilepsy, cognition, rg

results <- lapply(seq_len(nrow(pairs)), function(i) {
  snps <- read.delim(sprintf("lead_snps/%s_%s.txt", pairs$epilepsy[i], pairs$cognition[i]))
  rg   <- with(rg_table, rg[epilepsy == pairs$epilepsy[i] & cognition == pairs$cognition[i]])
  classify_pair(snps, rg)
})
