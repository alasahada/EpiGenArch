# Association of polygenic scores with case status
#
# Input (placeholder file): one row per individual with
#   group     NDD_only | NDD_Epi | control_cancer | control_rare
#   solved    solved | unsolved (molecular diagnostic status)
#   rv_group  unsolved NDD-Epi only: HC_pLoF | LC_pLoF_or_missense | none   (rare-variant burden)
#   exposure  unsolved NDD-Epi only: none | any                            (acquired epilepsy risk factors)
#   age, sex, PC1-PC10 and the six scores

library(dplyr)
library(purrr)
library(broom)

d      <- read.delim("analysis_table.tsv")
scores <- c("PGS_EPI", "PGS_NAFE", "PGS_GGE", "PGS_IQ", "PGS_CP", "PGS_EA")
covars <- c("age", "sex", paste0("PC", 1:10))

# 1. Standardise each score to mean 0 and SD 1 (odds ratios are then per SD)
d <- mutate(d, across(all_of(scores), ~ as.numeric(scale(.x))))

# 2. Logistic regression of case status on each score, adjusted for age, sex and 10 PCs
test_scores <- function(cases, controls, covars) {
  dat <- bind_rows(mutate(cases, case = 1), mutate(controls, case = 0))
  map_dfr(scores, function(s) {
    fit <- glm(reformulate(c(s, covars), "case"), data = dat, family = binomial)
    tidy(fit) %>%
      filter(term == s) %>%
      transmute(score = s, OR = exp(estimate),
                lower = exp(estimate - 1.96 * std.error),
                upper = exp(estimate + 1.96 * std.error),
                p = p.value, n_cases = nrow(cases), n_controls = nrow(controls))
  })
}

# 3. Cases versus cancer controls: NDD-only, NDD-Epi, solved and unsolved NDD-Epi, and
#    rare-variant and acquired-risk-factor strata within unsolved NDD-Epi
controls <- filter(d, group == "control_cancer")
epi      <- filter(d, group == "NDD_Epi")
unsolved <- filter(epi, solved == "unsolved")

comparisons <- list(
  "NDD-only"                              = filter(d, group == "NDD_only"),
  "NDD-Epi"                               = epi,
  "NDD-Epi, unsolved"                     = unsolved,
  "NDD-Epi, solved"                       = filter(epi, solved == "solved"),
  "unsolved, HC pLoF"                     = filter(unsolved, rv_group == "HC_pLoF"),
  "unsolved, LC pLoF or damaging missense" = filter(unsolved, rv_group == "LC_pLoF_or_missense"),
  "unsolved, no qualifying rare variant"  = filter(unsolved, rv_group == "none"),
  "unsolved, no acquired risk factor"     = filter(unsolved, exposure == "none"),
  "unsolved, one or more risk factors"    = filter(unsolved, exposure == "any")
)

results <- imap_dfr(comparisons, function(cases, name) {
  mutate(test_scores(cases, controls, covars), comparison = name)
}) %>%
  group_by(comparison) %>%
  mutate(bonferroni_significant = p < 0.05 / n()) %>%    # Bonferroni across the scores tested
  ungroup()

# 4. Does the effect of the score depend on molecular diagnostic status?
#    Joint model across NDD-Epi cases and the expanded control cohort (cancer and
#    non-neurological rare-disease controls) with a score x solved-status interaction
joint <- d %>%
  filter(group %in% c("NDD_Epi", "control_cancer", "control_rare")) %>%
  mutate(case   = as.integer(group == "NDD_Epi"),
         solved = factor(solved, levels = c("unsolved", "solved")))
fit_interaction <- glm(reformulate(c("PGS_EPI * solved", covars), "case"),
                       data = joint, family = binomial)
tidy(fit_interaction) %>% filter(grepl(":", term))

# 5. Replication in DDD: the same models, DDD probands against the same GEL cancer controls,
#    adjusted for sex and PC1-PC10 (age was not available for all DDD probands). DDD scores are
#    harmonised across cohorts and scaled to the GEL control distribution (see Methods).
ddd         <- read.delim("ddd_table.tsv")
ddd_results <- test_scores(ddd, controls, setdiff(covars, "age"))
