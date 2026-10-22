# Polygenic transmission disequilibrium test (pTDT; Weiner et al. 2017)
#
# Do affected probands carry more (or fewer) risk alleles than expected from their unaffected
# parents? Each proband's score is compared with the mean of the two parental scores
# (mid-parent score); the deviation is tested against zero.
#
# Input (placeholder file): one row per parent-offspring trio with
#   child_<score>, mother_<score>, father_<score> for each score,
#   group (NDD_only | NDD_Epi), solved, asd (0/1), panel_ndd / panel_epilepsy / panel_ion_channel (0/1)

library(dplyr)
library(purrr)

trios     <- read.delim("trio_table.tsv")
scores    <- c("PGS_EPI", "PGS_NAFE", "PGS_GGE", "PGS_IQ", "PGS_CP", "PGS_EA")
cognitive <- c("PGS_IQ", "PGS_CP", "PGS_EA")

mid_parent <- function(d, s) (d[[paste0("mother_", s)]] + d[[paste0("father_", s)]]) / 2

# SD of the mid-parent score, calculated once over all trios so that deviations are
# on the same scale in every subgroup
sd_mid <- map_dbl(set_names(scores), ~ sd(mid_parent(trios, .x)))

ptdt <- function(d, s) {
  deviation <- (d[[paste0("child_", s)]] - mid_parent(d, s)) / sd_mid[[s]]
  t <- t.test(deviation, mu = 0)                      # two-sided one-sample t-test
  tibble(score = s, n_trios = length(deviation), mean_deviation = unname(t$estimate),
         lower = t$conf.int[1], upper = t$conf.int[2], p = t$p.value)
}

run <- function(groups, scs) {
  imap_dfr(groups, function(d, name) map_dfr(scs, ~ ptdt(d, .x)) %>% mutate(subgroup = name))
}

subgroups <- list(
  "NDD-Epi, unsolved"                = filter(trios, group == "NDD_Epi", solved == "unsolved"),
  "NDD-Epi, solved"                  = filter(trios, group == "NDD_Epi", solved == "solved"),
  "NDD-Epi, epilepsy-gene carriers"  = filter(trios, group == "NDD_Epi", panel_epilepsy == 1),
  "NDD-Epi, NDD-gene carriers"       = filter(trios, group == "NDD_Epi", panel_ndd == 1),
  "NDD-Epi, ion-channel-gene carriers" = filter(trios, group == "NDD_Epi", panel_ion_channel == 1),
  "NDD-only"                         = filter(trios, group == "NDD_only")
)

# Cognitive scores are positively associated with autism, so these analyses are repeated
# without probands with an ASD annotation (HPO HP:0000717, HP:0000729 and descendants)
no_asd <- map(subgroups, ~ filter(.x, asd == 0))
names(no_asd) <- paste(names(subgroups), "(no ASD)")

results <- bind_rows(run(subgroups, scores), run(no_asd, cognitive)) %>%
  mutate(bonferroni_significant = p < 0.05 / n())     # Bonferroni across scores and subgroups
