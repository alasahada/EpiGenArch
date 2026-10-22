# Synthetic data for figure generation

These CSV files are synthetic demonstration inputs supplied for the plotting examples. They are not manuscript result estimates and contain no participant-level records. Do not use them to interpret the study findings.

From the repository root, run:

```sh
Rscript 9_figures.R data synthetic_figures
```

The script requires the R packages `ggplot2`, `dplyr`, `patchwork`, and `tibble`. It writes demonstration figures to `synthetic_figures/`; the plots are labelled as synthetic examples.
