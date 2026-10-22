# DEMONSTRATION ONLY — Synthetic example results, not manuscript estimates.
# Optional: Rscript figures/plot_publication_figures.R [input_directory] [output_directory]

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)
input_dir <- if (length(args) >= 1L) args[1] else 'data'
out_dir <- if (length(args) >= 2L) args[2] else 'synthetic_figures'
read_example <- function(filename) {
  path <- file.path(input_dir, filename)
  if (!file.exists(path)) stop('Missing example input: ', path)
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

COL_GRAY <- '#888780'
COL_AMBER <- '#E9A400'
COL_GREEN <- '#147D64'
COL_BLUE <- '#4A90C4'
COL_TEAL <- '#1D9E75'

# Label only nominally significant results; asterisk indicates Bonferroni significance.
p_label <- function(p, n_tests) {
  ifelse(p < 0.05, paste0('p=', formatC(p, format = 'e', digits = 2),
                           ifelse(p < 0.05/n_tests, ' *', '')), NA_character_)
}

# Publication-friendly figure typography; units are points.
base_theme <- theme_minimal(base_size = 13) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(colour = 'grey89', linewidth = 0.35),
        axis.text = element_text(size = 13, colour = 'grey20', face = 'bold'),
        axis.title.x = element_text(size = 14, face = 'bold', colour = 'grey20',
                                    margin = margin(t = 8)),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 18, face = 'bold', hjust = 0,
                                  margin = margin(b = 12)),
        plot.subtitle = element_text(size = 12, colour = 'grey30'),
        legend.position = 'none',
        plot.caption = element_text(size = 9, colour = 'grey45'),
        plot.margin = margin(12, 14, 12, 12))

save_figure <- function(name, figure, width, height) {
  ggsave(file.path(out_dir, paste0(name, '.pdf')), figure,
         width = width, height = height, units = 'in', device = grDevices::pdf,
         limitsize = FALSE)
  ggsave(file.path(out_dir, paste0(name, '.png')), figure,
         width = width, height = height, units = 'in', dpi = 600,
         bg = 'white', limitsize = FALSE)
  message('Saved: ', name)
}

# Shared plotting utility: horizontal interval and centered p-values.
# y_sep: higher separation protects vertically adjacent p-value annotations.
forest_plot <- function(dat, y_breaks, y_labels, x_limits, x_breaks,
                        x_title = 'OR per SD increase in PGS (95% CI)',
                        title = '', subtitle = NULL, log_x = TRUE,
                        label_offset = 0.28, p_size = 4.0,
                        ref_value = 1, ref_width = 0.9,
                        show_p = TRUE) {
  p <- ggplot(dat, aes(x = estimate, y = y)) +
    geom_vline(xintercept = ref_value, linetype = 'dashed',
               colour = 'grey48', linewidth = ref_width) +
    geom_segment(aes(x = lo, xend = hi, yend = y, colour = colour,
                     linetype = interval_type), linewidth = 1.05) +
    geom_point(aes(colour = colour, shape = marker), size = 3.9) +
    scale_colour_identity() +
    scale_shape_identity() +
    scale_linetype_identity() +
    scale_y_continuous(breaks = y_breaks, labels = y_labels,
                       expand = expansion(add = c(0.6, 1.3))) +
    labs(x = x_title, title = title, subtitle = subtitle,
         caption = 'SYNTHETIC DATA — illustration only; not manuscript results') + base_theme
  if (log_x) {
    p <- p + scale_x_log10(limits = x_limits, breaks = x_breaks,
                           labels = as.character(x_breaks),
                           expand = expansion(mult = c(0.02, 0.02)))
  } else {
    p <- p + scale_x_continuous(limits = x_limits, breaks = x_breaks,
                                 expand = expansion(mult = c(0.02, 0.02)))
  }
  if (show_p) {
    p <- p + geom_text(aes(x = estimate, y = y + label_offset, label = p_text),
                       inherit.aes = TRUE, size = p_size, colour = 'black',
                       hjust = 0.5, vjust = 0, na.rm = TRUE,
                       show.legend = FALSE)
  }
  p
}

# ── FIGURE 1: Epilepsy-PGS case-control associations and pTDT ────────────
f1a <- read_example('figure1a_example.csv') %>% mutate(
  pgs = factor(pgs, levels = c('PGS EPI','PGS GGE','PGS NAFE')),
  group = factor(group, levels = c('NDD-only','NDD-Epi solved','NDD-Epi unsolved')),
  y = (as.integer(pgs)-1)*5 + as.integer(group)-1,
  colour = c(COL_GRAY,COL_BLUE,COL_AMBER)[as.integer(group)],
  marker = 16L, interval_type = 'solid', p_text = p_label(p, 9)
)

p1a <- forest_plot(f1a, y_breaks = c(1,6,11),
                   y_labels = c('PGS EPI','PGS GGE','PGS NAFE'),
                   x_limits = c(0.78,1.56),
                   x_breaks = c(0.8,0.9,1,1.1,1.2,1.3,1.4,1.5),
                   title = 'A', label_offset = 0.30)

# Figure 1B reads wholly synthetic illustrative pTDT values.
f1b <- read_example('figure1b_example.csv') %>% mutate(
 lo = estimate - 1.96*se, hi = estimate + 1.96*se,
 y = c(1,3.4,5.8),
 colour = c(COL_GRAY,COL_AMBER,COL_BLUE), marker = 16L,
 interval_type = 'solid', p_text = p_label(p,3)
)

p1b <- forest_plot(f1b, y_breaks = f1b$y,
                   y_labels = paste0(f1b$group, '\n(n=',f1b$n,')'),
                   x_limits = c(-0.050,0.057), x_breaks = c(-0.04,-0.02,0,0.02,0.04),
                   x_title = 'PGS EPI transmission disequilibrium',
                   title = 'B', log_x = FALSE, ref_value = 0,
                   label_offset = 0.30)

# The shared legend occupies the extra space in Panel B rather than a bottom row.
legend1 <- ggplot(data.frame(
 group = factor(c('NDD-only','NDD-Epi solved','NDD-Epi unsolved'),
                levels = c('NDD-only','NDD-Epi solved','NDD-Epi unsolved')),
 x = 1, y = seq_len(3)), aes(x,y)) +
 geom_point(aes(colour = group), size = 4) +
 scale_colour_manual(values = c('NDD-only'=COL_GRAY,
                                'NDD-Epi solved'=COL_BLUE,
                                'NDD-Epi unsolved'=COL_AMBER)) +
 guides(colour = guide_legend(title = NULL, ncol = 1,
                              override.aes = list(size = 4))) +
 theme_void() + theme(legend.position = 'left',
                      legend.text = element_text(size = 12, face='bold'))

# Place legend in a reserved region above Panel B's plot, leaving data intact.
# Using patchwork maintains clean alignment and avoids covering confidence intervals.
p1right <- legend1 + p1b + plot_layout(heights = c(0.28,0.72))
fig1 <- p1a | p1right
save_figure('fig1_combined', fig1, width = 15.8, height = 7.1)

# ── FIGURE 2: GEL and DDD cognitive-PGS replication ──────────────────────
f2gel <- read_example('figure2_gel_example.csv')
f2ddd <- read_example('figure2_ddd_example.csv')
prep_f2 <- function(df) {
 df %>% mutate(
  pgs = factor(pgs, levels = c('PGS CP','PGS EA','PGS IQ')),
  group = factor(group, levels = c('NDD-Epi','NDD-only')),
  y = (as.integer(pgs)-1)*4.5 + as.integer(group)-1,
  colour = ifelse(group=='NDD-only',COL_GRAY,COL_BLUE),
  marker = 16L, interval_type='solid', p_text=p_label(p,6))
}
p2a <- forest_plot(prep_f2(f2gel), c(0.5,5,9.5),
                   c('PGS CP','PGS EA','PGS IQ'), c(0.77,1.22),
                   c(0.8,0.9,1,1.1,1.2), title='A  GEL cohort',
                   subtitle='Synthetic GEL-like example')
p2b <- forest_plot(prep_f2(f2ddd), c(0.5,5,9.5),
                   c('PGS CP','PGS EA','PGS IQ'), c(0.77,1.22),
                   c(0.8,0.9,1,1.1,1.2), title='B  DDD replication',
                   subtitle='Synthetic DDD-like example')
fig2legend <- ggplot(data.frame(x=c(1,2),y=c(1,1),
                                  label=c('NDD-only','NDD-Epi'),
                                  colour=c(COL_GRAY,COL_BLUE)),
                     aes(x=x,y=y)) +
  geom_point(aes(colour=colour),size=4) +
  geom_text(aes(label=label),nudge_x=0.08,hjust=0,size=4.3,fontface='bold') +
  scale_colour_identity() +
  coord_cartesian(xlim=c(0.8,2.55),ylim=c(0.8,1.2),clip='off') +
  theme_void()
fig2 <- (p2a | p2b) / fig2legend + plot_layout(heights=c(1,0.09))
save_figure('fig2_combined', fig2, width=15.8, height=7.1)

# ── FIGURE 3: Bivariate polygenic positioning relative to controls ───────
# Group means (+/- 95% CI) of epilepsy PGS (x) and cognitive PGS (y), in SD
# units relative to the cancer control group, which defines the origin (0,0).
# Scores were residualised on age, sex and the first ten ancestry principal
# components before standardisation; see analysis/ for the derivation.

COL_GREEN_LIGHT <- '#7FD4BC'

ARM_COLOURS <- c('NDD-only unsolved' = COL_GREEN,
                 'NDD-only solved'   = COL_GREEN_LIGHT,
                 'NDD-Epi unsolved'  = COL_AMBER,
                 'NDD-Epi solved'    = COL_BLUE)

f3 <- read_example('figure3_example.csv') %>% mutate(
  key   = paste(group, status),
  key   = factor(key, levels = names(ARM_COLOURS)),
  label = sprintf('%s %s (n=%s)', group, status,
                  format(n, big.mark = ',', trim = TRUE)),
  label = factor(label, levels = label[order(key)]),
  x_lo  = x - 1.96 * x_se, x_hi = x + 1.96 * x_se,
  y_lo  = y - 1.96 * y_se, y_hi = y + 1.96 * y_se
)

f3_palette <- setNames(ARM_COLOURS[as.character(f3$key)],
                       as.character(f3$label))

# Both axes carry the same unit (control SD), so the aspect ratio must be 1:1
f3_pad    <- 0.04
f3_limits <- range(c(f3$x_lo, f3$x_hi, f3$y_lo, f3$y_hi)) + c(-f3_pad, f3_pad)

c_theme <- base_theme +
  theme(panel.grid.major.y = element_line(colour = 'grey89', linewidth = 0.35),
        axis.title.y = element_text(size = 14, face = 'bold', colour = 'grey20',
                                    margin = margin(r = 8)),
        legend.position = c(0.98, 0.02),
        legend.justification = c(1, 0),
        legend.background = element_rect(fill = 'white', colour = NA),
        legend.key.size = unit(10, 'pt'),
        legend.text = element_text(size = 11))

p3 <- ggplot(f3, aes(x, y, colour = label)) +
  geom_hline(yintercept = 0, linetype = 'dotted',
             colour = 'grey55', linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = 'dotted',
             colour = 'grey55', linewidth = 0.5) +
  geom_errorbarh(aes(xmin = x_lo, xmax = x_hi), height = 0, linewidth = 1.05) +
  geom_errorbar(aes(ymin = y_lo, ymax = y_hi), width = 0, linewidth = 1.05) +
  geom_point(aes(fill = label), shape = 21, size = 4.2,
             colour = 'black', stroke = 0.8) +
  scale_colour_manual(values = f3_palette, name = NULL) +
  scale_fill_manual(values = f3_palette, guide = 'none') +
  coord_fixed(ratio = 1, xlim = f3_limits, ylim = f3_limits) +
  labs(x = 'PGS EPI (SD relative to controls)',
       y = 'PGS CP (SD relative to controls)') +
  c_theme

save_figure('fig3', p3, width = 7.6, height = 7.6)

# ── FIGURE 4: Monogenic diagnosis and rare-variant burden ───────────────
f4a <- read_example('figure4a_example.csv') %>% mutate(
 group = factor(group, levels=c('Solved','HC pLoF','LC pLoF / D-miss','Unsolved')),
 y = (as.integer(group)-1)*4.6 + ifelse(pgs=='PGS EPI',1,0),
 colour = ifelse(pgs=='PGS EPI',COL_AMBER,COL_GREEN),
 marker = 16L, interval_type='solid', p_text=p_label(p,8)
)
group4a <- f4a %>% group_by(group) %>% summarise(y=mean(y),n=first(n),.groups='drop')
p4a <- forest_plot(f4a, group4a$y,
                   paste0(group4a$group,'\n(n=',group4a$n,')'),
                   c(0.61,1.69), c(0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6),
                   title='A  Stratified by monogenic cause', label_offset=0.30)

f4b <- read_example('figure4b_example.csv') %>% mutate(
 pgs = factor(pgs, levels=c('PGS EPI','PGS GGE','PGS NAFE','PGS IQ','PGS EA','PGS CP')),
 y=(as.integer(pgs)-1)*4.1 + ifelse(exposure=='No exposure',1.1,0),
 colour=ifelse(pgs %in% c('PGS EPI','PGS GGE','PGS NAFE'),COL_AMBER,COL_GREEN),
 marker=ifelse(exposure=='No exposure',15L,17L),
 interval_type=ifelse(exposure=='No exposure','solid','dashed'),
 p_text=p_label(p,12)
)
p4b <- forest_plot(f4b,
                   y_breaks=seq(0.55,by=4.1,length.out=6),
                   y_labels=levels(f4b$pgs),
                   x_limits=c(0.73,1.73),
                   x_breaks=c(0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6),
                   title='B  Stratified by rare-variant burden',
                   subtitle='Synthetic exposure groups',
                   label_offset=0.19, p_size=3.9)

# Legend: color = PGS family; line and shape = exposure status in Panel B.
# Render once below both panels, matching the final figure's visual grammar.
leg_df <- data.frame(
 x=1:4, y=1,
 key=factor(c('Epilepsy PGS','Cognitive PGS','No exposure','+1 exposure'),
            levels=c('Epilepsy PGS','Cognitive PGS','No exposure','+1 exposure'))
)
leg_plot <- ggplot(leg_df,aes(x=x,y=y))+
 geom_segment(aes(x=x-0.22,xend=x+0.22,yend=y,
                  linetype=key,colour=key),linewidth=1.1)+
 geom_point(aes(shape=key,colour=key),size=4)+
 scale_colour_manual(values=c('Epilepsy PGS'=COL_AMBER,'Cognitive PGS'=COL_GREEN,
                              'No exposure'='grey30','+1 exposure'='grey30'))+
 scale_shape_manual(values=c('Epilepsy PGS'=16,'Cognitive PGS'=16,
                             'No exposure'=15,'+1 exposure'=17))+
 scale_linetype_manual(values=c('Epilepsy PGS'='blank','Cognitive PGS'='blank',
                                'No exposure'='solid','+1 exposure'='dashed'))+
 scale_x_continuous(breaks=1:4,labels=levels(leg_df$key))+
 coord_cartesian(ylim=c(0.88,1.12),clip='off')+
 theme_void()+theme(axis.text.x=element_text(size=13,face='bold',colour='grey25'),
                    plot.margin=margin(0,30,0,30))
fig4 <- (p4a | p4b) / leg_plot + plot_layout(heights=c(1,0.12))
save_figure('fig4_combined',fig4,width=16.0,height=8.4)
message('Figures 1, 2 , 3 and 4 exported.')