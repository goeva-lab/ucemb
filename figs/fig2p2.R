source("figs/fig2p0.R")

scores.DT <- fread("data/sP_bench/full_scores.tsv")[, .SD[.N == nrow(conf.DT)], by = .(method, stat)][λ > 0.5]
perf.DT <- fread("data/sP_bench/perf.tsv")[, .SD[.N == nrow(conf.DT)], by = .(method)]
scores.pl <- ggplot(scores.DT[startsWith(stat, "label recovery")]) +
  aes(y = value, x = forcats::fct_relevel(method, meth.ord)) +
  geom_violin(scale = "width", show.legend = FALSE) +
  geom_point(aes(fill = as.factor(frac * 100)), position = position_dodge(0.65), alpha = 0.5, size = 3, shape = "circle filled") +
  facet_grid(vars(`number of genes` = n_genes), vars(λ), labeller = label_both, scales = "free") +
  scale_fill_viridis_d(option = "inferno", limits = rev) +
  labs(y = "score", x = NULL, fill = "% affected cells", title = "label recovery jaccard index") +
  guides(fill = guide_legend(override.aes = list(alpha = 1))) +
  lims(y = c(0, 1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

perf.pl <- wrap_plots(
  ggplot(perf.DT) +
    aes(y = time.sec / 3600, x = forcats::fct_relevel(method, meth.ord), fill = cpu.pct) +
    scale_fill_viridis_c(option = "inferno", limits = c(0, 110)) +
    labs(y = "wall time (hours, log-scale)", fill = "% CPU use"),
  ggplot(perf.DT) +
    aes(y = mem.mb / 1024, x = forcats::fct_relevel(method, meth.ord), fill = time.sec / 3600) +
    scale_fill_viridis_c(option = "cividis", limits = c(perf.DT[, min(time.sec) / 3600], perf.DT[, (max(time.sec) / 3600) * 4]), trans = scales::transform_log10(), labels = function(e) sprintf("%.02f", e)) +
    labs(y = "peak memory use (GB, log-scale)", fill = "wall time\n(hours)"),
  ncol = 1,
  axis_titles = "collect"
) &
  geom_violin(aes(fill = NULL), scale = "width") &
  geom_point(position = position_jitter(0.3), shape = "circle filled", alpha = 0.4, stroke = 0) &
  scale_y_log10() &
  labs(x = "method") &
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

if (!interactive()) {
  gps("figs/fig2_label_recovery", scores.pl, w = 5000, h = 3000, units = "px")
  gps("figs/fig2_perf", perf.pl, w = 2000, h = 3500, units = "px")
}
