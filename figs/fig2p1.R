source("figs/fig2p0.R")

inp.l <- local({
  base.id <- "puck_210218_07"
  ad <- reticulate::import("anndata", convert = FALSE)

  orig <- ad$read_h5ad(sprintf("data/%s.h5ad", base.id))
  bg <- ad$read_h5ad(sprintf("data/sP_bench/inp/bg/%s.h5ad", base.id))
  fg.r <- function(fg.id) ad$read_h5ad(sprintf("data/sP_bench/inp/fg/%s_%s.h5ad", base.id, fg.id))

  ids.oi <- conf.DT[frac == 1][n.genes == 89][abs(x.cent - 3332.3) < 1e-4][abs(y.cent - 3291.6) < 1e-4][, fg.id]

  c.nm <- reticulate::py_to_r(orig$obs_names$to_list())
  g.nm <- reticulate::py_to_r(orig$var_names$to_list())
  is.fg <- local({
    x <- rep_len(TRUE, length(c.nm))
    x[match(reticulate::py_to_r(bg$obs_names$to_list()), c.nm)] <- FALSE
    x
  })

  within.circle <- local({
    smp.oi <- fg.r(106)

    x <- rep_len(FALSE, length(c.nm))
    x[match(reticulate::py_to_r(smp.oi[smp.oi$uns["pert_conf"]["orig_cells"], ]$obs_names$to_list()), c.nm)] <- TRUE
    x
  })

  list(
    X = `dimnames<-`(reticulate::py_to_r(orig$X$T), list(g.nm, c.nm)),
    meta = cbind(
      as.data.table(reticulate::py_to_r(orig$obs), "cell"),
      setnames(as.data.table(reticulate::py_to_r(orig$obsm["spatial"])), new = c("x", "y")),
      is.fg = is.fg,
      within.circle = within.circle
    ),
    deltas = rbindlist(purrr::map(ids.oi, function(fg.id) {
      fg <- fg.r(fg.id)
      c.idx <- as.vector(reticulate::py_to_r(fg$uns$pert_conf$orig_cells) + 1)
      g.idx <- as.vector(reticulate::py_to_r(fg$uns$pert_conf$orig_genes) + 1)
      out <- as.matrix(reticulate::py_to_r(fg$uns$pert_conf["delta"]$T)[g.idx, c.idx])
      dimnames(out) <- list(reticulate::py_to_r(fg$var_names$to_list())[g.idx], reticulate::py_to_r(fg$obs_names$to_list())[c.idx])

      # sanity check
      stopifnot(all(
        (reticulate::py_to_r(
          fg$raw$to_adata()$T[rownames(out), colnames(out)]$X
        ) +
          out) ==
          reticulate::py_to_r(
            fg$T[rownames(out), colnames(out)]$X
          )
      ))

      melt(
        as.data.table(out, "gene"),
        id.vars = "gene",
        variable.name = "cell"
      )[, let(λ = reticulate::py_to_r(fg$uns$pert_conf["λ"]))][]
    })),
    means = rbindlist(purrr::map(ids.oi, function(fg.id) {
      fg <- fg.r(fg.id)
      fg <- fg[fg$uns$pert_conf$orig_cells, ]
      data.table(
        expr.base = sparseMatrixStats::rowMeans2(reticulate::py_to_r(fg$raw$X$T)),
        expr.noise = sparseMatrixStats::rowMeans2(reticulate::py_to_r(fg$X$T)),
        λ = reticulate::py_to_r(fg$uns$pert_conf["λ"]),
        gene = reticulate::py_to_r(fg$var_names$to_list())
      )
    }))
  )
})

des.pl <- list(
  ggplot(inp.l[["meta"]]) +
    aes(x = x, y = y, colour = annotation) +
    geom_point(size = 0.2) +
    geom_hline(aes(yintercept = mean(max(y[!is.fg]), min(y[is.fg]))), linetype = "dotted", colour = "black", linewidth = 2) +
    scale_colour_brewer(palette = "RdBu", limits = function(e) c(e[e != "NA"], "NA")) +
    coord_fixed() +
    labs(colour = NULL, alpha = NULL) +
    guides(colour = guide_legend(override.aes = list(shape = 16, size = 2)), alpha = guide_legend(override.aes = list(shape = 16, size = 2))) +
    theme(axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()),
  ggplot(inp.l[["meta"]]) +
    geom_point(
      aes(x = collapse::fmedian(x, is.fg, TRA = 3), y = collapse::fmedian(y, is.fg, TRA = 3), colour = annotation),
      size = 0.2
    ) +
    ggforce::geom_circle(
      aes(x0 = x.cent, y0 = y.cent, r = rad),
      conf.DT[fg.id == 106][, .(rad, x.cent = x.cent - inp.l[["meta"]][(is.fg), median(x)], y.cent = y.cent - inp.l[["meta"]][(is.fg), median(y)], is.fg = TRUE)],
      colour = "red",
      linetype = "dashed",
      linewidth = 1.5
    ) +
    facet_wrap(vars(forcats::fct_rev(fifelse(is.fg, "foreground", "background"))), ncol = 1) +
    scale_colour_brewer(palette = "RdBu", limits = function(e) c(e[e != "NA"], "NA")) +
    coord_fixed() +
    labs(colour = NULL, alpha = NULL, x = "x", y = "y") +
    guides(colour = guide_legend(override.aes = list(shape = 16, size = 2)), alpha = guide_legend(override.aes = list(shape = 16, size = 2))) +
    theme(axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()),
  ggplot() +
    purrr::map(conf.DT[, sort(unique(λ))], function(λ) {
      geom_function(aes(colour = as.character(λ)), fun = function(x) (λ^x) * exp(-λ) / factorial(x), n = 1000, linewidth = 2, alpha = 0.9)
    }) +
    scale_colour_viridis_d(option = "magma", direction = -1, end = 0.8) +
    lims(x = c(0, inp.l[["deltas"]][, max(value)])) +
    labs(x = "synthetic perturbation noise value", y = "poisson density", colour = "λ") +
    theme(axis.title.x = element_text(margin = margin(t = 0.6, unit = "lines"))),
  inp.l[["means"]][λ == min(λ)][sample(.N)] |>
    ggplot() +
    aes(x = expr.base, y = expr.noise) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_point(colour = "#133767", alpha = 0.8, size = 2) +
    coord_fixed() +
    ggrepel::geom_label_repel(aes(label = gene), function(e) e[expr.noise > expr.base], max.overlaps = 40, size = 5, show.legend = FALSE, seed = 8) +
    guides(colour = guide_legend(override.aes = list(size = 3))) +
    labs(x = "mean counts in circle (control)", y = "mean counts in circle (perturbed)"),
  local({
    base.DT <- inp.l[["meta"]]
    base.expr <- inp.l[["X"]][, base.DT[, cell]]
    noise.expr <- split(inp.l[["deltas"]][value > 0], by = "λ") |>
      purrr::map(function(noise.expr) {
        with(
          noise.expr,
          Matrix::sparseMatrix(
            x = value,
            i = match(gene, rownames(base.expr)),
            j = match(as.character(cell), colnames(base.expr)),
            dims = dim(base.expr),
            dimnames = dimnames(base.expr)
          )
        ) +
          base.expr
      })

    lef <- function(x, gene) {
      SummarizedExperiment::assay(
        scrapper::normalizeRnaCounts.se(SummarizedExperiment::SummarizedExperiment(list(counts = x))),
        "logcounts"
      )[gene, ]
    }
    rbindlist(
      append(
        list(original = cbind(base.DT, Camk1d = lef(base.expr, "Camk1d"))),
        setNames(purrr::map(noise.expr, function(e) cbind(base.DT, Camk1d = lef(e, "Camk1d"))), sprintf("λ=%s", names(noise.expr)))
      ),
      idcol = "grp"
    )
  }) |>
    ggplot() +
    aes(x = collapse::fmedian(x, is.fg, TRA = 3), y = collapse::fmedian(y, is.fg, TRA = 3), colour = Camk1d) +
    geom_point(size = 0.01) +
    scale_colour_viridis_c(option = "magma", direction = -1) +
    facet_grid(vars(forcats::fct_relevel(fifelse(is.fg, "foreground", "background"), "background")), vars(grp), space = "free") +
    coord_fixed() +
    theme(axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
) |>
  purrr::map(function(pl) {
    pl +
      theme(
        panel.background = element_rect(fill = "transparent", colour = "black"),
        strip.background = element_rect(fill = "transparent", colour = "black"),
        legend.key = element_blank()
      )
  })

if (!interactive()) {
  list(
    des.pl,
    seq_along(des.pl),
    list(
      c(w = 2000, h = 2000),
      c(w = 2200, h = 3000),
      c(w = 1800, h = 1500),
      c(w = 1500, h = 1500),
      c(w = 4800, h = 1800)
    )
  ) |>
    purrr::pwalk(function(e, i, s) gps(sprintf("figs/fig2_%s", letters[i]), e, units = "px", w = s[["w"]], h = s[["h"]]))
}
