if (!exists("PRELUDE.HAS.RUN")) {
  suppressPackageStartupMessages({
    library(ggplot2)
    library(patchwork)
    library(stringr)
    library(data.table)
  })

  conf.DT <- local({
    ad <- reticulate::import("anndata", convert = FALSE)
    o <- purrr::map(
      list.files("data/sP_bench/inp/fg", pattern = ".*\\.h5ad", full.names = TRUE),
      function(p) {
        o <- ad$read_h5ad(p, backed = "r")
        on.exit(o$file$close(), add = TRUE)

        col.oi <- c("frac", "rad", "x_cent", "y_cent", "λ")
        append(
          purrr::map(
            setNames(col.oi, str_replace(col.oi, fixed("_"), ".")),
            function(k, o) reticulate::py_to_r(o[k]),
            o = o$uns["pert_conf"]
          ),
          list(
            n.genes = o$var$perturbed$sum() |> reticulate::py_to_r() |> as.integer(),
            fg.id = p |> str_split_i(fixed("."), 1) |> str_split_i(fixed("_"), -1) |> as.integer()
          )
        )
      }
    )
    as.data.table(purrr::map(transpose(append(o, list(nm = names(o[[1]]))), make.names = "nm"), unlist))[order(fg.id)]
  })

  meth.ord <- c("pca", "cpca", "ccpca", "gcpca", "pcpca", "cocost", "concord", "cvi", "sp")

  gps <- function(filename, plot, ..., formats = c("pdf" = cairo_pdf, "png" = png)) {
    purrr::iwalk(formats, function(dev, fmt) ggsave(sprintf("%s.%s", filename, fmt), plot, device = dev, ...))
  }

  update_theme(
    text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5),
    legend.key.spacing.y = unit(0.2, "lines"),
    panel.grid = element_blank()
  )

  PRELUDE.HAS.RUN <- TRUE
}
