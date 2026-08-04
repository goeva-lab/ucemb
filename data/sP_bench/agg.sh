#!/bin/bash

set -euxo pipefail

rg --no-filename -or '$1' -e '^::::\[\[ ZSH TIME \]\]:::: (.+) ::::\[\[ ZSH TIME \]\]::::' 'out/'*'/'*'/run.log' > './run.log'

Rscript /dev/stdin <<'EOF'
library(stringr)
library(data.table)
lines <- str_split(readLines("run.log"), fixed(" | "))
fwrite(
  as.data.table(lapply(
    transpose(append(lines, list(nm = str_to_lower(str_split_i(lines[[1]], fixed("="), 1)))), make.names = "nm"),
    str_split_i,
    pattern = fixed("="),
    i = 2
  ))[, local({
    x <- str_split(str_split_i(str_split_i(fg_h5ad_path, fixed("/"), -1), fixed("."), 1), fixed("_"))
    x <- matrix(unlist(x), byrow = 1, ncol = length(x[[1]]))

    .(
      method = as.factor(method_name),
      base.data = as.factor(apply(x[, -ncol(x)], 1, str_c, collapse = "_")),
      fg.id = as.integer(x[, ncol(x)]),
      emb.dim = as.integer(emb_dim),
      c.param = as.integer(c_param),
      time.sec = as.numeric(str_remove(time, regex("s$"))),
      cpu.pct = as.numeric(str_remove(pct, regex("%$"))),
      mem.mb = as.numeric(mem)
    )
  })],
  "perf.tsv"
)

invisible(lapply(c("per_j_scores", "full_scores"), function(e) {
  fwrite(
    rbindlist(
      lapply(
        list.files("out", pattern = sprintf("%s.tsv", e), recursive = TRUE, full.names = TRUE), 
        function(e) suppressWarnings(fread(e, sep = "\t")[, -c('orig_cells', 'orig_genes')])
      )
    ),
    sprintf("%s.tsv", e),
    sep = "\t"
  )
}))
EOF
