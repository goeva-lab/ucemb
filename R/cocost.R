argv <- commandArgs(trailingOnly = TRUE)

bg.ad <- anndataR::read_h5ad(argv[[1]])
fg.ad <- anndataR::read_h5ad(argv[[2]])
k <- as.integer(argv[[3]])
f_name <- file.path(argv[[4]], "emb.npz")
c_param <- as.numeric(argv[[5]])

so <- SeuratObject:::merge.Seurat(
  SeuratObject::AddMetaData(
    anndataR:::as_Seurat(bg.ad, x_mapping = "counts"),
    "bg",
    "cond"
  ),
  SeuratObject::AddMetaData(
    anndataR:::as_Seurat(fg.ad, x_mapping = "counts"),
    "fg",
    "cond"
  )
) |>
  Seurat::NormalizeData() |>
  SeuratObject::JoinLayers() |>
  Seurat::FindVariableFeatures() |>
  Seurat::ScaleData(split.by = "cond") |>
  Seurat::SplitObject(split.by = "cond")


bg <- t(as.matrix(SeuratObject::LayerData(so[["bg"]], layer = "scale.data")))
fg <- t(as.matrix(SeuratObject::LayerData(so[["fg"]], layer = "scale.data")))

kern <- kernlab::laplacedot(sigma = 0.5)

ccst.out <- CoCoST::CoCoST(
  fg,
  kernlab::kernelMatrix(kern, fg)@.Data,
  bg,
  kernlab::kernelMatrix(kern, bg)@.Data,
  para = c_param,
  Dim = k
)

z <- rbind(ccst.out[["bgComponents"]], ccst.out[["fgComponents"]])

library(reticulate)
nps <- import("numpy", convert = FALSE)$save
pop <- import_builtins()$open
with(pop(f_name, "wb") %as% f, {
  nps(f, z)
})
