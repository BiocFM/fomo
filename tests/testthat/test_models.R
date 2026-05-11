library(anndataR)
library(SpatialExperiment)
library(fomo)

test_that("test novae", {
  
  # read spe
  spe <- readRDS(
    system.file("extdata", "CosMx1k_MouseBrain1_100tx_100cl.rds", package = "fomo")
  )
  spatialCoords(spe) |> as.data.frame() -> coords_df
  colData(spe)$x_coord <- coords_df[, 1]
  colData(spe)$y_coord <- coords_df[, 2]
  reducedDim(spe, "spatial") <- as.matrix(spatialCoords(spe))
  
  # anndata
  adata <- as_AnnData(spe)
  colnames(adata$obsm$spatial) <- c("x_coord", "y_coord")
  adata$write_h5ad(
    tp <- tempfile(fileext = ".h5ad"), 
    mode = "w"
  )
  
  # run novae
  novae_data <- Run_novae(tp, accelerator = "cpu")
  novae_data$obsm$novae_latent |> head()
})