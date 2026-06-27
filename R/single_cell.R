#' scGPT
#' 
#' This function runs scGPT.
#' 
#' @param h5ad_file file that points to the location of a h5ad file with the
#'   count data
#' @param model_dir directory with the model. This should be a directory that
#'   the scGPT model, consisting of a `.pt` file with the model, as well as the
#'   `vocab.json` and `args.json` file.
#' @param gene_col The gene column name
#' @param batch_size The batch size passed to `scg.tasks.embed_data` (default 64)
#' @param device torch device ("cpu", "cuda", "mps" etc.)
#' 
#' @details
#' This function runs the `scg.tasks.embed_data` function
#'
#' @export
Run_scGPT <- function(h5ad_file, model_dir, gene_col, batch_size = 64L, device = "cpu") {
  proc <- basilisk::basiliskStart(.scgpt)
  on.exit(basilisk::basiliskStop(proc))
  basilisk::basiliskRun(proc, function(h5ad_file, model_dir, gene_col, batch_size, device) {
    
    # macOS fix: os.sched_getaffinity is Linux-only; patch it if missing
    os <- reticulate::import("os")
    if (!reticulate::py_has_attr(os, "sched_getaffinity")) {
      reticulate::py_run_string(
        "import os; os.sched_getaffinity = lambda x: set(range(os.cpu_count()))"
      )
    }
    
    # macOS fix: default multiprocessing start method is "spawn" on macOS, which
    # cannot pickle locally-defined classes used by scGPT's DataLoader; force "fork"
    mp <- reticulate::import("multiprocessing")
    start_method <- mp$get_start_method(allow_none = TRUE)
    if (is.null(start_method) || start_method != "fork") {
      mp$set_start_method("fork", force = TRUE)
    }
    sg <- reticulate::import("scgpt")
    sc <- reticulate::import("scanpy")
    
    # read data into an anndata
    adata <- sc$read_h5ad(h5ad_file)
    
    # run scgpt
    ref_embed_adata <- sg$tasks$embed_data(
      adata,
      model_dir,
      gene_col = gene_col,
      batch_size = as.integer(batch_size),
      device = device
    )
    
    # return embedding
    return(ref_embed_adata$obsm[["X_scGPT"]])
  }, h5ad_file = h5ad_file, model_dir = model_dir, gene_col = gene_col, batch_size = batch_size, device = device)
}