#' Run_novae
#' 
#' This function runs novae.
#' 
#' @param adata_path Path to the input anndata file.
#' @param accelerator Accelerator to use, either "cpu" or "cuda". Default is "cpu".
#' 
#' @export
Run_novae <- function(adata_path = NULL, accelerator = "cpu") {
  proc <- basilisk::basiliskStart(.novae)
  on.exit(basilisk::basiliskStop(proc))
  basilisk::basiliskRun(proc, function(adata_path, accelerator) {
    
    # libraries
    os    <- reticulate::import("os")
    os$environ[["HF_HOME"]]      <- ".cache/huggingface"
    os$environ[["MPLCONFIGDIR"]] <- ".cache/matplotlib"
    novae <- reticulate::import("novae")
    ad <- reticulate::import("anndata")
    np <- reticulate::import("numpy")
    sp <- reticulate::import("scipy")
    
    # read anndata
    adata <- ad$read_h5ad(adata_path)
    
    # sparse matrix
    adata$X <- sp$sparse$csr_matrix(adata$layers$counts)

    # spatial neighbors
    novae$spatial_neighbors(adata)
    
    # novae
    model <- novae$Novae$from_pretrained("MICS-Lab/novae-human-0")
    n_valid_cells <-  as.integer(adata$n_obs())
    model$swav_head$num_prototypes <- min(
      model$swav_head$num_prototypes,
      n_valid_cells %/% 2L
    )

    model$compute_representations(adata, zero_shot = TRUE, accelerator = accelerator, num_workers = 4L)
    model$assign_domains(adata)
    reticulate::py_to_r(adata)
  }, adata_path = adata_path, accelerator = accelerator)
}

#' Run_nimbus
#' 
#' This functions runs Nimbus.
#' 
#' @export
Run_nimbus <- function() {
  proc <- basilisk::basiliskStart(.nimbus)
  on.exit(basilisk::basiliskStop(proc))
  basilisk::basiliskRun(proc, function() {
    sg <- reticulate::import("nimbus_inference")
    message("nimbus_inference was loaded!")
    return(TRUE)
  })
}


#' Run_nicheformer
#' 
#' This functions runs NicheFormer.
#' 
#' @param adata_path Path to the input anndata file.
#' @param technology Technology used for the spatial transcriptomics data. One of "cosmx", "visium", "xenium", "merfish", "iss", or "dissociated".
#' 
#' @export
Run_nicheformer <- function(adata_path = NULL,
                            technology = c("cosmx", "visium", "xenium", "merfish", "iss", "dissociated")) {
  technology <- match.arg(technology)
  
  technology_mean_path <- system.file(
    "extdata", "model_means",
    paste0(technology, "_mean_script.npy"),
    package = "fomo"
  )
  if (!nzchar(technology_mean_path)) {
    stop("No model mean file found for technology: ", technology)
  }
  
  proc <- basilisk::basiliskStart(.nicheformer)
  on.exit(basilisk::basiliskStop(proc))
  basilisk::basiliskRun(proc, function(adata_path, technology_mean_path) {
    
    os <- reticulate::import("os")
    os$environ[["HF_HOME"]] <- ".cache/huggingface"
    os$environ[["MPLCONFIGDIR"]] <- ".cache/matplotlib" 
    
    transformers <- reticulate::import("transformers")
    AutoModelForMaskedLM <- transformers$AutoModelForMaskedLM
    AutoTokenizer <- transformers$AutoTokenizer
    ad <- reticulate::import("anndata")
    np <- reticulate::import("numpy")
    torch <- reticulate::import("torch")
    sp <- reticulate::import("scipy")
    builtins <- reticulate::import_builtins()
    
    # Load model and tokenizer
    model <- AutoModelForMaskedLM$from_pretrained("aletlvl/Nicheformer", trust_remote_code=TRUE)
    tokenizer <- AutoTokenizer$from_pretrained("aletlvl/Nicheformer", trust_remote_code=TRUE)
    
    # Set technology mean for HF tokenizer
    technology_mean <- np$load(technology_mean_path)
    tokenizer$`_load_technology_mean`(technology_mean)
    
    adata <- ad$read_h5ad(adata_path)
    adata$X <- sp$sparse$csr_matrix(adata$layers$counts)
    
    #########################
    
    # Device selection
    device <- if (torch$cuda$is_available()) torch$device("cuda") else torch$device("cpu")
    #device <- torch$device("cpu")
    cat(sprintf("Using device: %s\n", device$type))
    
    torch$cuda$empty_cache()
    
    # half() and to() return self, so chaining works
    if (torch$cuda$is_available()) {
      model <- model$half()$to(device)
    } else {
      model <- model$to(device)
    }
    
    model$eval()
    
    # Tokenize
    inputs <- tokenizer(adata)
    
    # Truncate to max_len. Python's tensor[:, :max_len] becomes narrow(dim, start, length).
    max_len <- 653L                                     # 650 + 3
    input_ids      <- inputs[["input_ids"]]$narrow(1L, 0L, max_len)
    attention_mask <- inputs[["attention_mask"]]$narrow(1L, 0L, max_len)
    
    # Patch positional index. $<- works for Python attribute assignment.
    model$nicheformer$pos <- model$nicheformer$pos$narrow(0L, 0L, max_len)
    
    # Batch processing
    batch_size     <- 64L
    n_samples      <- input_ids$size(0L)                # safer than $shape[[...]]
    all_embeddings <- list()
    
    # reticulate's with() dispatches to Python context managers
    with(torch$no_grad(), {
      for (start in seq(0L, n_samples - 1L, by = batch_size)) {
        cur <- as.integer(min(batch_size, n_samples - start))
        
        batch_ids  <- input_ids$narrow(0L, as.integer(start), cur)$to(device)
        batch_mask <- attention_mask$narrow(0L, as.integer(start), cur)$to(device)
        
        emb <- model$get_embeddings(
          input_ids       = batch_ids,
          attention_mask  = batch_mask,
          layer           = -1L,
          with_context    = TRUE
        )
        all_embeddings[[length(all_embeddings) + 1L]] <- emb$cpu()
        torch$cuda$empty_cache()
      }
    })
    
    embeddings <- torch$cat(all_embeddings, dim = 0L)
    
    shape_str <- paste(unlist(builtins$list(embeddings$shape)), collapse = ", ")
    cat(sprintf("Embeddings shape: [%s]\n", shape_str))
    
    # Return something R-friendly. Cast away from fp16 first, numpy's fp16
    # support is patchy and you usually don't want half-precision downstream.
    return(embeddings$float()$cpu()$numpy())
    
  }, adata_path = adata_path, technology_mean_path = technology_mean_path)
}
