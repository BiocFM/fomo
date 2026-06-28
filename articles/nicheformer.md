# Running Nicheformer foundational model

Abstract

This package takes an h5ad file with the gene count matrix and spatial
coordinates as input, runs Nicheformer creating the Python environment
automatically, and returns a matrix with the embeddings.

Loading the necessary packages to run `Nicheformer` in R

``` r

library(anndataR)
library(fomo)
library(SpatialExperiment)
```

    ## Loading required package: SingleCellExperiment

    ## Loading required package: SummarizedExperiment

    ## Loading required package: MatrixGenerics

    ## Loading required package: matrixStats

    ## 
    ## Attaching package: 'MatrixGenerics'

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
    ##     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
    ##     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
    ##     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
    ##     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
    ##     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
    ##     colWeightedMeans, colWeightedMedians, colWeightedSds,
    ##     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
    ##     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
    ##     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
    ##     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
    ##     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
    ##     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
    ##     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
    ##     rowWeightedSds, rowWeightedVars

    ## Loading required package: GenomicRanges

    ## Loading required package: stats4

    ## Loading required package: BiocGenerics

    ## Loading required package: generics

    ## 
    ## Attaching package: 'generics'

    ## The following objects are masked from 'package:base':
    ## 
    ##     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
    ##     setequal, union

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
    ##     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
    ##     get, grep, grepl, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, saveRDS, table, tapply, unique,
    ##     unsplit, which.max, which.min

    ## Loading required package: S4Vectors

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following object is masked from 'package:utils':
    ## 
    ##     findMatches

    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname

    ## Loading required package: IRanges

    ## Loading required package: Seqinfo

    ## Loading required package: Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## 
    ## Attaching package: 'Biobase'

    ## The following object is masked from 'package:MatrixGenerics':
    ## 
    ##     rowMedians

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     anyMissing, rowMedians

Preparing the AnnData object

``` r

spe <- readRDS(
  system.file("extdata", "CosMx1k_MouseBrain1_100tx_100cl.rds", package = "fomo")
)

adata <- as_AnnData(spe)

# add coordinates to metadata
spatialCoords(spe) |> as.data.frame() -> coords_df
adata$obs$x_coord <- coords_df[, 1]
adata$obs$y_coord <- coords_df[, 2]

# add coordinates as spatial
coords <- as.matrix(spatialCoords(spe))
rownames(coords) <- adata$obs_names
adata$obsm$spatial <- coords
colnames(adata$obsm$spatial) <- c("x_coord", "y_coord")
```

Here, we are using a mouse cosmx-based example with gene symbols as
variable names. However, `Nicheformer` requires human Ensembl IDs.
Therefore, this chunk retrieves and replaces the orthologues

``` r

library(homologene)
library(org.Hs.eg.db)
```

    ## Loading required package: AnnotationDbi

    ## 

``` r

library(AnnotationDbi)

gene_symbols <- adata$var_names

# mouse (taxid 10090) -> human (taxid 9606) ortholog symbols
orthologs <- homologene(gene_symbols, inTax = 10090, outTax = 9606)
# columns: "10090" (mouse symbol), "9606" (human symbol)

# human symbol -> human Ensembl ID
human_ensembl <- mapIds(
  org.Hs.eg.db,
  keys      = orthologs[["9606"]],
  column    = "ENSEMBL",
  keytype   = "SYMBOL",
  multiVals = "first"
)
```

    ## 'select()' returned 1:many mapping between keys and columns

``` r

orthologs$ensembl_id <- human_ensembl[orthologs[["9606"]]]

# match back to original gene order
new_names <- ifelse(
  gene_symbols %in% orthologs[["10090"]],
  orthologs$ensembl_id[match(gene_symbols, orthologs[["10090"]])],
  gene_symbols
)

# keep only valid Ensembl IDs (start with ENSG) and not duplicated
is_ensembl  <- grepl("^ENSG", new_names)
is_unique   <- !duplicated(new_names)
keep        <- is_ensembl & is_unique

message(sum(keep), " genes kept out of ", length(new_names))
```

    ## 896 genes kept out of 960

``` r

adata <- adata[, keep]$as_InMemoryAnnData()
rownames(adata$var) <- new_names[keep]
```

Generating the `h5ad` file and running `Nicheformer` in R

``` r

# write anndata to tempfile 
adata$write_h5ad(
  tp <- tempfile(fileext = ".h5ad"), 
  mode = "w"
)
```

    ## Warning: Matrix column names cannot be written to a <HDF5AnnData> object, they will be
    ## lost
    ## ℹ To write column names for obsm[['spatial']], store it as <data.frame> instead
    ##   of a double matrix
    ## ℹ NOTE: obs_names and var_names are stored separately

``` r

nicheformer_data <- Run_nicheformer(adata_path = tp,
                                    technology = "cosmx")
```

    ## Installing pyenv ...
    ## Done! pyenv has been installed to '/home/runner/.local/share/r-reticulate/pyenv/bin/pyenv'.
    ## Using Python: /home/runner/.pyenv/versions/3.10.0/bin/python3.10
    ## Creating virtual environment '/home/runner/.cache/R/basilisk/1.24.0/fomo/0.1.0/nicheformer' ...

    ## + /home/runner/.pyenv/versions/3.10.0/bin/python3.10 -m venv /home/runner/.cache/R/basilisk/1.24.0/fomo/0.1.0/nicheformer

    ## Done!
    ## Installing packages: pip, wheel, setuptools

    ## + /home/runner/.cache/R/basilisk/1.24.0/fomo/0.1.0/nicheformer/bin/python -m pip install --upgrade pip wheel setuptools

    ## Installing packages: 'transformers==4.57.6', 'tiktoken==0.9.0', 'sentencepiece==0.2.1', 'git+https://github.com/theislab/nicheformer.git@485cadbc5caa15119adfd54228f8a8af835fcabc'

    ## + /home/runner/.cache/R/basilisk/1.24.0/fomo/0.1.0/nicheformer/bin/python -m pip install --upgrade --no-user 'transformers==4.57.6' 'tiktoken==0.9.0' 'sentencepiece==0.2.1' 'git+https://github.com/theislab/nicheformer.git@485cadbc5caa15119adfd54228f8a8af835fcabc'

    ## Virtual environment '/home/runner/.cache/R/basilisk/1.24.0/fomo/0.1.0/nicheformer' successfully created.
    ## Using provided technology mean array with shape (20310,)
    ## [1] "cpu"
    ## Using device: cpu
    ## Embeddings shape: [100, 512]

``` r

reducedDim(spe, "nicheformer") <- nicheformer_data
```

Downstream analysis using the embeddings to calculate clusters using a
graph-based approach

``` r

library(scrapper)

seed <- 1234
set.seed(seed)

g <- buildSnnGraph(t(nicheformer_data), 
                   num.neighbors = 30, 
                   weight.scheme = "jaccard")
leiden_clusters <- clusterGraph(g, 
                                method = "leiden", 
                                leiden.resolution = 1.0, 
                                leiden.objective = "modularity")
spe$leiden_clusters <- leiden_clusters$membership

message(sprintf("Leiden found %d clusters", 
                length(unique(leiden_clusters))))
```

    ## Leiden found 2 clusters

The resulting embedding of Nicheformer will still be of high dimensions,
we can further reduce the dimensionality with PCA and UMAP.

``` r

# reduce with PCA
result_pca <- runPca(t(nicheformer_data), 
                     number = 30, 
                     scale = TRUE, 
                     seed = seed)
result_pca <- result_pca$components
reducedDim(spe, "nicheformer_PCA") <- t(result_pca)

# reduce with UMAP
result_umap <- runUmap(result_pca, 
                       optimize.seed = seed)
reducedDim(spe, "nicheformer_UMAP") <- result_umap
```

Let us visualize the UMAP and spatial observations

``` r

library(ggspavis)
```

    ## Loading required package: ggplot2

``` r

library(patchwork)
library(scater)
```

    ## Loading required package: scuttle

    ## 
    ## Attaching package: 'scuttle'

    ## The following objects are masked from 'package:scrapper':
    ## 
    ##     aggregateAcrossCells, normalizeCounts

``` r

p_umap <- plotUMAP(spe, dimred = "nicheformer_UMAP")

spe$in_tissue <- 1
spe$x_centroid <- spe$y_centroid <- NULL
p_spatial <- plotCoords(spe, annotate = "leiden_clusters", point_size = 2)

p_umap + p_spatial
```

![](nicheformer_files/figure-html/vis-1.png)

``` r

sessionInfo()
```

    ## R version 4.6.1 (2026-06-24)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 24.04.4 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
    ##  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
    ##  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
    ## [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
    ## 
    ## time zone: UTC
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] scater_1.40.1               scuttle_1.22.0             
    ##  [3] patchwork_1.3.2             ggspavis_1.18.0            
    ##  [5] ggplot2_4.0.3               scrapper_1.6.3             
    ##  [7] org.Hs.eg.db_3.23.1         AnnotationDbi_1.74.0       
    ##  [9] homologene_1.4.68.19.3.27   SpatialExperiment_1.22.0   
    ## [11] SingleCellExperiment_1.34.0 SummarizedExperiment_1.42.0
    ## [13] Biobase_2.72.0              GenomicRanges_1.64.0       
    ## [15] Seqinfo_1.2.0               IRanges_2.46.0             
    ## [17] S4Vectors_0.50.1            BiocGenerics_0.58.1        
    ## [19] generics_0.1.4              MatrixGenerics_1.24.0      
    ## [21] matrixStats_1.5.0           fomo_0.1.0                 
    ## [23] anndataR_1.2.0             
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] DBI_1.3.0           gridExtra_2.3.1     rlang_1.2.0        
    ##  [4] magrittr_2.0.5      otel_0.2.0          compiler_4.6.1     
    ##  [7] RSQLite_3.53.2      dir.expiry_1.20.0   png_0.1-9          
    ## [10] systemfonts_1.3.2   vctrs_0.7.3         pkgconfig_2.0.3    
    ## [13] crayon_1.5.3        fastmap_1.2.0       magick_2.9.1       
    ## [16] XVector_0.52.0      labeling_0.4.3      rmarkdown_2.31     
    ## [19] ggbeeswarm_0.7.3    ragg_1.5.2          purrr_1.2.2        
    ## [22] bit_4.6.0           xfun_0.59           cachem_1.1.0       
    ## [25] beachmat_2.28.0     jsonlite_2.0.0      blob_1.3.0         
    ## [28] rhdf5filters_1.24.0 DelayedArray_0.38.2 Rhdf5lib_2.0.0     
    ## [31] BiocParallel_1.46.0 irlba_2.3.7         parallel_4.6.1     
    ## [34] R6_2.6.1            bslib_0.11.0        RColorBrewer_1.1-3 
    ## [37] reticulate_1.46.0   jquerylib_0.1.4     Rcpp_1.1.1-1.1     
    ## [40] knitr_1.51          Matrix_1.7-5        tidyselect_1.2.1   
    ## [43] viridis_0.6.5       abind_1.4-8         yaml_2.3.12        
    ## [46] codetools_0.2-20    lattice_0.22-9      tibble_3.3.1       
    ## [49] withr_3.0.3         KEGGREST_1.52.2     S7_0.2.2           
    ## [52] evaluate_1.0.5      desc_1.4.3          Biostrings_2.80.1  
    ## [55] pillar_1.11.1       filelock_1.0.3      scales_1.4.0       
    ## [58] glue_1.8.1          tools_4.6.1         BiocNeighbors_2.6.0
    ## [61] ScaledMatrix_1.20.0 ggside_0.4.1        fs_2.1.0           
    ## [64] rhdf5_2.56.0        grid_4.6.1          basilisk_1.24.0    
    ## [67] beeswarm_0.4.0      BiocSingular_1.28.0 vipor_0.4.7        
    ## [70] cli_3.6.6           rsvd_1.0.5          rappdirs_0.3.4     
    ## [73] textshaping_1.0.5   viridisLite_0.4.3   S4Arrays_1.12.0    
    ## [76] dplyr_1.2.1         gtable_0.3.6        sass_0.4.10        
    ## [79] digest_0.6.39       SparseArray_1.12.2  ggrepel_0.9.8      
    ## [82] rjson_0.2.23        farver_2.1.2        memoise_2.0.1      
    ## [85] htmltools_0.5.9     pkgdown_2.2.0       lifecycle_1.0.5    
    ## [88] httr_1.4.8          bit64_4.8.2
