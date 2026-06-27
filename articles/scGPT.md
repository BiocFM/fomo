# Running scGPT with fomo

## Introduction

This vignette demonstrates how to use `Run_scGPT` to embed single-cell
RNA-seq data using the scGPT model(Cui et al. 2024).

scGPT is a foundation model for single-cell mRNA-Seq data. scGPT can be
used in many different context and strategies (see their official
[github](https://github.com/bowang-lab/scGPT) repository). In
particular, scGPT provides multiple possible pretrained models which are
linked to in their [github
README](https://github.com/bowang-lab/scGPT/blob/main/README.md). The
models differ based on the training data used. The most common choice of
model is the whole-human
[`scGPT_human`](https://drive.google.com/drive/folders/1oWh_-ZRdhtoGQ2Fw24HP41FgLoomVo-y),
which is trained 33 million normal human cells across a range of tissue
types. Tissue-type specific models are available. There is also a
“continual pretrained model”,
[`scGPT_CP`](https://drive.google.com/drive/folders/1_GROJTzXiAV8HB4imruOTk6PEGuNOcgB)
that inherits the pre-trained scGPT whole-human model, and was further
supervised by extra cell type labels (using the Tabula Sapiens dataset).

In `Run_scGPT` we run a basic “zero-shot” embedding of the data via the
`scg.tasks.embed_data()` function in the `scgpt` python package.
“Zero-shot” means that the input data is run through the pretrained
model, and the output embeddings are returned. Depending on which model
is used, this could correspond to the [zero-shot tutorials of
scGPT](https://github.com/bowang-lab/scGPT/tree/main/tutorials/zero-shot).
Those tutorials also provide strategies for how to use these to perform
integration of datasets.

### Huggingface

scGPT also provides access to their model [huggingface
model](https://huggingface.co/tdc/scGPT). We do not currently provide
implementation of this model, since scGPT offers many more models than
just that available on huggingface. We have focused on an implementation
that allows users the ability to work with and compare these many
different models. However, this means users must download the models and
store them in a local directory for `Run_scGPT`.

## Prerequisites

`Run_scGPT` takes a `.h5ad` filename containing the cell transcripts
data and a pretrained scGPT model directory as input, and returns a
matrix of cell embeddings.

### Input data

The input `.h5ad` file used in this vignette is the example data
provided by scGPT (batch_covid_subsampled_train.h5ad
)\[<https://github.com/bowang-lab/scGPT/blob/main/tutorials/zero-shot/Tutorial_ZeroShot_Reference_Mapping.ipynb>\],
but we have subsetted to the first 100 cells so that it will run
instantaneously even on a CPU. The subsetted 100 cells can be downloaded
from:
<https://drive.google.com/file/d/1Do7CXaaSTwEySGWHMKAkpN5G_g-OY8y2/view?usp=drive_link>

### Model weights

The pretrained scGPT model directory that we use is `scGPT_human` model
which can be downloaded from:
<https://drive.google.com/drive/folders/1oWh_-ZRdhtoGQ2Fw24HP41FgLoomVo-y>

The model directory should contain a `.pt` file with the model weights,
along with `vocab.json` and `args.json`.

## Running scGPT

Set the paths to the input data and model directory:

``` r

dir <- "/path/to/Hackathon/"
h5ad_file <- file.path(dir, "scHuman_covid/batch_covid_subsampled_test_100cells.h5ad")
model_dir  <- file.path(dir, "scGPT_model")
```

Run `Run_scGPT` to embed the cells. The result is a matrix with one row
per cell and one column per embedding dimension:

``` r

library(fomo)

result <- Run_scGPT(
    h5ad_file = h5ad_file,
    model_dir  = model_dir,
    gene_col   = "gene_name"
)
```

## Basic Python script

The basic python script that underlies this function is:

    import scanpy as sc
    import scgpt as scg
    adata = sc.read_h5ad(h5ad_file)

    ref_embed_adata = scg.tasks.embed_data(
        adata,
        model_dir,
        gene_col=gene_col,
        batch_size=64,
    )

There is additional code added to ensure that it runs on MacOS.

## Session Info

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
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] BiocStyle_2.40.0
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] digest_0.6.39       desc_1.4.3          R6_2.6.1           
    ##  [4] bookdown_0.47       fastmap_1.2.0       xfun_0.59          
    ##  [7] cachem_1.1.0        knitr_1.51          htmltools_0.5.9    
    ## [10] rmarkdown_2.31      lifecycle_1.0.5     cli_3.6.6          
    ## [13] sass_0.4.10         pkgdown_2.2.0       textshaping_1.0.5  
    ## [16] jquerylib_0.1.4     systemfonts_1.3.2   compiler_4.6.1     
    ## [19] tools_4.6.1         ragg_1.5.2          bslib_0.11.0       
    ## [22] evaluate_1.0.5      yaml_2.3.12         BiocManager_1.30.27
    ## [25] otel_0.2.0          jsonlite_2.0.0      rlang_1.2.0        
    ## [28] fs_2.1.0

## References

Cui, Haotian, Chloe Wang, Hassaan Maan, et al. 2024. “scGPT: Toward
Building a Foundation Model for Single-Cell Multi-Omics Using Generative
AI.” *Nature Methods* 21: 1470–80.
<https://doi.org/10.1038/s41592-024-02201-0>.
