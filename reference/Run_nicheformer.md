# Run_nicheformer

This functions runs NicheFormer.

## Usage

``` r
Run_nicheformer(
  adata_path = NULL,
  technology = c("cosmx", "visium", "xenium", "merfish", "iss", "dissociated"),
  device = "cpu"
)
```

## Arguments

- adata_path:

  Path to the input anndata file.

- technology:

  Technology used for the spatial transcriptomics data. One of "cosmx",
  "visium", "xenium", "merfish", "iss", or "dissociated".

- device:

  torch device ("cpu", "cuda", "mps" etc.)
