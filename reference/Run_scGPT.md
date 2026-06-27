# scGPT

This function runs scGPT.

## Usage

``` r
Run_scGPT(h5ad_file, model_dir, gene_col, batch_size = 64L)
```

## Arguments

- h5ad_file:

  file that points to the location of a h5ad file with the count data

- model_dir:

  directory with the model. This should be a directory that the scGPT
  model, consisting of a \`.pt\` file with the model, as well as the
  \`vocab.json\` and \`args.json\` file.

- gene_col:

  The gene column name

- batch_size:

  The batch size passed to \`scg.tasks.embed_data\` (default 64)

## Details

This function runs the \`scg.tasks.embed_data\` function
