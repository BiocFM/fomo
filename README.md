## fomo

This package allows deploying a set of pre-defined foundation models for various tasks in R. It provides an interface to easily access and utilize these models using:

-   `basilisk` (<https://bioconductor.org/packages/release/bioc/html/basilisk.html>)
-   `reticulate` (<https://rstudio.github.io/reticulate/>)

and for tasks such as generating embeddings of spatial or imaging datasets.

The package is designed to be user-friendly and efficient, making it accessible to both beginners and experienced users in the field of machine learning and deep learning.

Both `basilisk` and `reticulate` are used to manage Python dependencies and provide a seamless interface between R and Python. However, for foundation models that are deployed with `reticulate` only (e.g. KRONOS), example scripts for building and running environments can be found under `inst/scripts/venv_builds`.

Currently there are environments for the following models:

-   scGPT (<https://github.com/bowang-lab/scGPT>)
-   Novae (<https://github.com/prism-oncology/novae>)
-   KRONOS (<https://github.com/mahmoodlab/KRONOS>)
-   Nicheformer (<https://github.com/theislab/nicheformer>)
