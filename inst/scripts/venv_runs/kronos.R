library(fomo)
library(reticulate)

reticulate::use_virtualenv("kronos_env")
kr <- reticulate::import("kronos")