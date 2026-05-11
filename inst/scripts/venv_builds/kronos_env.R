library(fomo)
library(reticulate)

# setup for nicheformer
use_python(
  install_python(version=ver <- "3.13.0"),
  required = FALSE
)

# get requirements
req <- c(
  "setuptools==71.0.0",
  "git+https://github.com/mahmoodlab/KRONOS.git"
)

# create and install virtual env
virtualenv_create("kronos_env", python = ver)
virtualenv_install("kronos_env", packages = req)
