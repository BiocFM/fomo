library(fomo)
library(reticulate)

# setup for nicheformer
use_python(
  install_python(version=ver <- "3.13.0"),
  required = FALSE
)

# DOES NOT WORK NOW
# get requirements
req <- c(
  "setuptools==71.0.0",
  "git+https://github.com/theislab/nicheformer.git"
)

# create and install virtual env
virtualenv_create("nicheformer_env", python = ver)
virtualenv_install("nicheformer_env", packages = "setuptools==71.0.0")