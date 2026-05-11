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
  "git+https://github.com/prov-gigapath/prov-gigapath.git"
)

# create and install virtual env
conda_create("gigapath_env", 
             environment = system.file("scripts/venv_builds/environment.yml", 
                                       package = "fomo"))