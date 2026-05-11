#' @importFrom basilisk setBasiliskCheckVersions
.onLoad <- function(libname, pkgname) {
  basilisk::setBasiliskCheckVersions(FALSE)
}

#' @importFrom basilisk BasiliskEnvironment
.scgpt <- BasiliskEnvironment(
  pkgname="BiocFoundationGPU", 
  envname="scgpt",
  packages=c("python==3.12.13"),
  pip= c("scgpt==0.2.4",
         "torch==2.2.0",
         "ipython==9.12.0",
         "numpy==1.26.4"))

#' @importFrom basilisk BasiliskEnvironment
.novae <- BasiliskEnvironment(
  pkgname="BiocFoundationGPU", 
  envname="novae",
  packages=c("python==3.13.0"),
  pip= c("novae==1.0.4"))

#' @importFrom basilisk BasiliskEnvironment
.nimbus <- BasiliskEnvironment(
  pkgname="BiocFoundationGPU", 
  envname="nimbus",
  packages=c("python==3.10.0"),
  pip= c("Nimbus-Inference==0.0.5"))

  #' @importFrom basilisk BasiliskEnvironment
.nicheformer <- BasiliskEnvironment(
  pkgname="BiocFoundationGPU", 
  envname="nicheformer",
  packages=c("python==3.10.0"),
  pip= c("transformers==4.57.6","tiktoken==0.9.0","sentencepiece==0.2.1","git+https://github.com/theislab/nicheformer.git@485cadbc5caa15119adfd54228f8a8af835fcabc"))