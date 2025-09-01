#!/usr/bin/env Rscript
Sys.setenv(RENV_CONFIG_SANDBOX_ENABLED = "FALSE")
options(BiocManager.check_repositories = FALSE)
Sys.setenv(RENV_CONFIG_PROMPT = "FALSE")

bioc_packages <- c("edgeR", "DESeq2", "splatter", "SingleCellExperiment")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
library(renv)
renv::activate()
BiocManager::install(version = "3.17", ask = FALSE, update = TRUE, force = TRUE)
BiocManager::install("BiocVersion", version="3.17", ask=FALSE, force=TRUE)

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE, version = "3.17")
    message("✔ Installed Bioconductor package: ", pkg)
  } else {
    message("✔ Bioconductor package already installed: ", pkg)
  }
}

renv::snapshot(prompt = FALSE)
message("✔ renv snapshot complete")
