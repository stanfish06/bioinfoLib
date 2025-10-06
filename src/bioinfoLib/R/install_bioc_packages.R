#!/usr/bin/env Rscript
Sys.setenv(RENV_CONFIG_SANDBOX_ENABLED = "FALSE")
options(BiocManager.check_repositories = FALSE)
Sys.setenv(RENV_CONFIG_PROMPT = "FALSE")

bioc_version="3.21"
bioc_packages <- c("edgeR", "DESeq2", "splatter", "SingleCellExperiment")
r_project <- Sys.getenv("RENV_PROJECT", unset = getwd())
setwd(r_project)
print(getwd())
library(renv)
renv::activate(r_project)
BiocManager::install(version = bioc_version, ask = FALSE, update = TRUE, force = TRUE)
BiocManager::install("BiocVersion", version = bioc_version, ask = FALSE, force = TRUE)

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE, version = bioc_version)
    message("✔ Installed Bioconductor package: ", pkg)
  } else {
    message("✔ Bioconductor package already installed: ", pkg)
  }
}

renv::snapshot(prompt = FALSE)
message("✔ renv snapshot complete")
