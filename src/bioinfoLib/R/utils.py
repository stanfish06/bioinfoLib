import sys


def install_cran_packages():
    import rpy2.robjects as robjects

    # None will just install latest
    R_PACKAGES = {"tidyverse": None, "BiocManager": "3.17", "ggplot2": None}
    robjects.r("library(renv)")
    robjects.r('if (!file.exists("renv.lock")) renv::init(bare = TRUE)')
    robjects.r("renv::activate()")
    installed_packages = list(robjects.r("rownames(installed.packages())"))
    for pkg, version in R_PACKAGES.items():
        # some packages might be updated by biocmanager, so if they have been installed and work properly, just let them go.
        if pkg in installed_packages:
            continue
        try:
            if version:
                robjects.r(f'renv::install("{pkg}@{version}")')
            else:
                robjects.r(f'renv::install("{pkg}")')
            print(f"✔ R package {pkg} ({version}) installed")
        except Exception as e:
            print(f"✘ Failed to install {pkg}: {e}")
    robjects.r("renv::snapshot(prompt = FALSE)")


def start_r_session():
    import rpy2.robjects as robjects

    robjects.r("library(renv)")
    robjects.r("renv::activate()")
    python_interpretor = sys.executable
    robjects.r(f'renv::use_python("{python_interpretor}")')
    return robjects
