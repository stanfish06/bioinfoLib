import sys

import numpy as np
from anndata import AnnData


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


def splatter_helper(robjects):
    # Get the path to the Julia helper file
    import os

    current_dir = os.path.dirname(os.path.abspath(__file__))
    splatter_file = os.path.join(current_dir, "splatter_helper.R")

    # Load the Julia functions from file
    robjects.r(f'source("{splatter_file}")')

    splatter_helper_func = {
        "splatSimulateLoop": robjects.r["splatSimulateLoop"],
        "newSplatParams": robjects.r["newSplatParams"],
        "assay": robjects.r["assay"],
        "colData": robjects.r["colData"],
    }
    return splatter_helper_func


def splatter_simulate_loop(
    splatter_helper_func, robjects, bcv_common, batch_facLoc, batchCells, nSteps
):
    from rpy2.robjects import pandas2ri

    newSplatParams = {"bcv.common": bcv_common, "batch.facLoc": batch_facLoc}
    newSplatParams = splatter_helper_func["newSplatParams"](**newSplatParams)
    simulationParams = {"batchCells": robjects.IntVector(batchCells), "nSteps": nSteps}
    out = splatter_helper_func["splatSimulateLoop"](newSplatParams, **simulationParams)
    count_mat = np.array(splatter_helper_func["assay"](out, "counts")).T
    meta = robjects.r["as.data.frame"](splatter_helper_func["colData"](out))
    with (robjects.default_converter + pandas2ri.converter).context():
        obs = robjects.conversion.get_conversion().rpy2py(meta)
    out = AnnData(count_mat, obs=obs)
    return out
