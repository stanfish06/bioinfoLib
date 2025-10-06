import os
import subprocess
import sys

import numpy as np
from anndata import AnnData

base_dir = os.path.dirname(os.path.abspath(__file__))
r_home = os.environ["R_HOME"]


def install_cran_packages():
    cmd = (
        os.path.join(r_home, "bin", "Rscript"),
        "-e",
        ";".join(
            (
                f'setwd("{base_dir}")',
                'install.packages("renv")',
                "library(renv)",
                'if (!file.exists("renv.lock")) renv::init(bare = TRUE, bioconductor = "3.18")',
                f'renv::activate("{base_dir}")',
                "renv::snapshot(prompt = FALSE)",
            )
        ),
    )
    subprocess.check_call(cmd)
    import rpy2.robjects as robjects

    # None will just install latest
    R_PACKAGES = {
        "tidyverse": None,
        "BiocManager": "3.18",
        "ggplot2": None,
        "SimilarityMeasures": None,
    }
    robjects.r("library(renv)")
    robjects.r(f'renv::activate("{base_dir}")')
    installed_packages = list(robjects.r("rownames(installed.packages())"))
    print(installed_packages)
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


def install_bioc_packages():
    os.environ["RENV_PROJECT"] = base_dir
    cmd = (
        os.path.join(r_home, "bin", "Rscript"),
        os.path.join(base_dir, "install_bioc_packages.R"),
    )
    subprocess.check_call(cmd)


_np_cv_rules = None
_pd_cv_rules = None


def start_r_session():
    import rpy2.robjects as robjects
    from rpy2.robjects import default_converter, numpy2ri, pandas2ri

    global _np_cv_rules, _pd_cv_rules
    _np_cv_rules = default_converter + numpy2ri.converter
    _pd_cv_rules = default_converter + pandas2ri.converter

    robjects.r(f'setwd("{base_dir}")')
    robjects.r("library(renv)")
    robjects.r(f'renv::activate("{base_dir}")')
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


def SimilarityMeasures_helper(robjects):
    robjects.r("library(SimilarityMeasures)")
    SimilarityMeasures_helper_func = {
        "Frechet": robjects.r["Frechet"],
    }
    return SimilarityMeasures_helper_func


def trajectory_distance(curve1, curve2, type, similarity_helper_func):
    with _np_cv_rules.context():
        dist = similarity_helper_func[type](curve1, curve2)
    return dist


def splatter_simulate_loop(
    splatter_helper_func, robjects, bcv_common, batch_facLoc, batchCells, nSteps
):
    newSplatParams = {"bcv.common": bcv_common, "batch.facLoc": batch_facLoc}
    newSplatParams = splatter_helper_func["newSplatParams"](**newSplatParams)
    simulationParams = {"batchCells": robjects.IntVector(batchCells), "nSteps": nSteps}
    out = splatter_helper_func["splatSimulateLoop"](newSplatParams, **simulationParams)
    count_mat = np.array(splatter_helper_func["assay"](out, "counts")).T
    meta = robjects.r["as.data.frame"](splatter_helper_func["colData"](out))
    with _pd_cv_rules.context():
        obs = robjects.conversion.get_conversion().rpy2py(meta)
    out = AnnData(count_mat, obs=obs)
    return out
