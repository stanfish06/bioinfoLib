import scanpy as sc

from bioinfoLib.R.utils import splatter_helper, splatter_simulate_loop, start_r_session
from bioinfoLib.topology.containers import HomologyData

ro = start_r_session()
splat_func = splatter_helper(ro)
out = splatter_simulate_loop(splat_func, ro, 0.1, 1, [100], 50)
sc.pp.normalize_total(out, target_sum=1e4)
sc.pp.log1p(out)
sc.pp.pca(out)
sc.pp.neighbors(out)
sc.tl.diffmap(out)

dat = HomologyData(out.obsm["X_diffmap"])
dat.compute_homology(0.5)
dat.compute_bd_matrix(0.1)
dat.compute_loop_representatives(4, 4)

dat.boot(10, 0.1, 0.8, 0.8)
