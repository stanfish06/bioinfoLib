import ast
import pickle
import uuid
from dataclasses import dataclass, field
from itertools import product

import numpy as np
import pandas as pd
from juliacall import Main as julia
from rich.progress import Progress
from scipy.optimize import linear_sum_assignment
from scipy.stats import binom, false_discovery_control, gamma
from sklearn.metrics import pairwise_distances

from bioinfoLib.R.utils import SimilarityMeasures_helper

from .utils import (
    compute_geometric_similarity,
    compute_homological_equivalence,
    edge_idx_encode,
    trig_idx_encode,
)


# TODO: modify this data sturcture to enable cross species matching
@dataclass
class HomologyData:
    data: np.ndarray
    n_vertices: int = 0
    data_visualization: np.ndarray = field(default_factory=lambda: np.array([]))
    persistence_diagram: np.ndarray = field(default_factory=lambda: np.array([]))
    loops_eidx: list[list[np.ndarray]] = field(default_factory=list)
    loops_coords: list[list[np.ndarray]] = field(default_factory=list)
    bd_mat: tuple = ()
    bd_column_birth_t: np.ndarray = field(
        default_factory=lambda: np.array([], dtype=float)
    )
    bd_row_id: np.ndarray = field(default_factory=lambda: np.array([], dtype=int))
    persistence_diagram_boot: list[np.ndarray] = field(default_factory=list)
    loops_eidx_boot: list[list[np.ndarray]] = field(default_factory=list)
    loops_coords_boot: list[list[np.ndarray]] = field(default_factory=list)
    matching_df: list[pd.DataFrame] = field(default_factory=list)
    tracks: dict = field(default_factory=dict)  # e.g. (1,2): [(2,3),...]
    tracks_pvals: dict = field(default_factory=dict)
    n_booted: int = 0
    loop_rank: pd.DataFrame = field(default_factory=pd.DataFrame)
    parameters: dict = field(default_factory=dict)
    metadata: dict = field(default_factory=dict)

    def __post_init__(self):
        if len(self.data_visualization) != len(self.data):
            self.data_visualization = self.data
        self.n_vertices = self.data.shape[0]

    def compute_homology(self, thresh):
        dist_mat = pairwise_distances(self.data)
        dist_mat = (dist_mat + dist_mat.T) / 2
        filt = julia.Rips(dist_mat, sparse=True, threshold=thresh)
        # don't compute representatives if only persistence diagram is needed
        result_cycle = julia.ripserer(filt, reps=False)
        birth_t = np.array([i[1] for i in result_cycle[1]])
        death_t = np.array([i[2] for i in result_cycle[1]])
        self.persistence_diagram = np.vstack([birth_t, death_t]).T
        self.parameters["filtration_threshold_homology"] = thresh

    def compute_boundary_matrix(self, thresh, reduced=True):
        if self.parameters["filtration_threshold_homology"]:
            thresh = self.parameters["filtration_threshold_homology"]
        dist_mat = pairwise_distances(self.data)
        dist_mat = (dist_mat + dist_mat.T) / 2
        filt = julia.Rips(dist_mat, sparse=True, threshold=thresh)
        edges, trigs, birth_t = julia.boundary_mat_d2(filt)
        bd = np.array(np.array(edges).tolist()) - 1
        trigs = np.array(np.array(trigs).tolist()) - 1
        trig_idx = [trig_idx_encode(i, j, k, self.n_vertices) for (i, j, k) in trigs]
        _, trigs_keep = np.unique(trig_idx, return_index=True)
        edges_keep = np.array(
            [[3 * i, 3 * i + 1, 3 * i + 2] for i in trigs_keep]
        ).flatten()
        bd = bd[edges_keep, :]
        birth_t = np.array(birth_t)[trigs_keep]

        edge_idx = [edge_idx_encode(i, j, self.n_vertices) for (i, j) in bd]
        self.bd_row_id = np.unique(edge_idx)
        one_ridx_A = np.searchsorted(self.bd_row_id, edge_idx).astype(int)
        nrow_A = len(self.bd_row_id)
        ncol_A = int(len(one_ridx_A) / 3)
        one_cidx_A = np.repeat(np.arange(ncol_A), 3).astype(int)

        self.bd_column_birth_t = birth_t
        self.bd_mat = (one_ridx_A, one_cidx_A, nrow_A, ncol_A)
        self.parameters["filtration_threshold_bd_matrix"] = thresh

    def compute_loop_representatives(
        self,
        n_top,
        n_each,
        life_pct=0.1,
        n_force_deviate=4,
        n_reps_per_loop=8,
        loop_lower_pct=5,
        loop_upper_pct=95,
        n_max_cocycles=10,
    ):
        dist_mat = pairwise_distances(self.data)
        dist_mat = (dist_mat + dist_mat.T) / 2

        filt = julia.Rips(
            dist_mat,
            sparse=True,
            threshold=self.parameters["filtration_threshold_homology"],
        )
        cocycles = julia.ripserer(filt, reps=1)
        n_total_loops = len(cocycles[1])
        n_compute = min(n_total_loops, n_top)
        self.loops_coords = []
        self.loops_eidx = []
        if n_compute > 0:
            for i in range(n_compute):
                reps = julia.reconstruct_n_loop_representatives(
                    cocycles,
                    filt,
                    i,
                    n_each,
                    life_pct,
                    n_force_deviate,
                    n_reps_per_loop,
                    loop_lower_pct,
                    loop_upper_pct,
                    n_max_cocycles,
                )
                # julia to python
                reps = [list(lp) for lp in reps[0]]
                reps_eidx = []
                reps_coords = []
                for k in range(len(reps)):
                    rep_i_idx = [j - 1 for j in reps[k]]
                    rep_i_idx.append(rep_i_idx[0])
                    rep_i_coords = []
                    rep_i_eidx = []

                    for j in range(1, len(rep_i_idx)):
                        v1 = rep_i_idx[j - 1]
                        v2 = rep_i_idx[j]
                        edge_idx = edge_idx_encode(v1, v2, self.n_vertices)
                        if edge_idx in self.bd_row_id:
                            rep_i_eidx.append(
                                np.where(edge_idx == self.bd_row_id)[0][0]
                            )
                        rep_i_coords.append(self.data[v1, :])

                    reps_eidx.append(np.array(rep_i_eidx))
                    reps_coords.append(np.array(rep_i_coords))
                self.loops_coords.append(reps_coords)
                self.loops_eidx.append(reps_eidx)

            for i in range(len(self.loops_eidx)):
                sloop_birth_t = self.persistence_diagram[i, 0]
                if f"(0,{i})" not in self.tracks:
                    self.tracks[f"(0,{i})"] = {
                        "birth_t": sloop_birth_t,
                        "loops": [(0, i)],
                    }

    # Thoughts: for approx, use permutation test to have better match, but this will be computationally intense for sure
    def boot(
        self,
        n,
        thresh,
        max_frechet_dist,
        max_hamming_dist,
        n_reps_per_loop=4,
        rep_life_pct=0.1,
        n_nearest_loops=20,
        regression_mode="exact",
        ridge_coef_a=0.1,
        ridge_coef_b=1,
        do_approximation=True,
        n_neighbors=1,
        fresh_start=True,
        n_force_deviate=4,
        _n_reps_per_loop=8,
        loop_lower_pct=5,
        loop_upper_pct=95,
        n_max_cocycles=10,
    ):
        import rpy2.robjects as ro

        similarity_func = SimilarityMeasures_helper(ro)
        if not self.bd_mat:
            raise ValueError("compute boundary matrix first")
        if not self.loops_eidx:
            raise ValueError("compute original loops first")
        if fresh_start:
            self.clean_boot()
        source_loop_eidx_pool = self.loops_eidx.copy()
        source_loop_coords_pool = self.loops_coords.copy()
        source_loop_key = []
        source_loop_birth_t = []
        source_loop_death_t = []
        for i in range(len(self.loops_eidx)):
            sloop_birth_t = self.persistence_diagram[i, 0]
            sloop_death_t = self.persistence_diagram[i, 1]
            if f"(0,{i})" not in self.tracks:
                self.tracks[f"(0,{i})"] = {"birth_t": sloop_birth_t, "loops": [(0, i)]}
            source_loop_key.append(f"(0,{i})")
            source_loop_birth_t.append(sloop_birth_t)
            source_loop_death_t.append(sloop_death_t)
        # if there are tracks/heads in tracks, send them to source loop
        if self.tracks:
            for sid in self.tracks.keys():
                # boot idx start from 1, 0 is the original batch
                batch_id, loop_id = ast.literal_eval(sid)
                if batch_id > 0:
                    source_loop_eidx_pool.append(
                        self.loops_eidx_boot[batch_id - 1][loop_id]
                    )
                    source_loop_coords_pool.append(
                        self.loops_coords_boot[batch_id - 1][loop_id]
                    )
                    source_loop_key.append(sid)
        with Progress() as progress:
            task_boot = progress.add_task("[bold #FFA500]Booting", total=n)
            task_find_loop = progress.add_task("[bold green]Find loops")
            task_frechet = progress.add_task("[bold green]Calculate Frechet distance")
            task_homology = progress.add_task(
                "[bold green]Assess Homological equivalence"
            )
            n_success = 0
            for i in range(n):
                boot_idx = np.random.choice(
                    self.data.shape[0],
                    size=self.data.shape[0],
                    replace=True,
                )
                x_boot = self.data[boot_idx]
                x_boot = x_boot + np.random.normal(scale=0.0001, size=self.data.shape)
                dist_mat = pairwise_distances(x_boot)
                dist_mat = (dist_mat + dist_mat.T) / 2
                filt = julia.Rips(dist_mat, sparse=True, threshold=thresh)
                cocycles = julia.ripserer(filt, reps=1)
                birth_t = np.array([i[1] for i in cocycles[1]])
                death_t = np.array([i[2] for i in cocycles[1]])
                reps_eidx_boot = []
                reps_coord_boot = []
                progress.reset(task_find_loop, totol=len(cocycles[1]))
                for nc in range(len(cocycles[1])):
                    try:
                        reps = julia.reconstruct_n_loop_representatives(
                            cocycles,
                            filt,
                            nc,
                            n_reps_per_loop,
                            rep_life_pct,
                            n_force_deviate,
                            _n_reps_per_loop,
                            loop_lower_pct,
                            loop_upper_pct,
                            n_max_cocycles,
                        )
                        reps = [list(lp) for lp in reps[0]]
                        reps_eidx = []
                        reps_coords = []
                        for i in range(len(reps)):
                            rep_i_idx = [j - 1 for j in reps[i]]
                            rep_i_idx.append(rep_i_idx[0])
                            rep_i_coords = []
                            rep_i_eidx = []

                            for j in range(1, len(rep_i_idx)):
                                v1 = boot_idx[rep_i_idx[j - 1]]
                                v2 = boot_idx[rep_i_idx[j]]
                                edge_idx = edge_idx_encode(v1, v2, self.n_vertices)
                                if edge_idx in self.bd_row_id:
                                    rep_i_eidx.append(
                                        np.where(edge_idx == self.bd_row_id)[0][0]
                                    )
                                rep_i_coords.append(self.data[v1, :])
                            reps_eidx.append(np.array(rep_i_eidx))
                            reps_coords.append(np.array(rep_i_coords))
                        reps_eidx_boot.append(reps_eidx)
                        reps_coord_boot.append(reps_coords)
                    except ValueError:
                        continue
                    progress.update(task_find_loop, advance=1)
                pairs = [
                    ((i, j), source_loop_death_t[i])
                    for i, j in product(
                        range(len(source_loop_eidx_pool)), range(len(reps_eidx_boot))
                    )
                ]
                progress.reset(task_frechet, totol=len(pairs))
                result = []
                for (i, j), sloop_death_t in pairs:
                    result.append(
                        compute_geometric_similarity(
                            source_loops_coords=source_loop_coords_pool[i],
                            target_loops_coords=reps_coord_boot[j],
                            max_frechet_dist=max_frechet_dist,
                            similarity_func=similarity_func,
                            similarity_type="Frechet",
                        )
                    )
                    progress.update(task_frechet, advance=1)
                pairs_filt = []
                for si in range(len(source_loop_eidx_pool)):
                    n_best_matches = []
                    for j in range(len(reps_eidx_boot)):
                        k = si * len(reps_eidx_boot) + j
                        if len(n_best_matches) < n_nearest_loops:
                            if result[k] is not None:
                                n_best_matches.append([k, result[k]])
                                n_best_matches.sort(key=lambda x: x[1])
                        else:
                            if result[k] is not None:
                                if result[k] < n_best_matches[-1][1]:
                                    n_best_matches[-1] = [k, result[k]]
                                    n_best_matches.sort(key=lambda x: x[1])
                    pairs_filt.extend(
                        [(pairs[k], frech_dist) for k, frech_dist in n_best_matches]
                    )
                if pairs_filt:
                    result = []
                    progress.reset(task_homology, total=len(pairs_filt))
                    for ((i, j), sloop_death_t), frech_dist in pairs_filt:
                        result.append(
                            compute_homological_equivalence(
                                source_loops_edges=source_loop_eidx_pool[i],
                                target_loops_edges=reps_eidx_boot[j],
                                one_ridx_A=self.bd_mat[0],
                                one_cidx_A=self.bd_mat[1],
                                nrow_A=self.bd_mat[2],
                                ncol_A=self.bd_mat[3],
                                ridge_coef_a=ridge_coef_a,
                                ridge_coef_b=ridge_coef_b,
                                do_approximation=do_approximation,
                                n_neighbors=n_neighbors,
                                bd_column_birth_t=self.bd_column_birth_t,
                                source_loop_death_t=sloop_death_t,
                                regression_mode=regression_mode,
                            )
                        )
                        progress.update(task_homology, advance=1)

                    df = pd.DataFrame(pairs_filt, columns=["pair", "frechet_dist"])
                    df["hamming_dist"] = result
                    self.matching_df.append(df.copy())
                    df = df[
                        np.logical_and(
                            df["hamming_dist"] < max_hamming_dist,
                            ~np.isnan(df["hamming_dist"]),
                        )
                    ]
                    if df.empty:
                        continue
                    self.loops_eidx_boot.append(reps_eidx_boot)
                    self.loops_coords_boot.append(reps_coord_boot)
                    self.persistence_diagram_boot.append(
                        np.vstack([birth_t, death_t]).T
                    )
                    df[["source", "target"]] = np.array([p for p, _ in df["pair"]])
                    df["cost"] = df["hamming_dist"] * df["frechet_dist"]
                    cost_matrix_sub = df.pivot(
                        index="source", columns="target", values="cost"
                    ).fillna(np.inf)
                    cost_matrix = np.full(
                        [
                            cost_matrix_sub.shape[0],
                            cost_matrix_sub.shape[0] + cost_matrix_sub.shape[1],
                        ],
                        max_hamming_dist * max_frechet_dist,
                    )
                    cost_matrix[
                        : cost_matrix_sub.shape[0], : cost_matrix_sub.shape[1]
                    ] = cost_matrix_sub
                    row_ind, col_ind = linear_sum_assignment(cost_matrix)
                    mask = col_ind >= cost_matrix_sub.shape[1]
                    row_ind = cost_matrix_sub.index[row_ind[~mask]]
                    col_ind = cost_matrix_sub.columns[col_ind[~mask]]
                    for j in range(len(row_ind)):
                        skey = source_loop_key[row_ind[j]]
                        self.tracks[skey]["loops"].append(
                            (self.n_booted + 1, col_ind[j])
                        )
                    n_success += 1
                    self.n_booted += 1
                progress.update(
                    task_boot,
                    advance=1,
                    description=f"[bold #FFA500]Booting (DONE={n_success}/{n})",
                )
                progress.reset(task_find_loop, completed=0)
                progress.reset(task_frechet, completed=0)
                progress.reset(task_homology, completed=0)

    def clean_boot(self):
        self.tracks = {}
        self.n_booted = 0
        self.loops_coords_boot = []
        self.loops_eidx_boot = []
        self.persistence_diagram_boot = []
        self.matching_df = []

    def write_pkl(self, fname=None):
        if fname is not None:
            with open(f"{fname}.pkl", "wb") as f:
                pickle.dump(self, f)
        else:
            fname = str(uuid.uuid4().hex)
            with open(f"{fname}.pkl", "wb") as f:
                pickle.dump(self, f)

    def rank_loops(self):
        if not self.n_booted:
            raise ValueError("do bootstrapping first")
        presence_probs = [
            (
                ast.literal_eval(src_loop),
                len(track["loops"]),
                (len(track["loops"]) - 1) / self.n_booted,
            )
            for src_loop, track in self.tracks.items()
        ]
        mean_presence_prob = np.mean([p for _, _, p in presence_probs])
        presence_pvals = [
            (src_loop, p, 1 - binom.cdf(n, self.n_booted, mean_presence_prob))
            for src_loop, n, p in presence_probs
        ]
        # null distribution for lifetime
        lifetimes = [self.persistence_diagram[:, 1] - self.persistence_diagram[:, 0]]
        for boot_ph in self.persistence_diagram_boot:
            lifetimes.append(boot_ph[:, 1] - boot_ph[:, 0])
        lifetimes_full = np.concatenate(lifetimes)
        lifetimes_full = lifetimes_full[lifetimes_full < np.inf]
        params = gamma.fit(lifetimes_full, floc=0)
        persistence_pvals = [
            (
                lifetimes[src_loop[0]][src_loop[1]],
                1
                - gamma.cdf(
                    lifetimes[src_loop[0]][src_loop[1]],
                    params[0],
                    loc=params[1],
                    scale=params[2],
                ),
            )
            for src_loop, _, _ in presence_probs
        ]
        self.parameters["gamma_fit"] = params
        self.parameters["binom_fit"] = mean_presence_prob
        df = pd.DataFrame(
            presence_pvals, columns=["src_loop", "prob_presence", "pval_presence"]
        )
        df[["src_lifetime", "pval_persistence"]] = persistence_pvals
        df["pval_presence_adjust"] = false_discovery_control(df["pval_presence"])
        df["pval_persistence_adjust"] = false_discovery_control(df["pval_persistence"])
        self.loop_rank = df
