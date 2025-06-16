import ast
import pickle
import uuid
from dataclasses import dataclass, field
from itertools import product
from typing import Tuple

import numpy as np
import pandas as pd
from juliacall import Main as julia
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist
from scipy.stats import binom, false_discovery_control, gamma
from sklearn.metrics import pairwise_distances
from tqdm import tqdm

from .utils import edge_idx_encode, evaluate_match


@dataclass
class HomologyData:
    data_original: np.ndarray
    sample_label: str = ""
    n_true_loops: int = 0
    loop_size: list[float] = field(default_factory=list)
    noise_level: float = 0.0
    n_points: int = 0
    persistence_diagram_original: np.ndarray = field(
        default_factory=lambda: np.array([])
    )
    loops_eidx_original: list[np.ndarray] = field(default_factory=list)
    loops_coords_original: list[np.ndarray] = field(default_factory=list)
    bd_mat: Tuple = ()
    bd_column_birth_t: np.ndarray = field(
        default_factory=lambda: np.array([], dtype=float)
    )
    bd_row_id: np.ndarray = field(default_factory=lambda: np.array([], dtype=int))
    data_boot: list[np.ndarray] = field(default_factory=list)
    persistence_diagram_boot: list[np.ndarray] = field(default_factory=list)
    loops_eidx_boot: list[np.ndarray] = field(default_factory=list)
    loops_coords_boot: list[np.ndarray] = field(default_factory=list)
    matching_df: list[pd.DataFrame] = field(default_factory=list)
    tracks: dict = field(default_factory=dict)  # e.g. (1,2): [(2,3),...]
    tracks_pvals: dict = field(default_factory=dict)
    n_booted: int = 0
    loop_rank: pd.DataFrame = field(default_factory=pd.DataFrame)
    parameters: dict = field(default_factory=dict)

    def compute_original_homology(self, thresh):
        dist_mat = pairwise_distances(self.data_original)
        dist_mat = (dist_mat + dist_mat.T) / 2
        filt = julia.Rips(dist_mat, sparse=True, threshold=thresh)
        # don't compute representatives if only persistence diagram is needed
        result_cycle = julia.ripserer(filt, reps=False)
        birth_t = np.array([i[1] for i in result_cycle[1]])
        death_t = np.array([i[2] for i in result_cycle[1]])
        self.persistence_diagram_original = np.vstack([birth_t, death_t]).T
        self.parameters["max_filtration_PH_original"] = thresh

    def compute_bd_matrix(self, thresh):
        dist_mat = pairwise_distances(self.data_original)
        dist_mat = (dist_mat + dist_mat.T) / 2
        filt = julia.Rips(dist_mat, sparse=True, threshold=thresh)
        bd, birth_t = julia.boundary_mat_d2(filt)
        bd = np.array(np.array(bd).tolist()) - 1
        self.bd_column_birth_t = np.array(birth_t)

        edge_idx = [edge_idx_encode(i, j) for (i, j) in bd]
        self.bd_row_id = np.unique(edge_idx)
        one_ridx_A = np.searchsorted(self.bd_row_id, edge_idx).astype(int)
        nrow_A = len(self.bd_row_id)
        ncol_A = int(len(one_ridx_A) / 3)
        one_cidx_A = np.repeat(np.arange(ncol_A), 3).astype(int)
        self.bd_mat = (one_ridx_A, one_cidx_A, nrow_A, ncol_A)
        self.parameters["filtration_bd_matrix"] = thresh

    def compute_original_loops(self, thresh):
        if not self.bd_mat:
            raise ValueError("compute boundary matrix first")
        else:
            self.loops_eidx_original = []
            self.loops_coords_original = []
            self.parameters["max_filtration_PH_representative"] = thresh
            dist_mat = pairwise_distances(self.data_original)
            dist_mat = (dist_mat + dist_mat.T) / 2
            filt = julia.Rips(dist_mat, sparse=True, threshold=thresh)
            result_cycle = julia.ripserer(filt, reps=1)
            for n in range(len(result_cycle[1])):
                rep = result_cycle[1][n]
                rep = julia.Ripserer.reconstruct_cycle(filt, rep, rep.birth)
                rep_ridx = []
                for i, j in rep:
                    edge_idx = edge_idx_encode(i - 1, j - 1)
                    if edge_idx in self.bd_row_id:
                        rep_ridx.append(np.where(edge_idx == self.bd_row_id)[0][0])
                if len(rep_ridx) > 0:
                    self.loops_eidx_original.append(np.array(rep_ridx))
                    loop_i_coords = []
                    for i, (v1, v2) in enumerate(rep):
                        if i == 0:
                            loop_i_coords.extend(
                                self.data_original[[v1 - 1, v2 - 1], :]
                            )
                        elif i == 1:
                            e0 = loop_i_coords[:2]
                            dist_e0_e1 = cdist(
                                e0, self.data_original[[v1 - 1, v2 - 1], :]
                            )
                            i0, j0 = np.where(dist_e0_e1 == 0)
                            if i0 == 0:
                                loop_i_coords = loop_i_coords[::-1]
                            if j0 == 1:
                                loop_i_coords.append(self.data_original[v1 - 1, :])
                            else:
                                loop_i_coords.append(self.data_original[v2 - 1, :])
                        else:
                            e0 = np.expand_dims(loop_i_coords[-1], 0)
                            dist_e0_e1 = cdist(
                                e0, self.data_original[[v1 - 1, v2 - 1], :]
                            )
                            i0, j0 = np.where(dist_e0_e1 == 0)
                            if j0 == 1:
                                loop_i_coords.append(self.data_original[v1 - 1, :])
                            else:
                                loop_i_coords.append(self.data_original[v2 - 1, :])
                    self.loops_coords_original.append(np.array(loop_i_coords))

    def clean_boot(self):
        self.tracks = {}
        self.n_booted = 0
        self.data_boot = []
        self.loops_coords_boot = []
        self.loops_eidx_boot = []
        self.persistence_diagram_boot = []

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
        lifetimes = [
            self.persistence_diagram_original[:, 1]
            - self.persistence_diagram_original[:, 0]
        ]
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

    def boot(
        self,
        n,
        thresh,
        max_frechet_dist,
        max_hamming_dist,
        n_nearest_loops=20,
        n_search=4,
        ridge_coef_a=0.1,
        ridge_coef_b=0.1,
        do_approximation=True,
        n_neighbors=1,
        fresh_start=False,
        verbose=False,
    ):
        if not self.bd_mat:
            raise ValueError("compute boundary matrix first")
        if not self.loops_eidx_original:
            raise ValueError("compute original loops first")
        if fresh_start:
            self.clean_boot()
        source_loop_eidx_pool = self.loops_eidx_original.copy()
        source_loop_coords_pool = self.loops_coords_original.copy()
        source_loop_key = []
        source_loop_birth_t = []
        for i in range(len(self.loops_eidx_original)):
            sloop_birth_t = self.persistence_diagram_original[i, 0]
            if f"(0,{i})" not in self.tracks:
                self.tracks[f"(0,{i})"] = {"birth_t": sloop_birth_t, "loops": [(0, i)]}
            source_loop_key.append(f"(0,{i})")
            source_loop_birth_t.append(sloop_birth_t)
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

        for i in range(n):
            boot_idx = np.random.choice(
                self.data_original.shape[0],
                size=self.data_original.shape[0],
                replace=True,
            )
            x_boot = self.data_original[boot_idx]
            x_boot = x_boot + np.random.normal(
                scale=0.0001, size=self.data_original.shape
            )
            dist_mat = pairwise_distances(x_boot)
            dist_mat = (dist_mat + dist_mat.T) / 2
            filt = julia.Rips(dist_mat, sparse=True, threshold=thresh)
            result_cycle = julia.ripserer(filt, reps=1)
            birth_t = np.array([i[1] for i in result_cycle[1]])
            death_t = np.array([i[2] for i in result_cycle[1]])
            loop_eidx_boot = []
            rep_all = []
            success = False
            try:
                for nc in range(len(result_cycle[1])):
                    rep = result_cycle[1][nc]
                    rep = julia.Ripserer.reconstruct_cycle(filt, rep, rep.birth)
                    rep_ridx = [
                        np.where(
                            edge_idx_encode(boot_idx[i - 1], boot_idx[j - 1])
                            == self.bd_row_id
                        )[0][0]
                        for (i, j) in rep
                        if edge_idx_encode(boot_idx[i - 1], boot_idx[j - 1])
                        in self.bd_row_id
                    ]
                    loop_eidx_boot.append(rep_ridx)
                    rep_raw = np.vstack(
                        [
                            np.vstack(
                                [
                                    x_boot[[i - 1, j - 1], :],
                                    np.repeat(np.nan, x_boot.shape[1]),
                                ]
                            )
                            for (i, j) in rep
                        ]
                    )
                    rep_raw = []
                    for i, (v1, v2) in enumerate(rep):
                        if i == 0:
                            rep_raw.extend(x_boot[[v1 - 1, v2 - 1], :])
                        elif i == 1:
                            e0 = x_boot[:2]
                            dist_e0_e1 = cdist(e0, x_boot[[v1 - 1, v2 - 1], :])
                            i0, j0 = np.where(dist_e0_e1 == 0)
                            if i0 == 0:
                                rep_raw = rep_raw[::-1]
                            if j0 == 1:
                                rep_raw.append(x_boot[v1 - 1, :])
                            else:
                                rep_raw.append(x_boot[v2 - 1, :])
                        else:
                            e0 = np.expand_dims(rep_raw[-1], 0)
                            dist_e0_e1 = cdist(e0, x_boot[[v1 - 1, v2 - 1], :])
                            i0, j0 = np.where(dist_e0_e1 == 0)
                            if j0 == 1:
                                rep_raw.append(x_boot[v1 - 1, :])
                            else:
                                rep_raw.append(x_boot[v2 - 1, :])
                    rep_all.append(np.array(rep_raw))
                success = True
            except ValueError:
                continue
            if success:

                def evaluate_match_tmp(p, do_regression, source_loop_birth_t=None):
                    return evaluate_match(
                        max_frechet_dist=max_frechet_dist,
                        n_search=n_search,
                        ncol_A=self.bd_mat[3],
                        nrow_A=self.bd_mat[2],
                        one_cidx_A=self.bd_mat[1],
                        one_ridx_A=self.bd_mat[0],
                        ridge_coef_a=ridge_coef_a,
                        ridge_coef_b=ridge_coef_b,
                        source_loop_coords=source_loop_coords_pool[p[0]],
                        target_loop_coords=rep_all[p[1]],
                        source_loop_edges=source_loop_eidx_pool[p[0]],
                        target_loop_edges=loop_eidx_boot[p[1]],
                        do_regression=do_regression,
                        do_approximation=do_approximation,
                        n_neighbors=n_neighbors,
                        bd_column_birth_t=self.bd_column_birth_t,
                        source_loop_birth_t=source_loop_birth_t,
                    )

                # print("start initial matching", flush=True)
                pairs = [
                    ((i, j), source_loop_birth_t[i])
                    for i, j in product(
                        range(len(source_loop_eidx_pool)), range(len(loop_eidx_boot))
                    )
                ]
                if verbose:
                    result = [evaluate_match_tmp(p, False) for p, _ in tqdm(pairs)]
                else:
                    result = [evaluate_match_tmp(p, False) for p, _ in pairs]
                # print("done initial matching", flush=True)
                pairs_filt = []
                for si in range(len(source_loop_eidx_pool)):
                    n_best_matches = []
                    for j in range(len(loop_eidx_boot)):
                        k = si * len(loop_eidx_boot) + j
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
                        [(pairs[k], haus_dist) for k, haus_dist in n_best_matches]
                    )
                self.loops_eidx_boot.append(loop_eidx_boot)
                self.loops_coords_boot.append(rep_all)
                self.persistence_diagram_boot.append(np.vstack([birth_t, death_t]).T)
                if pairs_filt:
                    if verbose:
                        result = [
                            evaluate_match_tmp(p[0], True, p[1])
                            for p, _ in tqdm(pairs_filt)
                        ]
                    else:
                        result = [
                            evaluate_match_tmp(p[0], True, p[1]) for p, _ in pairs_filt
                        ]
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
                        self.n_booted = self.n_booted + 1
                        continue
                    df[["source", "target"]] = np.array([p for p, _ in df["pair"]])
                    # seems problematic
                    # df = (
                    #     df.sort_values(by=["hamming_dist", "frechet_dist"])
                    #     .groupby("target", as_index=False)
                    #     .first()
                    #     .sort_values(by=["hamming_dist", "frechet_dist"])
                    #     .groupby("source", as_index=False)
                    #     .first()
                    # )
                    cost_matrix_sub = df.pivot(
                        index="source", columns="target", values="frechet_dist"
                    ).fillna(np.inf)
                    cost_matrix = np.full(
                        [
                            cost_matrix_sub.shape[0],
                            cost_matrix_sub.shape[0] + cost_matrix_sub.shape[1],
                        ],
                        max_frechet_dist,
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
                    # unmatched_targets = np.setdiff1d(
                    #     list(range(len(loop_eidx_boot))), df["target"]
                    # )
                    # for j in range(df.shape[0]):
                    #     skey = source_loop_key[df["source"][j]]
                    #     self.tracks[skey]["loops"].append(
                    #         (self.n_booted + 1, df["target"][j])
                    #     )
                    # for unmatched target, we make them new heads
                    # for j in unmatched_targets:
                    #     self.tracks[f"({self.n_booted + 1},{j})"] = {
                    #         "birth_t": birth_t[j],
                    #         "loops": [(self.n_booted + 1, j)],
                    #     }
                    #     source_loop_eidx_pool.append(loop_eidx_boot[j])
                    #     source_loop_coords_pool.append(rep_all[j])
                    #     source_loop_key.append(f"({self.n_booted + 1},{j})")
                    #     source_loop_birth_t.append(birth_t[j])
                self.n_booted = self.n_booted + 1
