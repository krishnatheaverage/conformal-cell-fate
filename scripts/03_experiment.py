"""Conformal fate sets on LARRY in vitro: the initial experiment.

Protocol as settled with Nianping Liu:
  fit on well-1 clones, calibrate and test on well-2 clones split at the clone
  level, label = F_obs (a distribution over 11 fates, read off the clonal
  sisters), and condition coverage on k*(clone) = the smallest number of fates
  whose F_obs mass reaches 1 - alpha, instead of on any hand-drawn branch region.

Two set constructions are compared at the same nominal level:
  * "posterior"  : smallest set whose PREDICTED mass reaches 1 - alpha. This is
                   what you do if you trust the model's own probabilities, and it
                   is the thing we claim undercovers.
  * "conformal"  : split conformal, calibrated on held-out well-2 clones.
                   Scores: THR (1 - p_y) and APS (cumulative, randomized).

The exchangeable unit is the clone. Sequencing kills the cell, so a day-2 cell's
realized label is the fate of a uniformly drawn clonal sister; each calibration
clone therefore contributes exactly one score.
"""
import json
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import NearestNeighbors

ROOT = "/Users/krishnaharish/conformal-cell-fate"
FATES = ["Undifferentiated", "Neutrophil", "Monocyte", "Baso", "Mast", "Meg",
         "Lymphoid", "Erythroid", "Eos", "Ccr7_DC", "pDC"]
ALPHA = 0.1
N_REPEATS = 40
N_HVG = 2000
N_PC = 50
K_NN = 25
RNG0 = 20260812


# ---------------------------------------------------------------- features
def build_features():
    X = sparse.load_npz(f"{ROOT}/data/day2_counts.npz").tocsc()
    X.data = np.log1p(X.data)
    n = X.shape[0]
    mean = np.asarray(X.mean(0)).ravel()
    sq = np.asarray(X.multiply(X).mean(0)).ravel()
    var = sq - mean ** 2
    hvg = np.argsort(var)[::-1][:N_HVG]
    Xh = np.asarray(X[:, hvg].todense())
    Xh = (Xh - Xh.mean(0)) / (Xh.std(0) + 1e-8)
    pcs = PCA(n_components=N_PC, random_state=0).fit_transform(Xh)
    pcs = pcs / pcs[:, 0].std()
    print(f"features: {n} day-2 cells, {N_HVG} HVG -> {N_PC} PCs")
    return pcs


# ---------------------------------------------------------------- labels
def load_clones():
    """Clone table, with day2_cells remapped from global cell ids to rows of the
    day-2 expression matrix (which is stored in sorted global-index order)."""
    cl = pd.read_json(f"{ROOT}/data/clones.json")
    F = np.load(f"{ROOT}/data/clone_F.npy")
    day2_global = np.loadtxt(f"{ROOT}/data/day2_rows_1based.txt", dtype=np.int64) - 1
    day2_global.sort()
    pos = {g: i for i, g in enumerate(day2_global)}
    cl["day2_cells"] = cl.day2_cells.map(lambda cs: [pos[c] for c in cs])
    return cl, F


def kstar_from_F(F, alpha):
    s = np.sort(F, axis=1)[:, ::-1]
    return (np.cumsum(s, axis=1) < (1 - alpha) - 1e-12).sum(1) + 1


# ---------------------------------------------------------------- models
def fit_predictor(pcs, cl, F, train_rows, kind):
    """Return a function mapping day-2 cell indices to predicted fate probs."""
    cells, labels = [], []
    for i in train_rows:
        for c in cl.day2_cells.iloc[i]:
            cells.append(c)
            labels.append(F[i])
    cells = np.array(cells)
    Y = np.array(labels)

    if kind == "knn":
        nn = NearestNeighbors(n_neighbors=min(K_NN, len(cells))).fit(pcs[cells])

        def predict(query_cells):
            _, idx = nn.kneighbors(pcs[query_cells])
            P = Y[idx].mean(1)
            return P / P.sum(1, keepdims=True)
    elif kind == "logreg":
        # soft labels -> replicate each cell once per fate, weighted by F_obs
        reps, ys, ws = [], [], []
        for j, y in enumerate(Y):
            for f in np.nonzero(y)[0]:
                reps.append(cells[j])
                ys.append(f)
                ws.append(y[f])
        clf = LogisticRegression(max_iter=3000, C=0.5)
        clf.fit(pcs[np.array(reps)], np.array(ys), sample_weight=np.array(ws))
        classes = clf.classes_

        def predict(query_cells):
            P = np.zeros((len(query_cells), len(FATES)))
            P[:, classes] = clf.predict_proba(pcs[query_cells])
            return P / P.sum(1, keepdims=True)
    return predict


# ---------------------------------------------------------------- sets
def posterior_sets(P, alpha):
    """Smallest set whose PREDICTED mass reaches 1 - alpha (trust the model)."""
    order = np.argsort(P, axis=1)[:, ::-1]
    cum = np.cumsum(np.take_along_axis(P, order, 1), axis=1)
    k = (cum < (1 - alpha) - 1e-12).sum(1) + 1
    S = np.zeros_like(P, dtype=bool)
    for i in range(len(P)):
        S[i, order[i, :k[i]]] = True
    return S


def thr_scores(P, y):
    return 1.0 - P[np.arange(len(y)), y]


def aps_scores(P, y, rng):
    order = np.argsort(P, axis=1)[:, ::-1]
    Ps = np.take_along_axis(P, order, 1)
    cum = np.cumsum(Ps, axis=1)
    rank = np.argmax(order == y[:, None], axis=1)
    below = cum[np.arange(len(y)), rank] - Ps[np.arange(len(y)), rank]
    return below + rng.random(len(y)) * Ps[np.arange(len(y)), rank]


def thr_sets(P, qhat):
    return P >= (1.0 - qhat)


def aps_sets(P, qhat, rng):
    order = np.argsort(P, axis=1)[:, ::-1]
    Ps = np.take_along_axis(P, order, 1)
    cum = np.cumsum(Ps, axis=1)
    u = rng.random((len(P), 1))
    keep = (cum - Ps * u) <= qhat
    keep[:, 0] = True
    S = np.zeros_like(P, dtype=bool)
    np.put_along_axis(S, order, keep, axis=1)
    return S


def conformal_q(scores, alpha):
    n = len(scores)
    lvl = min(1.0, np.ceil((n + 1) * (1 - alpha)) / n)
    return np.quantile(scores, lvl, method="higher")


# ---------------------------------------------------------------- experiment
def run(pcs, cl, F, min_sisters, model_kind, alpha=ALPHA):
    ks = kstar_from_F(F, alpha)
    keep = cl.n_late.values >= min_sisters
    tr = np.where((cl.well.values == 1) & keep)[0]
    ev = np.where((cl.well.values == 2) & keep)[0]
    predict = fit_predictor(pcs, cl, F, tr, model_kind)

    # one representative day-2 cell per evaluation clone (drawn per repeat)
    recs = []
    for rep in range(N_REPEATS):
        rng = np.random.default_rng(RNG0 + rep)
        perm = rng.permutation(len(ev))
        half = len(ev) // 2
        cal_i, test_i = ev[perm[:half]], ev[perm[half:]]

        def rep_cells(idx):
            return np.array([cl.day2_cells.iloc[i][rng.integers(len(cl.day2_cells.iloc[i]))]
                             for i in idx])

        Pcal = predict(rep_cells(cal_i))
        ycal = np.array([rng.choice(len(FATES), p=F[i]) for i in cal_i])
        Ptest = predict(rep_cells(test_i))
        Ftest = F[test_i]

        sets = {}
        sets["posterior"] = posterior_sets(Ptest, alpha)
        sets["conformal-THR"] = thr_sets(Ptest, conformal_q(thr_scores(Pcal, ycal), alpha))
        sets["conformal-APS"] = aps_sets(Ptest, conformal_q(aps_scores(Pcal, ycal, rng), alpha), rng)

        for name, S in sets.items():
            cov = (S * Ftest).sum(1)          # exact E[cover] over the sister draw
            size = S.sum(1)
            for j, ci in enumerate(test_i):
                recs.append(dict(rep=rep, method=name, clone=int(cl.clone.iloc[ci]),
                                 kstar=int(ks[ci]), coverage=float(cov[j]), size=int(size[j]),
                                 excess=int(size[j]) - int(ks[ci]),
                                 n_late=int(cl.n_late.iloc[ci]),
                                 top_fate=cl.top_fate.iloc[ci]))
    df = pd.DataFrame(recs)
    df["min_sisters"] = min_sisters
    df["model"] = model_kind
    return df, len(tr), len(ev)


if __name__ == "__main__":
    pcs = build_features()
    cl, F = load_clones()
    out = []
    for model_kind in ("knn", "logreg"):
        for ms in (1, 3, 5, 10):
            df, ntr, nev = run(pcs, cl, F, ms, model_kind)
            out.append(df)
            m = df.groupby("method")["coverage"].mean()
            print(f"[{model_kind} >= {ms} sisters] train clones {ntr}, eval clones {nev} | "
                  + " ".join(f"{k}={v:.3f}" for k, v in m.items()))
    allo = pd.concat(out)
    allo.to_csv(f"{ROOT}/results/per_clone.csv.gz", index=False)
    print("\nwrote results/per_clone.csv.gz", allo.shape)
