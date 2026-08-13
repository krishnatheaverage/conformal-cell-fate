"""Is k* predictable from the day-2 state, and does stratifying on it fix coverage?

k* is a property of the LABEL, so conditioning coverage on it is an oracle
stratification: at deployment you do not know a new cell's clonal fate
distribution. If the k*-conditional undercoverage is to be fixed rather than only
reported, k* has to be predictable from the day-2 transcriptome itself. This
script tests that, then runs Mondrian (group-wise) conformal on the PREDICTED
group and asks whether conditional coverage is restored.
"""
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

import sys
sys.path.insert(0, "/Users/krishnaharish/conformal-cell-fate/scripts")
from importlib import import_module
E = import_module("03_experiment")

ROOT = "/Users/krishnaharish/conformal-cell-fate"
ALPHA = 0.1
MIN_SIS = 3
N_REPEATS = 40


def main():
    pcs = E.build_features()
    cl, F = E.load_clones()
    ks = E.kstar_from_F(F, ALPHA)
    keep = cl.n_late.values >= MIN_SIS
    tr = np.where((cl.well.values == 1) & keep)[0]
    ev = np.where((cl.well.values == 2) & keep)[0]

    # ---- can we predict "unsettled" (k* >= 2) from the day-2 state? -------
    def cells_and_y(idx):
        c, y = [], []
        for i in idx:
            for cell in cl.day2_cells.iloc[i]:
                c.append(cell)
                y.append(int(ks[i] >= 2))
        return np.array(c), np.array(y)

    ctr, ytr = cells_and_y(tr)
    cev, yev = cells_and_y(ev)
    clf = LogisticRegression(max_iter=3000, C=0.5).fit(pcs[ctr], ytr)
    p_unsettled = clf.predict_proba(pcs[cev])[:, 1]
    auc = roc_auc_score(yev, p_unsettled)
    print(f"predicting k* >= 2 from the day-2 transcriptome:")
    print(f"  train {len(ctr)} day-2 cells ({ytr.mean():.1%} unsettled), "
          f"eval {len(cev)} ({yev.mean():.1%} unsettled)")
    print(f"  AUC = {auc:.3f}   (0.5 = the day-2 state carries no signal about settledness)")

    predict = E.fit_predictor(pcs, cl, F, tr, "logreg")
    grp_model = clf

    # ---- Mondrian conformal on the PREDICTED group -----------------------
    recs = []
    for rep in range(N_REPEATS):
        rng = np.random.default_rng(E.RNG0 + rep)
        perm = rng.permutation(len(ev))
        cal_i, test_i = ev[perm[:len(ev) // 2]], ev[perm[len(ev) // 2:]]

        def rep_cells(idx):
            return np.array([cl.day2_cells.iloc[i][rng.integers(len(cl.day2_cells.iloc[i]))]
                             for i in idx])

        ccal, ctest = rep_cells(cal_i), rep_cells(test_i)
        Pcal, Ptest = predict(ccal), predict(ctest)
        ycal = np.array([rng.choice(len(E.FATES), p=F[i]) for i in cal_i])
        gcal = (grp_model.predict_proba(pcs[ccal])[:, 1] >= 0.5).astype(int)
        gtest = (grp_model.predict_proba(pcs[ctest])[:, 1] >= 0.5).astype(int)
        scal = E.thr_scores(Pcal, ycal)

        q_pooled = E.conformal_q(scal, ALPHA)
        S_pooled = E.thr_sets(Ptest, q_pooled)
        S_mond = np.zeros_like(S_pooled)
        for g in (0, 1):
            m_cal, m_te = gcal == g, gtest == g
            q = E.conformal_q(scal[m_cal], ALPHA) if m_cal.sum() >= 10 else q_pooled
            S_mond[m_te] = E.thr_sets(Ptest[m_te], q)

        for name, S in (("pooled", S_pooled), ("mondrian", S_mond)):
            cov = (S * F[test_i]).sum(1)
            for j, ci_ in enumerate(test_i):
                recs.append(dict(rep=rep, method=name, clone=int(cl.clone.iloc[ci_]),
                                 kstar=int(ks[ci_]), coverage=float(cov[j]),
                                 size=int(S[j].sum())))
    d = pd.DataFrame(recs)
    d.to_csv(f"{ROOT}/results/mondrian.csv.gz", index=False)

    def boot(g, v="coverage", B=4000):
        x = g.groupby("clone")[v].mean().values
        r = np.random.default_rng(0)
        bs = x[r.integers(0, len(x), size=(B, len(x)))].mean(1)
        return f"{x.mean():.3f} [{np.quantile(bs,.025):.3f}, {np.quantile(bs,.975):.3f}]"

    print(f"\nsplit conformal (THR), pooled vs Mondrian on predicted settledness:")
    print(f"{'':22s} {'k*=1':>24s} {'k*=2':>24s} {'k*=3':>24s} {'marginal':>24s}   size")
    for name, g in d.groupby("method"):
        cells = [boot(g[g.kstar == k]) for k in (1, 2, 3)]
        print(f"  {name:20s} " + " ".join(f"{c:>24s}" for c in cells)
              + f" {boot(g):>24s}   {g.groupby('clone')['size'].mean().mean():.2f}")


if __name__ == "__main__":
    main()
