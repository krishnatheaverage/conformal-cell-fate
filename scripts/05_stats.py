"""Statistics with the clone as the sampling unit.

The per-repeat spread in 04_figures only reflects the calibration/test split and
the sister draw. It does NOT reflect the fact that a k* stratum can contain very
few clones. Every interval reported to a coauthor is bootstrapped over CLONES.
"""
import numpy as np
import pandas as pd

ROOT = "/Users/krishnaharish/conformal-cell-fate"
ALPHA = 0.1
B = 4000
df = pd.read_csv(f"{ROOT}/results/per_clone.csv.gz")


def clone_boot(d, value="coverage", B=B, seed=0):
    """Mean of `value`, bootstrapped over clones (averaging repeats within clone)."""
    per_clone = d.groupby("clone")[value].mean()
    v = per_clone.values
    if len(v) == 0:
        return dict(n=0, m=np.nan, lo=np.nan, hi=np.nan)
    rng = np.random.default_rng(seed)
    bs = v[rng.integers(0, len(v), size=(B, len(v)))].mean(1)
    return dict(n=len(v), m=v.mean(), lo=np.quantile(bs, 0.025), hi=np.quantile(bs, 0.975))


def table(model, ms, by="kstar", value="coverage"):
    d = df[(df.model == model) & (df.min_sisters == ms)]
    rows = []
    for (meth, k), g in d.groupby(["method", by]):
        rows.append(dict(method=meth, **{by: k}, **clone_boot(g, value)))
    return pd.DataFrame(rows).sort_values([by, "method"])


if __name__ == "__main__":
    for model in ("knn", "logreg"):
        for ms in (3, 5):
            print(f"\n{'='*74}\n{model}, clones with >= {ms} sisters, alpha={ALPHA}\n{'='*74}")
            cov = table(model, ms, "kstar", "coverage")
            siz = table(model, ms, "kstar", "size")
            m = cov.merge(siz[["method", "kstar", "m"]], on=["method", "kstar"],
                          suffixes=("", "_size"))
            m["coverage"] = m.apply(lambda r: f"{r.m:.3f} [{r.lo:.3f}, {r.hi:.3f}]", axis=1)
            m["set size"] = m.m_size.round(2)
            m["undercovers?"] = np.where(m.hi < 1 - ALPHA, "YES", "")
            print(m[["kstar", "n", "method", "coverage", "set size", "undercovers?"]]
                  .to_string(index=False))

            print("\n  marginal:")
            d = df[(df.model == model) & (df.min_sisters == ms)]
            for meth, g in d.groupby("method"):
                b = clone_boot(g)
                s = clone_boot(g, "size")
                print(f"    {meth:14s} coverage {b['m']:.3f} [{b['lo']:.3f}, {b['hi']:.3f}]"
                      f"   mean size {s['m']:.2f}   (n={b['n']} clones)")

    print(f"\n{'='*74}\nper-fate coverage (knn, >= 3 sisters)\n{'='*74}")
    pf = table("knn", 3, "top_fate", "coverage")
    pf["coverage"] = pf.apply(lambda r: f"{r.m:.3f} [{r.lo:.3f}, {r.hi:.3f}]", axis=1)
    print(pf[["top_fate", "n", "method", "coverage"]].to_string(index=False))
