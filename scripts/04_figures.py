"""Figures: coverage against k*, set size against k*, and the per-fate breakdown."""
import numpy as np
import pandas as pd
from plotnine import *

ROOT = "/Users/krishnaharish/conformal-cell-fate"
ALPHA = 0.1
df = pd.read_csv(f"{ROOT}/results/per_clone.csv.gz")

THEME = theme_bw() + theme(
    figure_size=(7.2, 3.6), panel_grid_minor=element_blank(),
    strip_background=element_rect(fill="#EBEBEB", color="#4D4D4D"),
    legend_position="bottom", legend_title=element_blank(),
)
ORDER = ["posterior", "conformal-THR", "conformal-APS"]
COLS = {"posterior": "#D55E00", "conformal-THR": "#0072B2", "conformal-APS": "#009E73"}


def ci(x):
    """mean and a normal-approx 95% CI over repeats."""
    return pd.Series(dict(m=np.mean(x), lo=np.mean(x) - 1.96 * np.std(x) / np.sqrt(len(x)),
                          hi=np.mean(x) + 1.96 * np.std(x) / np.sqrt(len(x))))


def clone_ci(g, value="coverage", B=4000, seed=0):
    """Bootstrap over CLONES, not repeats: a k* stratum can hold very few clones."""
    v = g.groupby("clone")[value].mean().values
    rng = np.random.default_rng(seed)
    bs = v[rng.integers(0, len(v), size=(B, len(v)))].mean(1)
    return pd.Series(dict(m=v.mean(), lo=np.quantile(bs, 0.025),
                          hi=np.quantile(bs, 0.975), n=len(v)))


# ---- fig 1: coverage vs k*, both models, clones with >= 3 sisters ---------
def coverage_vs_kstar(ms=3, out="fig1_coverage_vs_kstar"):
    d = df[df.min_sisters == ms]
    agg = (d.groupby(["model", "method", "kstar"])
             .apply(clone_ci, include_groups=False).reset_index())
    agg["method"] = pd.Categorical(agg.method, ORDER)
    agg["model"] = agg.model.map({"knn": "kNN fate propagation",
                                  "logreg": "multinomial logistic"})
    p = (ggplot(agg, aes("factor(kstar)", "m", color="method", group="method"))
         + geom_hline(yintercept=1 - ALPHA, linetype="dashed", color="#4D4D4D", size=0.4)
         + geom_line(size=0.7) + geom_point(size=2.4)
         + geom_errorbar(aes(ymin="lo", ymax="hi"), width=0.12, size=0.5)
         + facet_wrap("model")
         + scale_color_manual(values=COLS)
         + labs(x="k*  (smallest number of fates whose observed clonal mass reaches 1 - alpha)",
                y=f"coverage of a random clonal sister\n(nominal {1-ALPHA:.2f}, dashed)")
         + THEME + theme(figure_size=(8.0, 3.8)))
    p.save(f"{ROOT}/figures/{out}.png", dpi=200, verbose=False)
    return agg


# ---- fig 2: set size and excess over k* ----------------------------------
def size_vs_kstar(model="logreg", ms=3, out="fig2_size_vs_kstar"):
    d = df[(df.model == model) & (df.min_sisters == ms)]
    per_rep = d.groupby(["method", "kstar", "rep"])[["size", "excess"]].mean().reset_index()
    a = per_rep.groupby(["method", "kstar"])["size"].apply(ci).unstack().reset_index()
    a["method"] = pd.Categorical(a.method, ORDER)
    floor = pd.DataFrame(dict(kstar=sorted(d.kstar.unique())))
    p = (ggplot(a, aes("factor(kstar)", "m", color="method", group="method"))
         + geom_line(floor, aes("factor(kstar)", "kstar", group=1), color="#4D4D4D",
                     linetype="dashed", size=0.4, inherit_aes=False)
         + geom_line(size=0.7) + geom_point(size=2.4)
         + geom_errorbar(aes(ymin="lo", ymax="hi"), width=0.12, size=0.5)
         + scale_color_manual(values=COLS)
         + labs(x="k*", y="mean prediction-set size\n(dashed = k*, the irreducible floor)")
         + THEME)
    p.save(f"{ROOT}/figures/{out}.png", dpi=200, verbose=False)
    return a


# ---- fig 3: per-fate coverage -------------------------------------------
def per_fate(model="logreg", ms=3, out="fig3_per_fate"):
    d = df[(df.model == model) & (df.min_sisters == ms)]
    a = (d.groupby(["method", "top_fate"])
           .apply(clone_ci, include_groups=False).reset_index())
    n = d[d.method == "posterior"].groupby("top_fate").clone.nunique()
    a = a[a.top_fate.map(n).fillna(0) >= 3]
    a["lab"] = a.top_fate.map(lambda f: f"{f}\n(n={n.get(f,0)})")
    a["method"] = pd.Categorical(a.method, ORDER)
    a = a.sort_values("top_fate", key=lambda s: -s.map(n))
    a["lab"] = pd.Categorical(a.lab, a.lab.unique())
    p = (ggplot(a, aes("lab", "m", color="method"))
         + geom_hline(yintercept=1 - ALPHA, linetype="dashed", color="#4D4D4D", size=0.4)
         + geom_point(size=2.4, position=position_dodge(width=0.5))
         + geom_errorbar(aes(ymin="lo", ymax="hi"), width=0.2, size=0.5,
                         position=position_dodge(width=0.5))
         + scale_color_manual(values=COLS)
         + labs(x="clone's majority fate", y=f"coverage (nominal {1-ALPHA:.2f})")
         + THEME + theme(figure_size=(8.4, 3.8),
                         axis_text_x=element_text(size=7)))
    p.save(f"{ROOT}/figures/{out}.png", dpi=200, verbose=False)
    return a


if __name__ == "__main__":
    print("== coverage vs k* (both models, >=3 sisters, clone bootstrap) ==")
    print(coverage_vs_kstar().to_string(index=False))
    print("\n== set size vs k* ==")
    print(size_vs_kstar().to_string(index=False))
    print("\n== per fate ==")
    print(per_fate().to_string(index=False))
