"""Build the clone-level fate labels and the scDiffEq well split.

Design (settled with Nianping Liu, 2026-08-12):

  * scDiffEq split verbatim (Vinyard et al., Nat Mach Intell 2025): day 2 is a
    single well of 28,249 cells sampled before the culture was split; days 4 and 6
    are two independent replicate wells. Fit on well 1, evaluate on well 2, day 2
    shared, and drop any clone whose day-4/6 members span both wells.
  * You never observe a cell's own fate. Sequencing kills the cell, so the label
    for a day-2 cell is what its clonal sisters became at day 4/6. The label object
    is therefore F_obs, a distribution over the 11 annotations, not a label.
  * k*(clone) = the smallest number of fates whose F_obs mass reaches 1 - alpha.
    k* = 1 means fate was already settled for that clone, k* >= 2 means it was not.
    This replaces any hand-drawn "branch region": it is a property of the label,
    so neither of us has to draw a line.
"""
import numpy as np
import pandas as pd
from scipy import sparse, io
import gzip
import json

ROOT = "/Users/krishnaharish/conformal-cell-fate"

FATES = ["Undifferentiated", "Neutrophil", "Monocyte", "Baso", "Mast", "Meg",
         "Lymphoid", "Erythroid", "Eos", "Ccr7_DC", "pDC"]

md = pd.read_csv(f"{ROOT}/data/metadata.csv.gz")
with gzip.open(f"{ROOT}/data/stateFate_inVitro_clone_matrix.mtx.gz", "rb") as fh:
    C = io.mmread(fh).tocsc()  # cells x clones
print("clone matrix:", C.shape, "nnz:", C.nnz)

n_cells, n_clones = C.shape
assert n_cells == len(md)

fate_idx = {f: i for i, f in enumerate(FATES)}
ann = md["Cell type annotation"].map(fate_idx).values
tp = md["Time point"].values
well = md["Well"].values

Ccsr = C.tocsr()
cell_clone = -np.ones(n_cells, dtype=np.int64)
multi = 0
for i in range(n_cells):
    lo, hi = Ccsr.indptr[i], Ccsr.indptr[i + 1]
    if hi - lo == 1:
        cell_clone[i] = Ccsr.indices[lo]
    elif hi - lo > 1:
        multi += 1
print(f"cells with a clone barcode: {(cell_clone >= 0).sum():,} "
      f"(cells assigned to >1 clone, dropped: {multi})")

rows = []
for cl in range(n_clones):
    members = C.indices[C.indptr[cl]:C.indptr[cl + 1]]
    members = members[cell_clone[members] == cl]
    if len(members) == 0:
        continue
    d2 = members[tp[members] == 2.0]
    late = members[tp[members] > 2.0]
    if len(d2) == 0 or len(late) == 0:
        continue
    wells_late = set(np.unique(well[late]).tolist())
    if len(wells_late) != 1:
        continue  # clone spans both replicate wells -> dropped, per scDiffEq
    w = wells_late.pop()
    counts = np.bincount(ann[late], minlength=len(FATES)).astype(float)
    rows.append(dict(clone=cl, well=w, n_day2=len(d2), n_late=len(late),
                     day2_cells=d2.tolist(), counts=counts.tolist()))

cl_df = pd.DataFrame(rows)
counts = np.array(cl_df["counts"].tolist())
F = counts / counts.sum(1, keepdims=True)


def kstar(F, alpha=0.1):
    """Smallest number of fates whose F_obs mass reaches 1 - alpha."""
    s = np.sort(F, axis=1)[:, ::-1]
    return (np.cumsum(s, axis=1) < (1 - alpha) - 1e-12).sum(1) + 1


for a in (0.05, 0.1, 0.2):
    cl_df[f"kstar_a{a}"] = kstar(F, a)

cl_df["n_fates_observed"] = (counts > 0).sum(1)
cl_df["top_fate"] = [FATES[i] for i in F.argmax(1)]
np.save(f"{ROOT}/data/clone_F.npy", F)
cl_df.drop(columns=["counts"]).to_json(f"{ROOT}/data/clones.json", orient="records")

print("\n== clones kept (day-2 member + late members in exactly one well) ==")
print(cl_df.groupby("well").agg(clones=("clone", "size"),
                                day2_cells=("n_day2", "sum"),
                                late_cells=("n_late", "sum")))
print("\n== k* distribution at alpha=0.1 (well 2 = the evaluation well) ==")
print(pd.crosstab(cl_df["kstar_a0.1"], cl_df["well"]))
print("\n== clone size (late members) ==")
print(cl_df["n_late"].describe())
print("\n== top fate of each clone, well 2 ==")
print(cl_df[cl_df.well == 2]["top_fate"].value_counts())
print("\n== k*=1 clones by which fate is the settled one (well 2) ==")
w2 = cl_df[(cl_df.well == 2) & (cl_df["kstar_a0.1"] == 1)]
print(w2["top_fate"].value_counts())
