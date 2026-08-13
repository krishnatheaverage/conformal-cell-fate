"""Parse the filtered (day-2 only) MatrixMarket stream from stdin into a CSR npz.

Input lines are "row col val" with 1-based indices from the LARRY in vitro
normed-counts matrix (130887 cells x 25289 genes), already filtered by awk to the
28,249 day-2 rows. We only ever need day-2 expression: fate labels come from the
clone matrix plus the day-4/6 annotations, not from day-4/6 expression.
"""
import sys
import numpy as np
import pandas as pd
from scipy import sparse

DAY2_ROWS = "/Users/krishnaharish/conformal-cell-fate/data/day2_rows_1based.txt"
OUT = "/Users/krishnaharish/conformal-cell-fate/data/day2_counts.npz"
N_GENES = 25289

keep = np.loadtxt(DAY2_ROWS, dtype=np.int64)
keep.sort()
remap = -np.ones(keep.max() + 2, dtype=np.int64)
remap[keep] = np.arange(len(keep))

chunks = []
nnz = 0
reader = pd.read_csv(
    sys.stdin, sep=" ", header=None, names=["r", "c", "v"],
    dtype={"r": np.int64, "c": np.int32, "v": np.float32},
    chunksize=8_000_000, engine="c",
)
for i, ch in enumerate(reader):
    r = remap[ch["r"].values]
    assert (r >= 0).all(), "unexpected row id outside the day-2 set"
    chunks.append((r.astype(np.int32), ch["c"].values - 1, ch["v"].values))
    nnz += len(ch)
    print(f"  chunk {i}: {nnz:,} nnz", file=sys.stderr, flush=True)

rows = np.concatenate([c[0] for c in chunks])
cols = np.concatenate([c[1] for c in chunks])
vals = np.concatenate([c[2] for c in chunks])
del chunks

X = sparse.coo_matrix((vals, (rows, cols)), shape=(len(keep), N_GENES)).tocsr()
sparse.save_npz(OUT, X)
print(f"saved {OUT} shape={X.shape} nnz={X.nnz:,}", file=sys.stderr)
