# Initial experiment: conformal fate sets on LARRY in vitro

Run 2026-08-12. Protocol settled with Nianping Liu; every design choice below is his call
or a correction he made.

## Setup

**Data.** LARRY in vitro (Weinreb et al. 2020), from the Klein lab release. Day 2 is a single
well of 28,249 cells sampled before the culture was split; days 4 and 6 are two independent
replicate wells.

**Split (scDiffEq, Vinyard et al. 2025, verbatim).** Fit on well-1 clones, evaluate on well-2
clones, day 2 shared, drop any clone whose day-4/6 members span both wells. Of 5,864 clone
barcodes, 1,408 have both a day-2 member and day-4/6 members confined to one well: 642 in
well 1 (895 day-2 cells) and 766 in well 2 (1,107 day-2 cells). 408 clones are dropped for
spanning both wells.

**Calibration.** Well-2 clones are split randomly at the clone level into calibration and test,
because calibration and test have to be *exchangeable*, not maximally separated. The clone is
the exchangeable unit, so sibling cells never straddle the split.

**Label.** Sequencing kills the cell, so no cell's own fate is ever observed, only its clonal
sisters'. The label is `F_obs`, the distribution over the 11 annotations among the clone's
day-4/6 members. A test cell's realized label is the fate of a uniformly drawn sister, so
per-clone coverage is computed exactly as `sum over fates in C of F_obs(fate)`.

**k\*.** `k*(clone)` = the smallest number of fates whose `F_obs` mass reaches `1 - alpha`.
`k* = 1` means fate was settled, `k* >= 2` means it was not. This replaces any hand-drawn
branch region: it is a property of the label, so nobody has to draw a line. It also turns out
to be the natural floor for set size (the irreducible ambiguity), so `|C| - k*` is the excess.

**Models.** kNN fate propagation (k=25 over 50 PCs of 2,000 HVGs, the field-standard smoothed
clonal fate bias) and multinomial logistic regression on the same features, both trained on
well-1 clones only with `F_obs` as a soft label.

**Set constructions**, all at `alpha = 0.1`, 40 random calibration/test splits:
- `posterior`: smallest set whose *predicted* mass reaches `1 - alpha` (trust the model).
- `conformal-THR`: split conformal, score `1 - p_y`.
- `conformal-APS`: split conformal, adaptive prediction sets (randomized).

All intervals below are bootstrapped over **clones**, not over repeats.

## Result 1: the marginal guarantee holds, and it is cheaper than trusting the model

Multinomial logistic, clones with >= 3 sisters (n = 466 evaluation clones):

| method | marginal coverage | mean set size |
|---|---|---|
| posterior | 0.952 [0.939, 0.964] | 3.21 |
| conformal-THR | 0.907 [0.889, 0.924] | 2.57 |
| conformal-APS | 0.910 [0.894, 0.925] | 2.72 |

The naive posterior sets do **not** undercover marginally. They overcover, by being 25% larger
than they need to be. The original framing of this project ("the Bayesian posteriors the field
relies on undercover") is not what the data says at the margin.

## Result 2: coverage decays monotonically in k\*, and conformal is not exempt

Same setting, conditional on the clone's own ambiguity:

| k* | n clones | conformal-THR | posterior |
|---|---|---|---|
| 1 | 213 | 0.966 [0.945, 0.984] | 0.986 [0.971, 0.998] |
| 2 | 214 | **0.858 [0.828, 0.885]** | 0.928 [0.904, 0.948] |
| 3 | 39 | **0.856 [0.807, 0.902]** | 0.904 [0.863, 0.943] |

The same monotone decay appears with the kNN predictor (0.952 -> 0.866 -> 0.857), so it is not
an artifact of one model. Split conformal buys its marginal 0.90 by overcovering the clones
whose fate was already settled and underspending on the clones that were still undecided. The
miscoverage lands exactly on the biologically interesting cells.

Part of this is structural: with a marginal guarantee and a predictor that cannot tell which
clones are ambiguous, coverage *must* decay in k\*. That is what makes Result 3 the load-bearing
one.

## Result 3: the ambiguity is nearly invisible in the day-2 state, so the standard fix fails

Predicting `k* >= 2` (unsettled) from the day-2 transcriptome, trained on well-1 clones and
evaluated on well-2: **AUC = 0.594** (441 train cells, 668 evaluation cells).

Mondrian (group-wise) conformal calibrated on that predicted group is the textbook repair for
conditional undercoverage. It does essentially nothing here:

| | k*=1 | k*=2 | k*=3 | marginal | size |
|---|---|---|---|---|---|
| pooled | 0.966 | 0.854 | 0.847 | 0.904 | 2.55 |
| Mondrian | 0.964 | 0.867 | 0.849 | 0.910 | 2.68 |

This is the same wall Weinreb hit: the day-2 transcriptome does not carry the fate bias. The
undercoverage is therefore not a calibration bug that a better grouping fixes. It is the assay
failing to tell you where its own ambiguity lives.

## Result 4: per-fate, coverage collapses on the rare populations

Multinomial logistic, >= 3 sisters, conformal-THR, by the clone's majority fate:

| fate | n clones | coverage |
|---|---|---|
| Undifferentiated | 218 | 0.930 [0.913, 0.945] |
| Monocyte | 98 | 0.940 [0.902, 0.973] |
| Neutrophil | 80 | 0.942 [0.909, 0.970] |
| Baso | 46 | 0.840 [0.760, 0.912] |
| Mast | 11 | 0.775 [0.566, 0.942] |
| Meg | 7 | 0.319 [0.179, 0.460] |
| Eos | 3 | 0.415 [0.200, 0.545] |

Nianping predicted this before any of it was run. Mono and Neu are fine; Ery, Meg, Baso, Eos are
not. The bottom three rows have too few clones to carry weight on their own, but the direction
is unambiguous and it matches the prediction.

## Open issues

1. **k\* is confounded with clone size.** A clone with one observed sister has `k* = 1` by
   construction. Restricting to >= 3 sisters is the primary analysis for that reason;
   >= 1, >= 5 and >= 10 are in `results/stats.txt`. The k\*=1 fraction moves from 64% (>= 1
   sister) to 36% (>= 5 sisters), so the threshold matters.
2. **Small strata.** k\*=3 holds 39 clones, and several fates hold under 10. Intervals are
   bootstrapped over clones, but they are still thin.
3. **veloVI / cell2fate / LatentVelo are not in here.** They need spliced and unspliced layers,
   which this data release does not ship. Either we pull a velocyto-processed LARRY or the
   baseline changes.

## Reproducing

```
python scripts/01_parse_day2.py   # streamed from the Klein lab release, see script header
python scripts/02_build_clones.py
python scripts/03_experiment.py
python scripts/04_figures.py
python scripts/05_stats.py
python scripts/06_predict_kstar.py
```
