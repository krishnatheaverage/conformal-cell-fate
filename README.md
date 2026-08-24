# Conformal Cell Fate

**Distribution-free uncertainty for vector-field models of development and in-silico perturbation.**

A short research proposal. The goal is to give the vector-field cell-fate machinery (Dynamo, Spateo,
and the wider Aristotle ecosystem) something it currently lacks: a prediction that comes with a
finite-sample, distribution-free guarantee, instead of a point estimate or a model-internal
posterior that is known to be miscalibrated.

Authors: Krishna Harish and Nianping Liu (Stanford). Target: ICBINB-BIO at
NeurIPS 2026.

Phase-0 has been run. The headline below is the original proposal; what the data actually said is
in [RESULTS.md](RESULTS.md), and it is not what we expected. Read that first.

---

## The gap

Dynamo reconstructs a continuous vector field from single-cell data and then *uses* that field to
predict cell fate (integrate streamlines), find least-action reprogramming paths, and predict the
outcome of in-silico gene perturbations. Spateo does the analogous thing in space and time. These
are the load-bearing outputs of the whole pipeline, and they are delivered as point predictions.

The existing uncertainty work on this problem (cell2fate, veloVI, LatentVelo) is entirely Bayesian:
you get a posterior spread and a heuristic "confidence score." Two problems with that:

1. **No coverage guarantee.** A posterior interval is only as calibrated as the model is correct.
   In practice these posteriors are overconfident, and worst exactly where it matters (branch
   points, transition states, out-of-distribution cells / species).
2. **Not comparable across methods.** A confidence score from one velocity model does not mean the
   same thing as one from another, so you cannot use it to rank or trust predictions method-agnostically.

There is no method that says, with a finite-sample guarantee: *"this predicted fate (or this
predicted perturbation effect) is correct with probability at least 1 - alpha."* That primitive is
missing, and conformal prediction supplies it directly.

## The proposal

Wrap the vector-field cell-fate machinery in split-conformal prediction. Three concrete instruments:

**1. Conformal fate sets.** For a starting cell, output a *set* of candidate terminal states with
guaranteed marginal coverage 1 - alpha, calibrated on held-out cells with known fate (lineage-traced
or terminal-annotated). Nonconformity score = distance between the vector-field-integrated endpoint
and the true terminal state, or the negative predicted-fate probability. Small set = confident cell;
large set = the field genuinely does not resolve this cell's future, and now you know it.

**2. Conformal perturbation intervals.** For an in-silico perturbation, output a calibrated set /
interval on the predicted shift in cell state or fate probability. This lets you rank perturbations
by *reliable* effect size, and only trigger a wet-lab (e.g. CRISPR) follow-up when the predicted
effect is both large and confidently non-null. The payoff is direct: fewer wasted experiments.

**3. Coverage as a diagnostic, including across species.** Where does coverage break? Localizing
miscoverage tells you exactly where the vector field is untrustworthy. Pushed to evo-devo: calibrate
on species A, measure coverage on species B. The conformal drift is then a *quantitative* score for
how conserved a developmental vector field is, i.e. developmental constraint you can put a number on.

## Why it is novel

Conformal prediction has touched some genomics classification tasks, but not vector-field /
dynamical cell-fate models, not in-silico perturbation ranking, and not as a cross-species
conservation diagnostic. The uncertainty literature in this exact area is Bayesian and
known-miscalibrated. Distribution-free coverage is a clean missing piece, and it is a natural,
lightweight add-on to tools the field already runs.

## Phase-0 (the first self-contained result)

Everything below is buildable on public data with the lab's own tools.

1. Take a standard velocity benchmark with known terminal states (e.g. pancreatic endocrinogenesis,
   hematopoiesis). Fit the vector field with Dynamo.
2. Split-conformal calibration of fate sets. Show empirical coverage lands at 1 - alpha across
   alpha, while the Bayesian posterior "confidence" undercovers (report the miscoverage gap).
   **This step was run and the second half is refuted: the posteriors overcover (0.952 marginal),
   they do not undercover. See RESULTS.md.**
3. Show set size tracks biology: sets blow up at branch points and shrink along committed lineages.
4. Stretch: one in-silico perturbation ranked by conformal effect vs. by point effect, and one
   cross-species coverage number.

Step 2 alone is a clean, honest, self-contained result and the backbone of a workshop paper.

## Status

Phase-0 complete on LARRY in vitro (Weinreb 2020), scDiffEq split. 1,408 clones after the one-well
filter. The result is a negative one and it is not the one we proposed:

- The Bayesian posteriors do not undercover. They overcover, at 0.952 marginal, with sets about 25%
  larger than needed. Split conformal reaches 0.907 with smaller sets.
- What does fail is conditional: coverage decays monotonically in k*, the smallest number of fates
  whose observed clonal mass reaches 1 - alpha. At >=3 sisters, k*=1 gets 0.966 and k*=2 gets 0.858.
  The same decay shows up under a kNN fate-propagation predictor, so it is not an artifact of one
  model.
- It cannot be repaired by grouping. Predicting k*>=2 from the day-2 transcriptome runs at AUC
  0.594, and Mondrian conformal on the predicted group moves k*=2 coverage only 0.854 -> 0.867.

Numbers, tables, and caveats in [RESULTS.md](RESULTS.md). Scripts run in order, `scripts/01` to
`scripts/06`.
