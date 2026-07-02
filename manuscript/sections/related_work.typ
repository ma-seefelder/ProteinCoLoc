// Related Work — comparison matrix + prose positioning against Tapqir / Costes /
// Manders / Pearson-Spearman / ProteinCoLoc v1 (ADVI), across the seven capability axes
// (amortization, calibrated UQ/SBC, BF/model comparison, registration-UQ, spatial map,
// OOD flag, speed). Includes the mandatory four-axis Tapqir differentiation paragraph.
// Intent: draft the positioning early so downstream experiment phases are shaped by the
// claims they must support; frame honestly (calibration "under the simulator", OOD blind spot).
// ============================================================================
// TAPQIR CELL VERIFICATION (Plan 10-03 Task 1) — the three ASSUMED (⚠) matrix
// cells verified against the primary source before the prose was locked (D-09).
// Source: Ordabayev, Friedman, Gelles, Theobald. "Bayesian machine learning
// analysis of single-molecule fluorescence colocalization images." eLife 2022;
// 11:e73860. DOI 10.7554/eLife.73860 (https://elifesciences.org/articles/73860).
// Full-text methods/model read on 2026-07-02.
//
//   A1 — formal SBC / rank-uniformity calibration?
//     VERDICT: CONFIRMED (no formal SBC). Tapqir is Bayesian and validates via
//     posterior-predictive checking (Gelman 2013) and parameter-recovery on
//     simulated data with known ground truth; the article contains no
//     simulation-based-calibration rank-uniformity or coverage test. It does
//     report that p(specific) values "accurately reflect the uncertainty" — a
//     calibration-of-probabilities claim, not a formal SBC proof. Matrix cell
//     kept as "Partial": Bayesian per-spot probabilities, no formal SBC.
//     Source: eLife 73860, Results ("posterior predictive checking") + Methods.
//
//   A2 — explicit Bayes factor / model-comparison output?
//     VERDICT: CONFIRMED (no explicit BF). Tapqir emits per-spot / per-AOI
//     posterior probabilities p(specific) = p(m=1) of target-specific binding,
//     not a dataset-level Bayes factor for a colocalization-vs-null model
//     comparison. The term "Bayes factor" does not appear in the article.
//     Matrix cell kept as "Partial": per-spot classification probability, no
//     explicit BF. Source: eLife 73860, Results (p(specific)) + Methods.
//
//   A3 — registration/drift uncertainty propagated as an inferred latent?
//     VERDICT: CONFIRMED (external preprocessing). Channel mapping between
//     target/binder images and microscope-drift correction are PREPROCESSING
//     steps done by external software (Friedman & Gelles 2015) before Tapqir
//     runs; registration accuracy is treated as a fixed input that limits
//     localization precision, not as an inferred latent whose uncertainty
//     widens the colocalization posterior. (Tapqir does model within-AOI spot
//     xy position via an affine-beta prior — a distinct quantity from
//     cross-channel registration UQ.) Matrix cell kept as "No — external
//     channel mapping + drift correction". Source: eLife 73860, Fig. 1B legend
//     + Methods ("preprocessing ... drift correction"; registration accuracy).
//
// Outcome: all three ASSUMED cells are supported by the paper; none required
// correction. Prose below is framed conservatively per D-09 (overlapping-but-
// distinct data regimes, not blanket superiority).
// ============================================================================

// Compile-safe interim body (Task 2 replaces this with the matrix + four-axis
// prose). Kept citation-free so the skeleton stays green between task commits.
_Draft slot._ Related-work positioning (comparison matrix + Tapqir four-axis
differentiation) is authored in Plan 10-03 Task 2; the three Tapqir cells above
were verified against Ordabayev et al. 2022 (DOI 10.7554/eLife.73860) first.
