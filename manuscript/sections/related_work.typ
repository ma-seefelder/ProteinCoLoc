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

Quantitative colocalization analysis of two-channel fluorescence microscopy has
historically relied on threshold-dependent co-occurrence coefficients
@manders1993 and on correlation-with-randomization significance testing
@costes2004. These classical estimators return a single scalar or a frequentist
$p$-value and carry no posterior uncertainty, no model-comparison Bayes factor,
and no mechanism to signal when their own assumptions are violated. The original
ProteinCoLoc (v1) improved on this with a hierarchical Bayesian model and a
kernel-density Bayes factor, but its posterior is fit per dataset by ADVI and it
provides neither a calibration proof nor a registration- or misspecification-aware
output. The most closely related modern method is Tapqir @ordabayev2022, a
physics-based Bayesian model of single-molecule colocalization (CoSMoS) image
data. @tab:related-work positions ProteinCoLoc v2.0 against these methods across
seven capability axes; the three Tapqir cells were verified against the primary
source (see the provenance note in this file's source).

// Wide 8-column matrix — placed on its own landscape page (equivalent to lib/helpers.typ
// landscapetable) so the seven capability columns read without cramping, rather than being
// squeezed into the two-column body.
#page(flipped: true, columns: 1)[
  #figure(
    kind: table,
    caption: [Capability comparison of colocalization-analysis methods across
      seven axes. "Yes" = capability present; "No" = absent; "Partial" =
      present in a qualified or domain-restricted form (see cell text). The
      three Tapqir cells marked Partial/No for calibration, Bayes factor, and
      registration were verified against Ordabayev et al. 2022
      (DOI 10.7554/eLife.73860); ProteinCoLoc v2.0 entries flagged with a phase
      number are _reserved slots_ that later phases substantiate, not yet-proven
      results. v2.0 calibration is stated as holding _under the simulator_.],
    text(size: 9pt)[
      #table(
        columns: (1.5fr, 1fr, 1.5fr, 1.5fr, 1.3fr, 1.3fr, 1.1fr, 1.2fr),
        align: (left, left, left, left, left, left, left, left),
        table.header(
          [Method],
          [Amortized inference],
          [Calibrated UQ / SBC],
          [Bayes factor / model comparison],
          [Registration UQ as latent],
          [Spatial per-region map],
          [OOD / misspec. flag],
          [Speed],
        ),
        [ProteinCoLoc v2.0],
        [Yes — single forward pass],
        [Yes — formal SBC, under the simulator (Ph 5)],
        [Yes — 2-way (Ph 5) → 3-way (Ph 13)],
        [Yes — inferred latent (Ph 11)],
        [Yes — per-region $Delta rho$ map (Ph 12)],
        [Yes — summary-bounded (Ph 5)],
        [ms (forward pass)],

        [Tapqir @ordabayev2022],
        [No — per-dataset SVI],
        [Partial — Bayesian per-spot probabilities; posterior-predictive checks, no formal SBC],
        [Partial — per-spot $p$(specific); no explicit model-comparison BF],
        [No — external channel mapping + drift correction],
        [Partial — per-AOI, single-molecule regime],
        [No — no explicit misspecification flag],
        [Slow — per-dataset SVI],

        [Costes randomization @costes2004],
        [No],
        [No — $p$-value only],
        [No],
        [No],
        [No — global],
        [No],
        [Moderate — randomization test],

        [Manders M1/M2 @manders1993],
        [No],
        [No],
        [No],
        [No],
        [No — global],
        [No],
        [Fast],

        [Pearson / Spearman],
        [No],
        [No],
        [No],
        [No],
        [No — global],
        [No],
        [Fast],

        [ProteinCoLoc v1 (ADVI)],
        [No — per-dataset ADVI],
        [Partial — Bayesian posterior, no SBC proof],
        [Yes — KDE $Delta rho$ Bayes factor],
        [No — fixed registration],
        [Partial — patch-level, no spatial prior],
        [No],
        [~0.5 s – min per dataset],
      )
    ],
  ) <tab:related-work>
]

== Differentiation from Tapqir

Tapqir @ordabayev2022 is the most methodologically adjacent prior work, so we
state the delta explicitly along four axes. We frame it honestly: Tapqir and
ProteinCoLoc v2.0 target *overlapping but distinct data regimes* — Tapqir is
specialized to sparse single-molecule CoSMoS spots, whereas v2.0 addresses
general dense/diffuse two-channel patch-correlation microscopy — so this is a
statement of differing capabilities and domains, not of blanket superiority.

*Amortization.* This is the cleanest, least-contestable delta. Tapqir fits a
variational posterior *per dataset* by stochastic variational inference, so each
new dataset pays the full optimization cost. ProteinCoLoc v2.0 trains a neural
posterior/ratio estimator once and then infers each new dataset in a single
amortized forward pass (milliseconds), with the honest caveat that the headline
speedup is a forward-pass figure (median ~325×) that reduces to ~16× on the full
posterior-sampling workload against a realized ~0.5 s/pair ADVI baseline, at
matched RMSE.

*Calibration / SBC.* v2.0 provides a formal simulation-based-calibration
rank-uniformity and coverage proof @talts2018 — but we state this precisely as
holding *under the simulator*: it certifies that the amortized posterior is
calibrated for data drawn from the training forward model, not unconditionally.
Tapqir is itself Bayesian and validates through posterior-predictive checks and
parameter recovery on simulated data, and reports that its per-spot probabilities
reflect classification uncertainty; what it does not report is a formal SBC
rank-uniformity proof (verified against the paper). The defensible claim is
therefore "v2.0 adds a formal SBC coverage proof," not "Tapqir is uncalibrated."

*Registration uncertainty as a latent.* v2.0 (Phase 11) promotes sub-pixel
channel registration to an inferred latent, so the colocalization posterior
widens honestly when registration is uncertain. In Tapqir, cross-channel mapping
and microscope-drift correction are external preprocessing steps (performed by
separate software before inference); registration accuracy bounds localization
precision but is not propagated as an inferred latent into the colocalization
call (verified against the paper). Tapqir does model within-AOI spot position,
which is a distinct quantity from cross-channel registration uncertainty.

*Spatial per-region map.* We frame this as complementary domains rather than
superiority. v2.0 (Phase 12) produces an amortized per-region $Delta rho$ map
with per-region uncertainty over a dense correlation lattice, appropriate to
dense/diffuse fluorescence. Tapqir operates on sparse single-molecule areas of
interest and yields per-spot results in that regime; the two address overlapping
but different spatial problems.

Finally, the honesty posture is built into the positioning itself: v2.0's
calibration guarantee is *under the simulator*, and its out-of-distribution flag
is *summary-orthogonal* — because it is computed from the same fixed
patch-correlation summary statistic that feeds inference, misspecifications that
leave that summary unchanged are provably invisible to it. This blind spot is
named rather than hidden, and must be carried forward into the Discussion and the
OOD figure rather than quietly dropped.
