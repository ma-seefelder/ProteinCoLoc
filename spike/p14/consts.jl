#=
ProteinCoLoc: A Julia package for the analysis of protein co-localization in microscopy images
Copyright (C) 2023  Dr. rer. nat. Manuel Seefelder
E-Mail: manuel.seefelder@uni-ulm.de
Postal address: Department of Gene Therapy, University of Ulm, Helmholzstr. 8/1, 89081 Ulm, Germany

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
=#

# spike/p14/consts.jl --- Phase-14 Tier-1 pre-registration (D-07).
#
# THE ANTI-SNOOPING CONTRACT. EVERY SEED, COUNTER, SAMPLE SIZE, ALPHA AND BAR THAT ANY
# REPORTED PHASE-14 NUMBER WILL BE SCORED AGAINST IS LOCKED IN THIS ONE FILE AND
# COMMITTED BEFORE ANY PHASE-14 CODE RUNS AND BEFORE A SINGLE DECISION, ABSTENTION OR
# COVERAGE NUMBER EXISTS. Phase 14 turns a posterior and a three-way Bayes factor into
# an ACTIONABLE CALL -- {coloc / not / ABSTAIN} at a user-set FDR -- and the only thing
# that makes its eventual numbers citable is that the bars were set FIRST. Changing a
# value below after a run would be "tune until it passes" data-snooping. These constants
# ARE the committed pre-registration.
#
# TWO-TIER STRUCTURE -- READ THIS BEFORE EDITING ANYTHING.
#   TIER 1 is the single guard block in THIS file, committed before anything runs.
#   TIER 2 is INHERITED, NOT APPENDED HERE. The three-way cut threshold tau is
#     DELIBERATELY ABSENT from this file, and its absence is asserted by
#     spike/test/test_p14_consts.jl. D-07 says tau is INHERITED: Phase 14 LOADS it from
#     the Phase-13 artifact with its provenance asserted, and never re-derives it and
#     never copies its literal here. Copying the value is how two numbers silently
#     diverge; re-deriving guarantees a second, different cut on the same hypothesis
#     space. Phase 13 MEASURED that number (spike/p13/consts.jl, Tier-2 block) and
#     Phase 14 has no licence to have an opinion about it.
#     NOTE THE DIFFERENCE FROM PHASE 13. Phase 13's Tier 2 was a value it would MEASURE
#     and APPEND to its own file. Phase 14's "Tier 2" is a value it will LOAD from
#     someone else's file. There is therefore NO appended second guard block in this
#     file, and there must never be one.
#   NO CONSTANT IN THIS BLOCK IS EVER *EDITED*. Constants are only ever APPENDED, so the
#   git history of this file is itself the pre-registration audit trail. A diff that
#   MODIFIES a line below is, by construction, a pre-registration breach.
#
# SEED DISCIPLINE. P14_DEV_SEED and P14_FIX_SEED are FRESH Random123 streams, EXECUTABLY
# proven disjoint from every reserved stream in the repository. The forbidden inventory is
# NOT retyped: it is DERIVED by including spike/p13/consts.jl and RECOMPUTING
# `_p13_forbidden()` -- which itself recomputes both frozen ship-gate seed families
# (PROD_SEED, PROD_SEED_V2) rather than trusting a comment -- and then adding Phase 13's
# OWN two seeds, which Phase 13 could not forbid to itself. Beyond seed identity this file
# asserts the FULL PHILOX KEY-WORD CROSS-PRODUCT against the entire Phase-13 family:
# two seeds that differ are not enough if `seed xor salt xor counter` can still collide.
#
# DECOUPLING (hard constraint, CLAUDE.md and D-01): spike-local constants; touches no src/,
# adds no dependency, and reaches spike/p13/ only READ-ONLY through a guarded include.
#
# NOTE THE PATH DEPTH: from spike/p14/ the repo root is TWO levels up, not three
# (joinpath(@__DIR__, "..", "..", "src", ...)). Stated here because the same reminder is
# carried in three separate spike/p13/ files and the mistake is silent.
#
# Guarded as ONE Tier-1 block keyed on :P14_DECLARED_DEVIATIONS -- a name this file alone
# declares, declared UNCONDITIONALLY at the FOOT of the block -- so a re-include under a
# test file or a runner is a silent no-op, while a FOREIGN definition of some other
# Phase-14 name can never make this body skip.
#
# A SEED IS THE WORST POSSIBLE SENTINEL IN THIS PROJECT, AND THAT IS NOT A THEORY.
# Phase 13 was originally keyed on :P13_DEV_SEED. Because every Phase-N pre-registration
# is obliged to RE-DECLARE prior phases' seeds in order to assert stream disjointness,
# spike/validation/p12_consts.jl legitimately mirrors P13_DEV_SEED -- and that mirror made
# a full-suite run skip the entire Phase-13 body, leaving ~104 Tier-1 constants undefined.
# The defect was on the guard side. See .planning/CONVENTIONS.md C-01 and
# 13-D15-AMENDMENT.md section 5B. THIS FILE MIRRORS P13_DEV_SEED AND P13_FIX_SEED BELOW,
# which is exactly why its sentinel is not a seed.
#
# EVERY CALLER USES THE GUARDED INCLUDE:
#     isdefined(@__MODULE__, :P14_DECLARED_DEVIATIONS) || include(joinpath(@__DIR__, "consts.jl"))

# The Phase-13 Tier-1 (and Tier-2) pre-registration, loaded READ-ONLY so the forbidden-seed
# inventory can be RECOMPUTED rather than retyped. This sits OUTSIDE the guard block for the
# same reason p13/consts.jl puts its `module _GC` outside its own: p13/consts.jl declares a
# top-level `module`, and idempotency is supplied by the guarded-include idiom instead.
isdefined(@__MODULE__, :P13_DECLARED_DEVIATIONS) ||
    include(joinpath(@__DIR__, "..", "p13", "consts.jl"))

if !isdefined(@__MODULE__, :P14_DECLARED_DEVIATIONS)
    import Random123: Philox4x     # counter-based RNG; the fresh disjoint Phase-14 streams

    # =====================================================================================
    # A. FORBIDDEN seeds and the FRESH Phase-14 streams
    # =====================================================================================
    # DERIVED, NEVER RETYPED. `_p13_forbidden()` already contains: zero, the six named
    # reserved streams (DEFAULT_MASTER_SEED, NPE_MASTER_SEED, VAL_MASTER_SEED, VAL_FIX_SEED,
    # RATIO_PAIR_SEED, CORPUS_MASTER_SEED), P11_DEV_SEED, all 25 burned spike keys, and both
    # RECOMPUTED ship-gate families (PROD_SEED, PROD_SEED_V2). Re-typing that inventory here
    # would create a second copy that can drift; recomputing it cannot.
    #
    # PHASE 13'S OWN TWO SEEDS ARE ADDED, because Phase 13 could not forbid itself: they are
    # absent from `_p13_forbidden()` by construction (p13/consts.jl asserts exactly that), and
    # they are the streams a Phase-14 collision is MOST likely to hit, since Phase 14 consumes
    # Phase 13's net and reads its artifacts.
    "The Phase-14 forbidden-seed list in declaration order, so an accidental duplicate stays visible."
    _p14_forbidden_list() = vcat(
        collect(UInt64, _p13_forbidden()),
        UInt64[UInt64(P13_DEV_SEED), UInt64(P13_FIX_SEED)],
    )

    """
        _p14_forbidden() -> Set{UInt64}

    Every seed a Phase-14 stream must not be: the whole RECOMPUTED Phase-13 forbidden set
    plus Phase 13's own reported and fixture seeds. `P14_DEV_SEED` and `P14_FIX_SEED` are
    asserted absent from this set at the foot of this block AND, independently, in
    `spike/test/test_p14_consts.jl`.
    """
    _p14_forbidden() = Set{UInt64}(_p14_forbidden_list())

    # The salt inventory of the whole repository, MIRRORED from P13_REPO_SALTS and extended
    # with P13_SALT itself -- which P13_REPO_SALTS does not contain, for the same
    # cannot-forbid-itself reason as the seeds above. A repeated (seed, salt) pair is a
    # repeated Philox key, hence a repeated stream.
    const P14_REPO_SALTS = (P13_REPO_SALTS..., P13_SALT)

    # --- The FRESH Phase-14 streams -------------------------------------------------------
    const P14_DEV_SEED = 0x0000_0000_0B14_DE71  # "0B14" = phase 14, "DE71" = DEV-1 (mirrors P11/P13)
    const P14_FIX_SEED = 0x0000_0000_0B14_F1F7  # FIXTURE stream, mirroring the VAL_FIX_SEED and
                                                # P13_FIX_SEED discipline: fixtures must NEVER
                                                # consume (and so never pre-observe) a reported
                                                # stream.
    const P14_SALT     = 0x5851_F42D_4C95_7F2D  # the MMIX LCG multiplier; asserted below to be
                                                # none of the nine entries of P14_REPO_SALTS,
                                                # rather than assumed to be new.

    "The Phase-14 REPORTED RNG: `Philox4x` keyed by `(P14_DEV_SEED xor P14_SALT, counter)`."
    p14_rng(counter::Integer = 0) =
        Philox4x(UInt64, (UInt64(P14_DEV_SEED) ⊻ P14_SALT, UInt64(counter)))

    "The Phase-14 FIXTURE RNG: same construction off `P14_FIX_SEED`, a provably disjoint key."
    p14_fix_rng(counter::Integer = 0) =
        Philox4x(UInt64, (UInt64(P14_FIX_SEED) ⊻ P14_SALT, UInt64(counter)))

    # RESERVED counters, ONE PER ACTIVITY, so two activities can never silently share a
    # sub-stream. The fixture counter rides the FIXTURE seed and is never a reported counter.
    const P14_CAL_COUNTER     = 1   # split-conformal CALIBRATION draws (D-02 hedge)
    const P14_EVAL_COUNTER    = 2   # the SHARED reported evaluation set: SC1-b, SC1-d and
                                    # SC3-a/b/c all read THIS one set. Sharing is deliberate --
                                    # a risk-coverage curve and the coverage it is scored at
                                    # must come from the same draws or they are not the same
                                    # experiment -- and it is recorded here so that the
                                    # multiplicity is visible rather than discovered later.
    const P14_OOD_COUNTER     = 3   # the SC3-d misspecification arm (matched in-distribution
                                    # vs misspecified, P14_N_OOD per arm)
    const P14_REAL_COUNTER    = 4   # the D-03a six-TIFF illustration (read-only real substrate)
    const P14_FIXTURE_COUNTER = 99  # FIXTURES ONLY -- never a reported counter

    const P14_REPORTED_COUNTERS = (P14_CAL_COUNTER, P14_EVAL_COUNTER,
                                   P14_OOD_COUNTER, P14_REAL_COUNTER)
    const P14_ALL_COUNTERS      = (P14_REPORTED_COUNTERS..., P14_FIXTURE_COUNTER)
    # Phase 13's complete counter family, named here so the cross-product below is exhaustive
    # rather than a sample. Kept in the P14_ namespace so this file never shadows a P13 name.
    const P14_P13_COUNTERS = (P13_TAU_COUNTER, P13_DATAGEN_COUNTER, P13_GATE_COUNTER,
                              P13_ALPHA_COUNTER, P13_CONTINUITY_COUNTER, P13_FIXTURE_COUNTER)

    # THE KEY WORDS THEMSELVES, not merely the seeds. Two seeds that differ are not enough:
    # `seed xor salt xor counter` is what actually keys Philox in this project's runners
    # (spike/p13/run_three_way_gate.jl:155-165), so THAT is the quantity that must be proven
    # distinct. "Obviously distinct" is how two lanes end up sharing a Philox key.
    "Every Philox FIRST key word Phase 14 can construct, in declaration order."
    _p14_key_words() = vcat(
        UInt64[UInt64(P14_DEV_SEED) ⊻ P14_SALT ⊻ UInt64(c) for c in P14_ALL_COUNTERS],
        UInt64[UInt64(P14_FIX_SEED) ⊻ P14_SALT ⊻ UInt64(c) for c in P14_ALL_COUNTERS],
        UInt64[UInt64(P14_DEV_SEED) ⊻ P14_SALT, UInt64(P14_FIX_SEED) ⊻ P14_SALT],
    )

    "Every Philox FIRST key word Phase 13 can construct, including its two bare key words."
    _p14_p13_key_words() = vcat(
        UInt64[UInt64(P13_DEV_SEED) ⊻ P13_SALT ⊻ UInt64(c) for c in P14_P13_COUNTERS],
        UInt64[UInt64(P13_FIX_SEED) ⊻ P13_SALT ⊻ UInt64(c) for c in P14_P13_COUNTERS],
        UInt64[UInt64(P13_DEV_SEED) ⊻ P13_SALT, UInt64(P13_FIX_SEED) ⊻ P13_SALT],
    )

    # =====================================================================================
    # B. Sample sizes
    # =====================================================================================
    # 14-RESEARCH section B.4, and this arithmetic is EXACT rather than a judgement: split
    # conformal's realized coverage on a fresh evaluation set is Binomial, so at n = 2000 and
    # alpha = 0.10 its standard deviation is sqrt(0.10 * 0.90 / 2000) = 0.0067. That is what
    # pins the SC1-d band; the band is DERIVED from alpha and n by formula and is therefore
    # not a value anyone gets to choose. Larger n is strictly better and costs almost nothing
    # here, so 2000 is a floor chosen for headroom, not a compromise.
    const P14_N_CAL  = 2000   # conformal calibration draws
    const P14_N_EVAL = 2000   # the shared reported evaluation set
    const P14_N_OOD  = 1000   # PER ARM, matched in-distribution vs misspecified (SC3-d)

    # =====================================================================================
    # C. The two alphas -- and they must NEVER be conflated (D-07)
    # =====================================================================================
    # THEY SHARE THE VALUE 0.10, AND THAT IS PRECISELY WHY NEITHER MAY EVER BE ASSIGNED FROM
    # THE OTHER. A binding written as `alpha_fdr = P14_ALPHA_CONFORMAL` would be invisible in
    # every test that checks values, and would silently fuse a distribution-free hedge level
    # to a user-facing error-rate parameter. test_p14_consts.jl greps the comment-stripped
    # source of this file to forbid exactly that assignment.
    #
    # alpha_conformal: the miscoverage level of the split-conformal hedge (D-02). It is a
    # PRE-REGISTERED CONSTANT and lives HERE, ruled by the user on 2026-08-03.
    const P14_ALPHA_CONFORMAL = 0.10

    # alpha_FDR: NOT A CONSTANT OF THIS FILE. SC1 says "user-set Bayesian FDR", so alpha_FDR
    # is a USER PARAMETER supplied per call to the decision function (D-07). What IS frozen is
    # the GRID that SC1-b sweeps when it demonstrates FDR control across levels -- the sweep
    # has to be pre-registered or "we controlled FDR" is a claim about a grid chosen after the
    # fact. The grid is not a default and is not a bar.
    const P14_ALPHA_FDR_GRID = (0.01, 0.05, 0.10, 0.20)

    # The SC1-c prior-sensitivity sweep. REPORTED, NOT GATED -- there is no bar attached to it
    # anywhere in this file, and if none exists the sweep cannot accidentally become a gate.
    # 0.217 is the class mass the Phase-13 label cut puts on `coloc`; the other four rungs
    # bracket it by roughly a factor of four in each direction.
    const P14_PI_COLOC_GRID = (0.05, 0.10, 0.217, 0.40, 0.70)

    # =====================================================================================
    # D. The class-order contract, frozen so it cannot drift
    # =====================================================================================
    # FIVE DIFFERENT CLASS ORDERINGS COEXIST IN THIS CODEBASE: `three_way_probs`,
    # `three_way_log_bf`, `ThreeWayLogBF`, `confusion_class_order` and `pi_class_masses` do
    # not agree on positional order. ECE, AUC and FDR are all invariant to a CONSISTENT
    # relabelling, which is exactly what makes an inconsistent one so dangerous: it produces
    # numbers that look right. EVERY READ SITE IN PHASE 14 MUST BE BY NAME -- never [1], [2]
    # or [3] -- and this tuple is the single frozen statement of what the names are.
    const P14_CLASS_KEYS = (:coloc, :random, :exclusion)

    # =====================================================================================
    # E. THE FOUR SC3 BARS AND THE COVERAGE FLOOR -- JUDGEMENT CALLS WITH NO DERIVATION
    # =====================================================================================
    # READ THIS BEFORE QUOTING ANY SC3 NUMBER.
    #
    # THE FIVE CONSTANTS IMMEDIATELY BELOW ARE JUDGEMENT CALLS WITH NO DERIVATION. Every one
    # of them is tagged [ASSUMED] in 14-RESEARCH section C.3, and its Assumptions Log A3 names
    # all five as judgement calls and recommends putting them to the user rather than letting
    # an executor freeze them silently. Only two carry any stated reasoning at all: the
    # selective-skill statistic is NORMALIZED between the random and oracle risk-coverage
    # curves so the bar is immune to the base error rate (the normalization is derived, the
    # 0.60 is not), and the coverage floor has a recorded REASON (below ~10-20% coverage the
    # selective-risk denominator is too small to mean anything) but not a derivation of 0.20
    # specifically. The other three are bare suggestions.
    #
    # WHY THAT IS SAID SO LOUDLY. This project has recorded FOUR separate cases of an
    # underived bar measuring something other than what it named: SC1g's component-vs-total,
    # the n=2-against-271 real arm, the Wald-labelled-Wilson sizing (DEF-12-04), and the
    # Phase-12 Stage-1 control ceiling -- a gating ceiling with no derivation anywhere in its
    # own pre-registration. STATE.md records the consequence: the count is now high enough
    # that "we pre-registered it" no longer settles an argument on its own.
    #
    # USER RULING, 2026-08-03. These five values were put to the user as a BLOCKING decision
    # checkpoint BEFORE this file was written -- before spike/p14/ existed at all -- and the
    # user ruled ACCEPT-AS-PROPOSED. The ruling is recorded verbatim in
    # .planning/phases/14-decision-and-abstention-layer/14-01-SUMMARY.md, in a section
    # APPENDED below the BLOCKED record that proves the question was asked first. No proposed
    # value was changed by the ruling; what the ruling establishes is that these five bars are
    # the USER'S, ratified before any Phase-14 result existed, rather than five numbers an
    # agent chose.
    #
    # THE BINDING CONSEQUENCE, IN THE USER'S OWN TERMS AND NOT SOFTENED:
    # IF A BAR BELOW IS MISSED, THE RESULT IS **REPORTED** AND **NOT RE-TUNED**. The bar is
    # not relaxed, not re-derived, not re-scoped and not swapped for a neighbouring statistic.
    # This matches the Phase-11 and Phase-12 precedent, where a missed pre-registered bar was
    # closed as a negative-but-useful finding rather than repaired into a pass. This file may
    # not be amended to relax any of them, and P14_ITERATION_ALLOWANCE below cannot be spent
    # on them: its one pre-declared trigger concerns conformal coverage and none of these
    # five.
    const P14_SKILL_FLOOR      = 0.60  # SC3-a: selective skill, normalized between the random
                                       # and oracle risk-coverage curves
    const P14_SPEARMAN_FLOOR   = 0.95  # SC3-b: Spearman(coverage, selective risk)
    const P14_COVERAGE_FLOOR   = 0.20  # SC3-b: the coverage range the Spearman is computed on
    const P14_AUC_HARD_FLOOR   = 0.80  # SC3-c: abstention concentrates on would-have-been-wrong
                                       # items
    const P14_OOD_MARGIN_FLOOR = 0.50  # SC3-d: abstain rate under strong misspecification minus
                                       # the in-distribution abstain rate

    # Enumerated as SYMBOLS so the report can list the judgement calls MECHANICALLY rather
    # than by hand. A bar that is quietly dropped from the report is a bar that stops being
    # visible as a judgement call; a hand-maintained list is how that happens.
    const P14_JUDGEMENT_CALL_BARS = (:P14_SKILL_FLOOR, :P14_SPEARMAN_FLOOR, :P14_COVERAGE_FLOOR,
                                     :P14_AUC_HARD_FLOOR, :P14_OOD_MARGIN_FLOOR)

    # =====================================================================================
    # F. The vacuity LABEL -- a label, never a verdict (14-RESEARCH Pitfall 9)
    # =====================================================================================
    # A decision layer that abstains on essentially nothing has not demonstrated selective
    # prediction; it has demonstrated that its trigger never fires. But "abstains rarely" can
    # also be the CORRECT behaviour on well-specified data, so this must NOT be a failure
    # condition. Falling below this rate produces a LABEL on the result (`vacuous_hedge =
    # true`), exactly as `vacuous_pass` is designed at spike/p13/result.jl:284-298, and the
    # label travels with every quoted SC3 number.
    const P14_AMBIGUOUS_RATE_FLOOR = 0.01

    # =====================================================================================
    # G. The iteration allowance
    # =====================================================================================
    # ONE documented iteration is authorised, and it is authorised HERE, before any result
    # exists. A SECOND iteration is NOT authorised by this file and cannot be authorised by
    # amending this file. Phase 7 was amended TWICE after seeing results and the credibility
    # cost of that is the reason this constant exists.
    const P14_ITERATION_ALLOWANCE = 1
    const P14_ITERATION_TRIGGER = """
    THE ONE PRE-DECLARED CONDITION for spending P14_ITERATION_ALLOWANCE: the realized
    conformal coverage on the reported evaluation set falls BELOW the SC1-d band, which is
    DERIVED from P14_ALPHA_CONFORMAL and P14_N_EVAL by formula and is therefore not itself
    adjustable. The allowance is then spent on switching the nonconformity score from LAC
    (least-ambiguous set-valued classifier) to APS (adaptive prediction sets) and re-running
    ONCE. It is NEVER spent on relaxing a bar after a result, and specifically NEVER on any
    member of P14_JUDGEMENT_CALL_BARS -- a missed bar there is REPORTED, not re-tuned. A
    second iteration is not authorised by this file.
    """

    # =====================================================================================
    # H. DECLARED deviations
    # =====================================================================================
    const P14_SRC_UNTOUCHED = true
    # Naming the deviations IN ADVANCE is the discipline. The report must enumerate exactly
    # these four and no others; anything else that changes is an UNDECLARED deviation.
    const P14_DECLARED_DEVIATIONS = (
        "SC1 is AMENDED by D-02: conformal sets are met in SUBSTANCE by a hand-rolled split-conformal hedge, and the named ConformalPrediction.jl is NOT added, because adding it would break the byte-frozen spike environment guard",
        "SC2 is AMENDED by D-05: the literal OR of three abstention triggers is replaced by an ASYMMETRIC fusion, because cross-method disagreement alone is the Phase-9 product thesis and must DECIDE, not silence the tool",
        "tau is INHERITED per D-07: loaded from the Phase-13 artifact with provenance asserted, never re-derived and never written into this file as a literal",
        "the real-data substrate is SUBSTITUTED per D-03a: the six committed microscopy TIFFs under test/test_images/ replace D-03's corpus rows, which carry no unsealed physical ground truth and whose bytes are unfetched, precisely so that Phase 16's blind corpus stays blind",
    )

    # =====================================================================================
    # EXECUTABLE self-checks (cheap; no inference, no simulation, no file writes)
    # =====================================================================================
    # Seed freshness against the RECOMPUTED, never-retyped forbidden inventory.
    @assert !(UInt64(P14_DEV_SEED) in _p14_forbidden()) "P14_DEV_SEED collides with a forbidden seed"
    @assert !(UInt64(P14_FIX_SEED) in _p14_forbidden()) "P14_FIX_SEED collides with a forbidden seed"
    @assert UInt64(P14_DEV_SEED) != UInt64(P14_FIX_SEED) "the fixture stream must differ from the reported stream"
    @assert !(P14_SALT in P14_REPO_SALTS) "P14_SALT reuses an existing repository salt"
    # A duplicate inside the forbidden list would silently shrink the set; catch that.
    @assert length(_p14_forbidden()) == length(_p14_forbidden_list()) "duplicate entry in the Phase-14 forbidden-seed list"
    # THE P13 RECOMPUTE ACTUALLY RAN. A comment claiming the Phase-13 inventory was derived is
    # not evidence that the include happened; this is.
    @assert UInt64(P13_DEV_SEED) in _p14_forbidden() "the Phase-13 forbidden-seed recompute did not run"
    @assert UInt64(P13_FIX_SEED) in _p14_forbidden() "the Phase-13 fixture seed was not added to the Phase-14 forbidden set"
    @assert any(v -> UInt64(v) in _p14_forbidden(), values(_GC.PROD_SEED_V2)) "the recomputed ship-gate seed family did not reach the Phase-14 forbidden set"
    # Fixtures must never ride a counter a reported number rides on.
    @assert P14_FIXTURE_COUNTER ∉ P14_REPORTED_COUNTERS "the fixture counter is also a reported counter"
    @assert length(Set(P14_ALL_COUNTERS)) == length(P14_ALL_COUNTERS) "two Phase-14 activities share a counter"

    # THE FULL PHILOX KEY-WORD CROSS-PRODUCT. Every Phase-14 key word against every Phase-13
    # key word, and every Phase-14 key word against every other. A shared key between a
    # Phase-13 reported stream and a Phase-14 evaluation stream would mean Phase 14 scored
    # itself on draws a Phase-13 number already rode -- silently, and undetectably from the
    # artifacts alone.
    for kw14 in _p14_key_words(), kw13 in _p14_p13_key_words()
        @assert kw14 != kw13 "a Phase-14 Philox key word collides with a Phase-13 key word: 0x$(string(kw14, base = 16, pad = 16))"
    end
    @assert length(Set(_p14_key_words())) == length(_p14_key_words()) "two Phase-14 Philox key words collide with each other"

    # Sample sizes, alphas and grids.
    @assert P14_N_CAL > 0 && P14_N_EVAL > 0 && P14_N_OOD > 0 "a Phase-14 sample size is non-positive"
    @assert 0.0 < P14_ALPHA_CONFORMAL < 1.0 "P14_ALPHA_CONFORMAL is not a miscoverage level"
    @assert issorted(P14_ALPHA_FDR_GRID) && first(P14_ALPHA_FDR_GRID) > 0.0 &&
            last(P14_ALPHA_FDR_GRID) < 1.0 "P14_ALPHA_FDR_GRID is not an ordered grid of levels"
    @assert issorted(P14_PI_COLOC_GRID) && first(P14_PI_COLOC_GRID) > 0.0 &&
            last(P14_PI_COLOC_GRID) < 1.0 "P14_PI_COLOC_GRID is not an ordered grid of prior masses"
    # The two alphas are SEPARATE BINDINGS. They share a value; they must not share a name or
    # a derivation, and the grid must not have been silently reduced to the conformal level.
    @assert P14_ALPHA_CONFORMAL in P14_ALPHA_FDR_GRID "the shared 0.10 value is no longer shared, so the D-07 conflation warning above is stale and must be re-read"
    @assert length(P14_ALPHA_FDR_GRID) > 1 "P14_ALPHA_FDR_GRID collapsed to a single level and can no longer demonstrate control ACROSS levels"

    # The class-order contract.
    @assert P14_CLASS_KEYS === (:coloc, :random, :exclusion) "the frozen class-key order changed"
    @assert length(Set(P14_CLASS_KEYS)) == 3 "P14_CLASS_KEYS contains a duplicate class name"

    # The judgement-call bars: in range, enumerated, and all five actually defined.
    @assert 0.0 < P14_SKILL_FLOOR < 1.0 "P14_SKILL_FLOOR is not a normalized skill"
    @assert 0.0 < P14_SPEARMAN_FLOOR <= 1.0 "P14_SPEARMAN_FLOOR is not a correlation"
    @assert 0.0 < P14_COVERAGE_FLOOR < 1.0 "P14_COVERAGE_FLOOR is not a coverage"
    @assert 0.5 < P14_AUC_HARD_FLOOR < 1.0 "P14_AUC_HARD_FLOOR is at or below chance"
    @assert 0.0 < P14_OOD_MARGIN_FLOOR < 1.0 "P14_OOD_MARGIN_FLOOR is not an abstain-rate margin"
    @assert 0.0 < P14_AMBIGUOUS_RATE_FLOOR < 1.0 "P14_AMBIGUOUS_RATE_FLOOR is not a rate"
    @assert length(P14_JUDGEMENT_CALL_BARS) == 5 "the judgement-call enumeration no longer lists five bars"
    @assert all(s -> isdefined(@__MODULE__, s), P14_JUDGEMENT_CALL_BARS) "a name in P14_JUDGEMENT_CALL_BARS is not defined"

    # The iteration allowance is pre-declared, non-trivial and names its ONE remedy.
    @assert P14_ITERATION_ALLOWANCE == 1 "the iteration allowance is no longer exactly one"
    @assert P14_ITERATION_TRIGGER isa AbstractString && !isempty(strip(P14_ITERATION_TRIGGER)) "the iteration trigger is empty"
    @assert occursin("APS", P14_ITERATION_TRIGGER) "the iteration trigger no longer names the APS remedy"

    # Decoupling and the fixed-length declared-deviation enumeration.
    @assert P14_SRC_UNTOUCHED == true "Phase 14 does not touch src/ (D-01)"
    @assert length(P14_DECLARED_DEVIATIONS) == 4 "the declared-deviation enumeration changed length"
    @assert all(d -> d isa AbstractString && !isempty(strip(d)), P14_DECLARED_DEVIATIONS) "a declared deviation is empty"

    # TIER 2 IS INHERITED, NOT DECLARED HERE. D-07 requires tau to be LOADED from the
    # Phase-13 artifact with its provenance asserted. A Phase-14 literal for it -- under any
    # name -- is the exact divergence D-07 exists to prevent, so its absence is asserted here
    # and, independently, in spike/test/test_p14_consts.jl by BOTH `isdefined` and a
    # comment-stripped source grep. Note that the Phase-13 value IS in scope, because this
    # file includes p13/consts.jl above; that is the point. It is read from there, at the one
    # place it was measured, and never mirrored into this namespace.
    # THE NEEDLE IS ASSEMBLED AT RUN TIME FROM FRAGMENTS, and that is deliberate rather than
    # coy. The freeze's own acceptance criterion is a source grep over the comment-stripped
    # text of this file, so a contiguous literal spelling of the forbidden name -- written
    # here in order to FORBID it -- would match, and this assertion would fail the file it is
    # protecting. Same construction, same reason, as the three needles
    # spike/test/test_p12_decoupling.jl builds by concatenation so its source scan can scan
    # itself.
    @assert !isdefined(@__MODULE__, Symbol("P14_", "TAU")) "a Phase-14 tau literal exists; D-07 requires tau to be inherited from the Phase-13 artifact by loading, never declared here"
end
