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

# spike/p13/real_images.jl --- D-15 (AMENDED) real-image ingestion. REPORTED, NOT GATED. QUALITATIVE ONLY.
#
# (a) WHAT THIS IS. THE SIX COMMITTED MICROSCOPY TIFFS UNDER test/test_images/ --
#     positive_c{1,2,3}.tif AND negative_c{1,2,3}.tif, 1028 x 1376, RGB{N0f16} ON DISK -- ARE THE
#     AMENDED D-15 SUBSTRATE. THEY CARRY NO COLOCALIZATION GROUND-TRUTH LABEL. ANYTHING READ
#     THROUGH THIS FILE IS A BEHAVIOUR CHECK, NEVER A CORRECTNESS CHECK. THE ARM CAN SHOW THAT
#     THREE-WAY VERDICTS BEHAVE SENSIBLY ON REAL MICROSCOPY; IT CANNOT SHOW THAT THEY ARE
#     CORRECT. NO REAL-IMAGE QUANTITY CARRIES A PASS/FAIL THRESHOLD ANYWHERE IN THIS PHASE
#     (P13_REAL_IS_GATED = false, P13_REAL_QUALITATIVE_ONLY = true). LABELLED REAL SEGREGATION
#     VALIDATION IS DEFERRED TO THE PHASE-16 BLIND EVALUATION.
#
# (b) THE NAMING CORRECTION (P13_REAL_NAMING_CORRECTION, STATED HERE AND NOT ONLY IN THE REPORT).
#     THE FOLDER NAMES positive/ AND negative/ ARE THE ORIGINAL PACKAGE'S *BIOLOGICAL* TEST
#     CONDITIONS, NOT COLOCALIZATION LABELS. THE "negative" PAIR MEASURES MEAN PATCH CORRELATION
#     +0.2481 AND IS THEREFORE A POSITIVELY CORRELATED PAIR, NOT AN ANTI-CORRELATED ONE. NOTHING
#     IN THIS PHASE MAY TREAT negative/ AS AN EXCLUSION EXAMPLE, AND A READER WHO ASSUMES
#     OTHERWISE WILL MISREAD EVERY REAL-IMAGE FIGURE IN THE PHASE.
#
# (b2) THE CHANNEL-PAIR CAVEAT, RECORDED HONESTLY BECAUSE IT LANDED AFTER THE PRE-REGISTRATION.
#     P13_REAL_CHANNEL_PAIR IS FROZEN AT c1/c2. test/runtests.jl:105 records the fixture channels
#     as ["blue", "green", "red"], so c1 is the DAPI/Hoechst NUCLEAR COUNTERSTAIN, not a target
#     protein -- and spike/simulator/ghat.jl now records (second hand-patch entry, 2026-07-25)
#     that the c1/c2 figures are SUPERSEDED for colocalization purposes by the c2/c3 = green/red
#     pair (unmasked positive mu = 0.4603, negative mu = 0.3815). THE CONSEQUENCE FOR THIS FILE IS
#     NARROW AND MUST NOT BE OVERSTATED: the c1/c2 anchor regression below is a READ-CHAIN
#     IDENTITY check -- it proves this ingestion is bit-for-bit the same load path the frozen
#     ghat calibration used -- and it is NOT a claim that c1/c2 measures colocalization. Every
#     entry point here takes `channels` as a keyword precisely so the reported run can read the
#     coloc pair as well; the frozen pre-registration is honoured, not silently reinterpreted.
#
# (c) THE TARGET-SUBSTITUTION RECORD (P13_REAL_SUBSTITUTION_RECORD, 11-10-PLAN.md HOUSE SHAPE).
#     D-15 ORIGINALLY NAMED THE PHYSICAL MITOCHONDRIA ANCHOR OF THE FROZEN PROVENANCE MANIFEST AS
#     THE REAL-DATA CHECK. VERIFICATION SHOWED BOTH physical-primary ROWS ARE
#     sha256 = "PENDING-FETCH", bytes = 0 AND split = sealed_holdout -- UNFETCHED, AND RESERVED
#     BEHIND THE ANTI-SNOOPING CONTROL FOR THE PHASE-16 BLIND EVALUATION. CONSUMING ONE HERE WOULD
#     IRREVERSIBLY BURN PHASE 16 ON THE VERY HYPOTHESIS PHASE 16 EXISTS TO EVALUATE. THE
#     SANCTIONED SUBSTITUTE IS THE SIX COMMITTED MICROSCOPY TIFFS UNDER test/test_images/. THE
#     SEALED-HOLDOUT ACCESSOR IS NOT CALLED -- NOT IN A SCRIPT, NOT IN A TEST, NOT BEHIND A FLAG
#     -- AND NO SEAL-BREAK ESCAPE HATCH EXISTS OR MAY BE RE-INTRODUCED.
#
# (d) READ-ONLY DISCIPLINE. test/test_images/ IS READ-ONLY INPUT. THIS FILE OPENS THE TIFFS AND
#     NOTHING ELSE: IT NEVER WRITES, MOVES OR MODIFIES ANYTHING UNDER test/, NO test/ PATH APPEARS
#     IN ANY PHASE-13 files_modified LIST, AND verify_real_readonly_digest() BELOW MAKES THAT
#     STRUCTURAL RATHER THAN MERELY INTENDED.
#
# (e) THE LINEAGE CLAIM. THESE ARE NOT ARBITRARY TIFFS. spike/simulator/calibration.jl:156-164
#     `_real_anchor()` ALREADY LOADS THESE FIXTURES THROUGH THE SAME src/LoadImages.jl LOADER, AND
#     spike/simulator/ghat.jl:58-64 RECORDS THE RESULTING c1/c2 ANCHORS AS positive mu = 0.3292,
#     negative mu = 0.2481. THIS FILE REPRODUCES THOSE NUMBERS (MEASURED 0.32916 / 0.24805), WHICH
#     IS WHAT MAKES THE REAL SUBSTRATE A LINEAGE-CARRYING REFERENCE POINT IN THIS PROJECT'S OWN
#     SIMULATOR CALIBRATION RATHER THAN AN ARBITRARY SET OF FILES.
#
# (f) DECOUPLING (hard constraint, CLAUDE.md): spike-local. src/ is reached ONLY read-only,
#     transitively through spike/contract.jl's include chain, for the frozen
#     MultiChannelImage(name, paths, channels) -> load_tiff loader and the unchanged
#     patch()/correlation() summary. No src/ byte, no manifest byte and no dependency is touched;
#     no path under the sealed provenance directory is constructed, opened or referenced by any
#     executable line; the seal is never opened by this phase. The literal include line is
#     `include(joinpath(@__DIR__, "..", "contract.jl"))`, kept below with the house guard idiom.
#
# NO NEW LOADER, NO CROP, NO RESIZE, NO PAD, NO num_patches CHANGE, NO MASKING, NO BACKGROUND
# SUBTRACTION. A channel-sum, a max-projection or a hand-rolled RGB-to-luminance reduction would
# silently move off the frozen anchor. `_apply_mask!` zeroes sub-threshold pixels and
# `_exclude_zero` then DELETES them, which would change the very statistic the net was trained
# on -- so neither the masking helper nor any resize helper is called anywhere below, and a
# source grep asserts their absence. Otsu is used ONLY to define the D-16 object mask, never to
# modify pixels that feed `patch_summary`.
#
# ===== END OF EXPLANATORY HEADER -- SEAL SELF-CHECK SENTINEL, DO NOT REWORD =====

using Statistics     # mean (the mask fraction)
using SHA            # sha256 of the git index listing (Julia standard library, no dependency)

# ORDER MATTERS, and it is the documented one. The pre-registration first, then spike/contract.jl
# -- which puts `using StatsBase; using Statistics; using Images` in scope BEFORE it
# read-only-includes src/LoadImages.jl and src/colocalization.jl, the load-order landmine
# documented at spike/contract.jl:34-39 -- then the 128-row encoding the m-bar statistic is read
# out of, then the substrate-agnostic alpha transform. Guarded for idempotency under a runner.
isdefined(@__MODULE__, :P13_REAL_IS_GATED) || include(joinpath(@__DIR__, "consts.jl"))
isdefined(@__MODULE__, :build_mci)         || include(joinpath(@__DIR__, "..", "contract.jl"))
isdefined(@__MODULE__, :encode_d01)        || include(joinpath(@__DIR__, "..", "data", "encode.jl"))
isdefined(@__MODULE__, :alpha_segregate)   || include(joinpath(@__DIR__, "alpha_series.jl"))

# --- (0) The MACHINE-CHECKED seal assertion ---------------------------------------------------
# The 11-10-PLAN.md:191-195 mechanism: this file reads its OWN text, strips the explanatory
# header above (everything up to and including the sentinel comment line), and asserts the
# remainder names neither the sealed provenance directory nor the sealed-holdout accessor. That
# turns "we intended not to reach the seal" into "the include() FAILS if we did".
#
# BOTH FORBIDDEN NAMES ARE ASSEMBLED FROM FRAGMENTS AT THE ASSERTION SITE, ON PURPOSE, so the
# assertion does not itself embed the literals it forbids -- otherwise the check would trip on
# its own text and would have to be weakened into uselessness.
let src = read(@__FILE__, String)
    marker = "END OF EXPLANATORY HEADER"
    hit  = findfirst(marker, src)
    body = hit === nothing ? src : src[(last(hit) + 1):end]
    forbidden_dir = join(('c', 'o', 'r', 'p', 'u', 's'))
    forbidden_fn  = "open" * "_sealed" * "_holdout"
    occursin(forbidden_dir, body) && error(
        "real_images.jl: the sealed provenance directory is named in executable code. The " *
        "amended D-15 substitutes test/test_images/ precisely so that Phase 16's blind " *
        "evaluation is not burned here; remove the reference.")
    occursin(forbidden_fn, body) && error(
        "real_images.jl: the sealed-holdout accessor is named in executable code. It is never " *
        "called -- not in a script, not in a test, not behind a flag -- and no seal-break " *
        "escape hatch may be re-introduced.")
end

# --- (1) Repo-root resolution WITHOUT mutating the process CWD --------------------------------

"""
    P13_REPO_ROOT

The repository root, resolved statically from this file's own location.

The depth from `spike/p13/` is IDENTICAL to the depth from `spike/simulator/`, so this is
literally the idiom `spike/simulator/calibration.jl:59` already proves (`const _REPO_ROOT =
normpath(joinpath(@__DIR__, "..", "..")))`) and which `_real_anchor()` at
`spike/simulator/calibration.jl:156-164` already uses to read these exact fixtures.

The process working directory is DELIBERATELY never mutated. `test/runtests.jl:21` changes
directory before loading the fixtures; copying that here would make an include-time side effect
out of a path lookup, which the harness does not expect and which would leak into every other
suite running in the same session. A source grep asserts the CWD-mutating call is absent.
"""
const P13_REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

"""
    real_tif(cond::AbstractString, i::Integer) -> String

The absolute path of one committed fixture channel, `test/test_images/<cond>/<cond>_c<i>.tif`,
anchored at [`P13_REPO_ROOT`](@ref) and therefore independent of the caller's working directory.

`cond` is validated against the pre-registered `P13_REAL_CONDITIONS` and `i` against the three
channels the fixtures carry, BEFORE the filesystem is touched -- the explicit-`ArgumentError`-at-
the-boundary house style of `src/amortized/api.jl:111-120`. A missing file raises an
`ArgumentError` naming the expected absolute path rather than surfacing a bare `SystemError`
from deep inside the image loader.

READ-ONLY. This function returns a path to open; nothing in this phase writes, moves or modifies
anything under `test/`.
"""
function real_tif(cond::AbstractString, i::Integer)
    cond in P13_REAL_CONDITIONS || throw(ArgumentError(
        "real_tif: condition must be one of the pre-registered P13_REAL_CONDITIONS " *
        "$(P13_REAL_CONDITIONS), got \"$cond\". The amended D-15 substrate is exactly the six " *
        "committed microscopy TIFFs under test/test_images/; no other real substrate is " *
        "reachable from this phase."))
    i in 1:3 || throw(ArgumentError(
        "real_tif: channel index must be 1, 2 or 3 (the fixtures carry three channels, " *
        "recorded at test/runtests.jl:105 as [\"blue\", \"green\", \"red\"]), got $i"))
    p = joinpath(P13_REPO_ROOT, "test", "test_images", String(cond), "$(cond)_c$(i).tif")
    isfile(p) || throw(ArgumentError(
        "real_tif: expected the committed fixture at $p, but no file is there. The six TIFFs " *
        "are READ-ONLY INPUT committed to this repository; if one is missing the working tree " *
        "is not the shipped tree."))
    return p
end

# --- (2) Ingestion through the FROZEN constructor ----------------------------------------------

"""
    load_real(cond; channels = P13_REAL_CHANNEL_PAIR) -> MultiChannelImage

Load one condition's channels through the FROZEN convenience constructor
`MultiChannelImage(name, paths, channels)` (`src/LoadImages.jl:432-455`).

THE ARGUMENT ORDER IS `(name, paths, channels)`, NOT `(paths, channels, name)`. Channel names are
ALWAYS passed: an empty `channels` vector makes the constructor emit a `@warn` into the run log
and invent default names, which would put a spurious warning into every reported real-image run.

THE LOADER IS `Float64.(Images.Gray.(Images.load(path)))` (`src/LoadImages.jl:214-217`) and
nothing else. The on-disk element type is `RGB{N0f16}`; `Images.Gray` applies the Rec-601 luma
weighting; the result is a `Matrix{Float64}` in `[0, 1]` that is neither z-scored nor
background-subtracted. The frames are dim -- the first fixture channel measures min 6.1036e-5,
max 0.19794, mean 0.014986 -- so nothing sits anywhere near the top of the range.

NO MASKING AND NO BACKGROUND SUBTRACTION IS APPLIED ANYWHERE ON THE AMORTIZED PATH, here or
downstream. Otsu appears only where the shipped code already puts it (at construction, and in the
D-16 object-mask rule), never as a modification of the pixels that feed `patch_summary`.

Asserts the measured frame size `P13_REAL_IMSIZE` and `Float64` element type on every channel,
with errors naming the measured expectation, so an unexpected frame cannot slide silently into
the fixed grid.
"""
function load_real(cond::AbstractString; channels = P13_REAL_CHANNEL_PAIR)
    paths = String[real_tif(cond, i) for i in channels]
    names = String[string("ch", i) for i in channels]
    mci = MultiChannelImage(String(cond), paths, names)
    for (k, d) in enumerate(mci.data)
        size(d) == P13_REAL_IMSIZE || throw(ArgumentError(
            "load_real: channel $k of \"$cond\" has frame size $(size(d)), but the committed " *
            "fixtures measure $(P13_REAL_IMSIZE). NO CROP, RESIZE OR PAD STEP EXISTS OR MAY BE " *
            "ADDED here: adding one would be a new unregistered preprocessing choice and would " *
            "diverge from the frozen calibration read path."))
        eltype(d) == Float64 || throw(ArgumentError(
            "load_real: channel $k of \"$cond\" has element type $(eltype(d)), expected " *
            "Float64. The frozen loader is Float64.(Images.Gray.(Images.load(path))); a " *
            "different element type means a different read chain."))
    end
    return mci
end

"""
    real_pair(cond; channels = P13_REAL_CHANNEL_PAIR) -> Vector{Matrix{Float64}}

The two-channel matrix pair extracted from [`load_real`](@ref).

THIS IS THE INTERFACE BOUNDARY THAT MAKES ONE CODE PATH SERVE BOTH SUBSTRATES (D-15).
`simulate_pair` (`spike/simulator/forward.jl:103`) and `load_tiff.(paths)`
(`src/LoadImages.jl:433`) both already emit exactly this type, so `alpha_segregate` is typed on
`AbstractMatrix{Float64}` and never learns where its matrix came from. A provenance branch inside
the transform would make the simulated and real arms different experiments and their rungs
incomparable, which is why the substrate-dependent knobs live in [`real_provenance`](@ref)
instead.
"""
function real_pair(cond::AbstractString; channels = P13_REAL_CHANNEL_PAIR)
    length(channels) == 2 || throw(ArgumentError(
        "real_pair: the summary contract is a TWO-channel pair, got $(length(channels)) " *
        "channels ($channels). Read a redundancy arm as a second pair, never as a wider stack."))
    mci = load_real(cond; channels = channels)
    return Vector{Matrix{Float64}}([Matrix{Float64}(mci.data[1]), Matrix{Float64}(mci.data[2])])
end

# --- (3) The m-bar statistic (the SAME one alpha_ladder_summaries reports) ----------------------

"""
    real_mbar_of(pair) -> Float64

The mask-weighted mean of the PRESENT continuous summary rows of
`encode_d01(patch_summary(build_mci(pair)))` -- the SAME statistic `alpha_ladder_summaries`
(`spike/p13/alpha_series.jl`) reports per rung, so a rung and an anchor are directly comparable.

Rows 1:G^2 are the imputed per-patch correlations and rows G^2+1:2G^2 the parallel present/absent
mask (`spike/data/encode.jl:56-75`). THE WEIGHTING IS NOT COSMETIC: `encode_d01` imputes a
`missing` patch as 0.0, so an UNWEIGHTED mean would silently pull the statistic toward 0 exactly
where patches drop out. Returns `NaN` when every patch is missing, mirroring `induced_mu`'s
NaN-not-throw degeneracy idiom (`spike/contract.jl:104`).

Takes a raw pair rather than a condition so a transposed or otherwise derived pair can be
measured with the identical statistic (the A14 falsifier does exactly that).
"""
function real_mbar_of(pair)
    enc  = encode_d01(patch_summary(build_mci(Vector{Matrix{Float64}}(pair))))
    n    = length(enc) ÷ 2
    vals = @view enc[1:n]
    msk  = @view enc[(n + 1):(2n)]
    wsum = sum(msk)
    return wsum == 0.0 ? NaN : sum(vals .* msk) / wsum
end

"""
    real_mbar(cond; channels = P13_REAL_CHANNEL_PAIR, G = 8) -> Float64

The measured m-bar of one real condition, through the frozen ingestion.

THE FROZEN LINEAGE EXPECTATION IS `P13_REAL_ANCHOR_MBAR` -- positive +0.3292, negative +0.2481,
recorded at `spike/simulator/ghat.jl:58-64` from `_real_anchor()`'s unmasked c1/c2 read and
reproduced here to `P13_REAL_ANCHOR_TOL`. That regression is a READ-CHAIN IDENTITY check: it
proves this ingestion is the same load path the frozen calibration used, not a lookalike. See
the header caveat (b2) -- ghat.jl now records the c1/c2 figures as superseded *for
colocalization purposes* by the c2/c3 green/red pair, which does not weaken the identity claim
and does not license reinterpreting the c1/c2 number as a colocalization measurement.

`G` is accepted so the grid is stated at the call site rather than assumed by a reader, but the
frozen summary contract `patch_summary` (`spike/contract.jl:85`) IS the 8x8 D-10 grid and Phase
13 introduces no second grid; any other value raises rather than silently measuring something
else.
"""
function real_mbar(cond::AbstractString; channels = P13_REAL_CHANNEL_PAIR, G::Integer = 8)
    G == 8 || throw(ArgumentError(
        "real_mbar: the frozen summary contract patch_summary (spike/contract.jl:85) is the " *
        "8x8 D-10 patch grid and takes no grid argument; Phase 13 changes no num_patches and " *
        "introduces no second grid, so G = $G is not reachable from this arm."))
    return real_mbar_of(real_pair(cond; channels = channels))
end

# --- (4) Grid provenance: DERIVED truncation, never a literal ----------------------------------

"""
    real_grid_provenance(; G = 8, imsize = P13_REAL_IMSIZE) -> NamedTuple

The grid fit of the real frame, DERIVED by the same integer division `patch()` performs
(`src/colocalization.jl:37-43`), never transcribed as a literal. Fields: `imsize`, `patch_size`,
`rows_dropped`, `cols_dropped`, `px_per_patch`, `truncation_note`.

`patch()` TRUNCATES the frame to the largest exact multiple of the grid before patching, so a
reader who divides the row extent by the grid and expects a fractional patch size is reading the
wrong operation. The column axis divides exactly; the small row remainder is dropped; every patch
still holds three orders of magnitude more pixels than the 15-survivor floor of `_exclude_zero`,
and zero patches went `missing` on either fixture. `truncation_note` is the one sentence the
report must carry, assembled from the derived numbers so it cannot drift away from them.

NO CROP, RESIZE OR PAD STEP EXISTS ANYWHERE IN THIS ARM. Adding one would be a new unregistered
preprocessing choice and would diverge from the frozen calibration read path whose anchors this
ingestion regresses against.

TRANSPOSITION FOOTNOTE. `IMSIZE_SET` (`spike/data/seeding.jl:52-58`) contains this frame size
WITH THE AXES SWAPPED, at weight 0.03 -- the training distribution deliberately includes it. The
per-patch pixel count and the truncation are identical either way; only the column-major `vec()`
ordering of the grid corresponds to a transposed spatial layout, and the simulator's fields are
isotropic with a symmetric `SHIFT_PRIOR`, so this is statistically immaterial. It is recorded as
a footnote rather than left as a silent difference, and the A14 falsifier checks it in one line.
"""
function real_grid_provenance(; G::Integer = 8, imsize = P13_REAL_IMSIZE)
    G > 0 || throw(ArgumentError("real_grid_provenance: the grid must be positive, got G = $G"))
    rows, cols = imsize
    psx = rows ÷ G
    psy = cols ÷ G
    (psx > 0 && psy > 0) || throw(ArgumentError(
        "real_grid_provenance: a $(rows) x $(cols) frame cannot carry a $(G) x $(G) grid"))
    dropped_rows = rows - G * psx
    dropped_cols = cols - G * psy
    pct = round(100 * dropped_rows / rows; digits = 2)
    note = string("patch() truncates each frame to the largest exact multiple of the grid, so ",
                  "at G = ", G, " the last ", dropped_rows, " of ", rows, " rows (", pct,
                  "%) and the last ", dropped_cols, " of ", cols, " columns of each ", rows,
                  " x ", cols, " frame are not used; every patch holds ", psx * psy,
                  " px and zero patches go missing. No crop, resize or pad step exists.")
    return (imsize = (rows, cols), patch_size = (psx, psy),
            rows_dropped = dropped_rows, cols_dropped = dropped_cols,
            px_per_patch = psx * psy, truncation_note = note)
end

# --- (5) Substrate provenance: where every substrate-dependent knob lives -----------------------

"""
    real_provenance(cond; channels = P13_REAL_CHANNEL_PAIR) -> NamedTuple

The substrate provenance record for one real condition. Fields: `substrate`, `condition`,
`channels`, `paths`, `imsize`, `otsu_ch1`, `mask_fraction`, `background_floor`,
`source_zero_count`, `alpha0_mbar`, `is_gated`, `qualitative_only`, `read_only`.

THE THREE SUBSTRATE-DEPENDENT KNOBS LIVE HERE, IN PROVENANCE, NEVER IN THE ALGORITHM: the
background floor, the frame size, and what `alpha = 0` MEANS. On simulated substrate the ladder
runs random -> near-exclusion by construction; on this substrate `alpha = 0` is a MODERATELY
COLOCALIZED pair. The two arms are reported separately and are never averaged -- they are both
informative, and they are not the same experiment.

`source_zero_count` is `count(iszero, ...)` on the UNMODIFIED second channel and is reported
because it is the reference the corrected zero invariant compares against. The committed fixtures
already carry a handful of exact-zero pixels BEFORE any alpha is applied (2 px in `positive_c2`,
3 px in `negative_c2`, of 1 414 528), which is exactly why the pre-registered rule is
`P13_ALPHA_ZERO_INVARIANT = :count_iszero_preserved` -- "alpha introduces no NEW zero" -- and
never an all-pixels-strictly-positive form, whose naive assertion would fail a CORRECT algorithm
at `alpha = 0` where the output IS the input.

`is_gated`, `qualitative_only` and `read_only` are carried on every record so a downstream reader
that quotes a real-image number also carries the sentinel saying it is not a gate.
"""
function real_provenance(cond::AbstractString; channels = P13_REAL_CHANNEL_PAIR)
    mci  = load_real(cond; channels = channels)
    pair = Vector{Matrix{Float64}}([Matrix{Float64}(mci.data[1]), Matrix{Float64}(mci.data[2])])
    M    = ch1_object_mask(pair)
    return (substrate         = :real,
            condition         = String(cond),
            channels          = Tuple(channels),
            paths             = String[p for p in mci.path],
            imsize            = size(pair[1]),
            otsu_ch1          = mci.otsu_threshold[1],
            mask_fraction     = mean(M),
            background_floor  = alpha_background_floor(pair[2]),
            source_zero_count = count(iszero, pair[2]),
            alpha0_mbar       = real_mbar_of(pair),
            is_gated          = P13_REAL_IS_GATED,
            qualitative_only  = P13_REAL_QUALITATIVE_ONLY,
            read_only         = P13_REAL_READ_ONLY)
end

# --- (6) The real-substrate alpha ladder: substrate only, NO transform logic --------------------

"""
    real_alpha_ladder(cond; channels = P13_REAL_CHANNEL_PAIR, grid = P13_REAL_ALPHA_GRID)
        -> NamedTuple

Run the D-16 alpha ladder on one real condition through the SAME `alpha_segregate` /
`alpha_ladder` / `alpha_ladder_summaries` code path the simulated arm uses
(`spike/p13/alpha_series.jl`).

THIS FUNCTION ADDS NO TRANSFORM LOGIC OF ITS OWN. Its entire job is to supply the real substrate
and the substrate provenance; every arithmetic step, every guard and every invariant is the
shared transform's. That is what makes the real and simulated ladders comparable rung for rung
rather than two similar-looking curves, and `grid` defaults to `P13_REAL_ALPHA_GRID`, which the
pre-registration pins to the simulated `P13_ALPHA_LADDER`.

Two pre-registered guards fire BEFORE the ladder runs, as named errors rather than as an
uninterpretable curve:

  * the mask fraction must lie inside `P13_ALPHA_MASK_FRACTION_BOUNDS`. Outside that band the
    mask complement cannot absorb the redistributed mass and the boost blows up, so the
    "segregation" signal is confounded with a global rescale of a shrinking region. Measured on
    these fixtures: 13.49% and 22.92%, comfortably inside.
  * the background floor must be strictly positive. `_exclude_zero` (`src/colocalization.jl`)
    DROPS any pixel where either channel is zero, so a floor of 0.0 would DELETE the masked
    pixels from the correlation instead of anti-correlating them and the ladder would read flat
    for a reason that has nothing to do with segregation.

Returns `(condition, channels, ladder, rungs, mask_fraction, background_floor,
source_zero_count, provenance)` where `rungs` carries the per-rung
`(alpha, mbar, n_missing, max_value, sum_ratio, n_zero)` rows. REPORTED, NOT GATED: no threshold
is applied to any of it.
"""
function real_alpha_ladder(cond::AbstractString; channels = P13_REAL_CHANNEL_PAIR,
                           grid = P13_REAL_ALPHA_GRID)
    pair = real_pair(cond; channels = channels)
    M    = ch1_object_mask(pair)
    frac = mean(M)
    (first(P13_ALPHA_MASK_FRACTION_BOUNDS) <= frac <= last(P13_ALPHA_MASK_FRACTION_BOUNDS)) ||
        throw(ArgumentError(
            "real_alpha_ladder: the ch1 object-mask fraction of \"$cond\" is $frac, outside the " *
            "pre-registered band $(P13_ALPHA_MASK_FRACTION_BOUNDS). Outside that band the mask " *
            "complement cannot absorb the redistributed ch2 mass, the boost blows up, and the " *
            "resulting ladder is uninterpretable rather than merely noisy."))
    b = alpha_background_floor(pair[2])
    b > 0.0 || throw(ArgumentError(
        "real_alpha_ladder: the background floor of \"$cond\" is $b. It must be STRICTLY " *
        "POSITIVE because _exclude_zero (src/colocalization.jl) DROPS any pixel where either " *
        "channel is zero: a floor of 0.0 deletes the masked pixels from the correlation instead " *
        "of anti-correlating them, and the ladder then reads flat for the wrong reason."))
    rungs = alpha_ladder_summaries(pair; ladder = grid)
    return (condition         = String(cond),
            channels          = Tuple(channels),
            ladder            = grid,
            rungs             = rungs,
            mask_fraction     = frac,
            background_floor  = b,
            source_zero_count = count(iszero, pair[2]),
            provenance        = real_provenance(cond; channels = channels))
end

# --- (7) The structural read-only proof --------------------------------------------------------

"""
    verify_real_readonly_digest() -> String

The sha256 hex digest of `git ls-files -s test/test_images`, run with the repository root as the
subprocess working directory (`Cmd(...; dir = ...)`, so the calling process's own directory is
never mutated).

THIS IS THE STRUCTURAL READ-ONLY PROOF. `git ls-files -s` prints the index mode, blob hash, stage
and path of every tracked fixture, so the digest changes if any fixture's CONTENT, mode or name
changes -- it is not a timestamp check. The reported runner records it before and after the run
and asserts equality, which turns "test/test_images/ is read-only input" from a stated intention
into a checked property. Uses only `SHA` (a Julia standard library) and `Base.read` on a `Cmd`;
no dependency is added anywhere in this arm.
"""
function verify_real_readonly_digest()
    listing = read(Cmd(`git ls-files -s test/test_images`; dir = P13_REPO_ROOT), String)
    isempty(strip(listing)) && throw(ErrorException(
        "verify_real_readonly_digest: git listed no tracked files under test/test_images. The " *
        "six microscopy TIFFs are committed fixtures; an empty listing means this is not the " *
        "shipped tree, or git is unavailable."))
    return bytes2hex(SHA.sha256(listing))
end
