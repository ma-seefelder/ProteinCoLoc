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

# corpus/config.jl --- pre-declared corpus data-contract constants (D-04, D-08).
#
# Every enum, seed, schema tag, and threshold the external ground-truth corpus is
# built against is committed HERE, in this Wave-0 plan (08-01), BEFORE any manifest,
# any download, or any split assignment exists. This is the anti-data-snooping
# guarantee (D-08): the tier/split/role vocabulary and the master seed are FIXED up
# front, NOT tuned after seeing which datasets land where. Every downstream plan
# (hash, fetch, manifest, load, CBS, anchors) reads its vocabulary from this single
# source of truth so the contract cannot drift silently.
#
# DECOUPLING (hard constraint, CLAUDE.md): corpus-local declarative config. This file
# touches no src/, no root Project.toml. Pure `const` declarations -- no logic beyond
# a `normpath`/`joinpath` path derivation.

# --- Corpus master seed (D-08) --------------------------------------------------
# Random123 (Philox) master seed for every seeded corpus draw (CBS subset selection,
# dev/eval split assignment), so the corpus is bit-reproducible across runs and thread
# counts, exactly like the Phase-3 cache and the Phase-9 comparator harness.
# DISTINCT from every other master seed in the repo so the corpus draw streams can
# never collide with the spike validation / NPE / comparator streams:
#   spike VAL_MASTER_SEED = 0x5BC0FFEE, NPE_MASTER_SEED = 0xC0FFEE,
#   comparator MASTER_SEED = 0x00C0FFEE.
const CORPUS_MASTER_SEED = 0x0000000000C05EED  # corpus Philox master seed (distinct namespace)

# --- Manifest schema version (D-08) ---------------------------------------------
# Bump when MANIFEST_COLUMNS changes shape/meaning so a stale manifest is detectable.
const MANIFEST_SCHEMA_VERSION = 1

# --- Tier / split / role vocabulary (D-06, D-09) --------------------------------
# The physical anchors are `physical-primary`; the CBS benchmark is `simulated-secondary`.
# Tiers are structurally separated by the manifest accessors (D-06), never merely by a
# comment -- this enum is the single source that guard code validates against (SC2).
const TIERS = ("physical-primary", "simulated-secondary")

# `sealed_holdout` is the blind physical-anchor split opened ONLY at Phase-16 evaluation
# (D-09); `dev`/`eval` are the seeded CBS working splits. The sealed split is unreachable
# from every default manifest accessor by construction.
const SPLITS = ("sealed_holdout", "dev", "eval")

# Dataset roles within a tier: a 100%-coloc positive anchor, a segregated negative anchor,
# and the graded CBS benchmark rows.
const ROLES = ("positive", "negative", "benchmark")

# --- Manifest schema columns (D-08) ---------------------------------------------
# The 14-field data contract. Downstream manifest assembly maps one row per anchor/image
# into exactly these columns, in this order, so the committed manifest.csv is a stable,
# reviewable contract artifact.
const MANIFEST_COLUMNS = (
    :anchor_id,           # stable within-corpus identifier
    :tier,                # ∈ TIERS (physical-primary | simulated-secondary)
    :role,                # ∈ ROLES (positive | negative | benchmark)
    :truth_label,         # human-readable biological truth (e.g. "colocalized")
    :coloc_ground_truth,  # numeric truth degree asserted by biology/citation, NOT computed
    :accession,           # DOI / archive accession pinning provenance
    :source_url,          # retrievable byte source
    :citation,            # publication asserting the truth (non-circularity)
    :license,             # redistributable license tag
    :format,              # on-disk format (e.g. tiff)
    :channels,            # channel count / layout
    :sha256,              # SHA-256 content hash over data-defining bytes (integrity)
    :split,               # ∈ SPLITS (sealed_holdout | dev | eval)
    :bytes,               # byte size of the source artifact
)

# --- Disjoint split-stream salt (D-09) ------------------------------------------
# XOR-ed into the first Philox key word to carve a provably disjoint draw stream for the
# CBS dev/eval split, mirroring the datagen.jl HOLDOUT_SALT/FOLD_SALT idiom. Disjointness
# from every other stream holds because CORPUS_MASTER_SEED itself is a distinct namespace;
# the salt additionally separates the split stream from any future corpus draw stream.
const CBS_SPLIT_SALT = 0xD1B54A32D192ED03  # CBS dev/eval split Philox key salt

# --- Fetch / integrity constants (D-04, D-05) -----------------------------------
# Wall-clock ceiling for a single source download; a stalled fetch is killed rather than
# hanging CI (ASVS DoS mitigation). Applied by the 08-04 fetch layer.
const DOWNLOAD_TIMEOUT_S = 120

# NOTE (D-04/D-05): the download-integrity hash MUST be cryptographic SHA-256
# (`SHA.sha256` / `SHA256_CTX`), NEVER Base `hash`. A content-hash MISMATCH is a HARD
# error (tamper/corruption -> abort); an UNREACHABLE source is a graceful skip-with-flag.
# The two paths are deliberately distinct and must never share one try/catch.

# --- Summary-path validity floor (Pitfall 6) ------------------------------------
# `src/colocalization.jl:correlation` drops zero/NaN/missing pixels via `_exclude_zero`
# and marks any 8x8 patch with `length <= 15` valid pixels as `missing`. Documented here
# so fixture/anchor curation can guarantee dense-enough non-zero signal per patch.
const MIN_PATCH_VALID_PIXELS = 15

# --- Downloaded-bytes directory (D-04) ------------------------------------------
# Where fetched image bytes land. This tree is gitignored (see .gitignore corpus/data/*):
# large/licensed bytes never enter git; only the manifest/fetch-script/tests/fixture are
# tracked. `normpath`/`@__DIR__` keeps the path stable regardless of caller cwd.
const CORPUS_DATA_DIR = normpath(joinpath(@__DIR__, "data"))
