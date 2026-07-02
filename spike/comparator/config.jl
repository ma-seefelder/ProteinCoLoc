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

# spike/comparator/config.jl --- CMP-04/CMP-14: pre-declared thresholds (D-14).
#
# Every "knows when the classics are wrong" threshold the Phase-9 comparator scores
# against is committed HERE, in this Wave-0 plan (09-01), BEFORE any comparison table
# exists. This is the anti-data-snooping guarantee the roadmap flagged as a hard
# blocker (D-14): the divergence cutoffs and traffic-light bands are FIXED before the
# run, NOT tuned post-hoc to make a table look good. The artifact header (written by
# 09-05) quotes these consts verbatim as the pre-registration record.
#
# DECOUPLING (hard constraint, CLAUDE.md): spike-local declarative config; touches no
# src/, no root Project.toml. Pure `const` declarations — no logic beyond them.

# --- Costes randomization significance (D-05; RESEARCH Pattern 2) ---------------
# The Costes-p significance test (Costes et al. 2004, Biophys J 86:3993) scrambles the
# image into resolution-sized blocks, recomputes correlation on each scramble, and
# reports the fraction of scrambles whose correlation ≥ the observed one.
const COSTES_N_SCRAMBLE = 200            # scramble count — Costes et al. 2004 standard
const COSTES_BLOCK_PX   = 7              # block edge in px ≈ one resolution element for σ_psf = 1.3
const COSTES_SALT       = 0xA5A5A5A5DEADBEEF  # Philox key salt for the scramble stream;
                                              # DISTINCT from HOLDOUT_SALT (0x9E37…) and
                                              # FOLD_SALT (0xD1B5…) so the Costes null draws
                                              # never collide with the data/CV streams (D-10).
const COSTES_ALPHA      = 0.05           # significance cutoff for the colocalization verdict

# --- Divergence / traffic-light bands (D-10/D-14; RESEARCH Pattern 3) -----------
# The divergence indicator marks rows where a classical verdict diverges from the
# simulator ground-truth regime. Presented BayesInteractomics-style as a green/amber/red
# traffic light rather than a bare number. Bands are pre-declared here (D-14): amber at
# ≥ WARN, red at ≥ FAIL. WARN < FAIL is enforced by the testset.
const DIVERGENCE_WARN = 0.30             # amber band lower edge (moderate disagreement)
const DIVERGENCE_FAIL = 0.60             # red band lower edge (severe disagreement)

# --- Harness master seed (D-12) -------------------------------------------------
# Random123 (Philox) master seed for the whole comparator harness, so the emitted
# table is bit-reproducible across runs and thread counts, exactly like the Phase-3 cache.
const MASTER_SEED = 0x00000000_00C0FFEE  # comparator harness Philox master seed

# --- Tapqir sanity-anchor value (D-07/D-08; resolved by 09-04) ------------------
# The Tapqir bridge reproduces a PUBLISHED Tapqir tutorial dataset (eLife 2022 11:e73860,
# Part II `cosmos` sample data set) and checks the recovered quantity against the published
# value (a sanity anchor only, NOT a comparator column on our simulator images — D-07). The
# published reference value is captured ONCE by running the tutorial in an isolated Conda
# env; the capture records published ground truth and is NOT a tuned threshold (D-14).
#
# CLEAN-SKIP PATH TAKEN (09-04, this machine — Windows). No Conda backend is available
# (conda/mamba/micromamba absent from PATH; JULIA_CONDAPKG_BACKEND unset), and materializing
# tapqir 1.1.19 (2023) would require a large, likely-failing network install (the system
# Python is 3.14, unsupported by that Tapqir release). Per D-08 the anchor is optional and
# must never gate the classical battery, so we DELIBERATELY leave TAPQIR_PUBLISHED_VALUE as
# NaN. This is the documented clean-skip, NOT a silently guessed number (D-14 discipline):
# with NaN the bridge cannot assess an anchor and returns status=:skipped, which the 09-06
# D-13 gate accepts (status ∈ {:passed, :skipped}). To pin a real value later, materialize
# spike/comparator/tapqir_env/ on a machine with Conda and re-run 09-04 Task 2.
const TAPQIR_PUBLISHED_VALUE = NaN       # documented clean-skip (no Conda backend); pin a real value by re-running 09-04 Task 2 with a materialized env
const TAPQIR_TOL             = 0.05      # provisional acceptance tolerance for the anchor when a value is pinned
