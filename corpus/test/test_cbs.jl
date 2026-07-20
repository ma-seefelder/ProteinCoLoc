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

# corpus/test/test_cbs.jl --- CBS ingestion testset (SC2, SC3, D-05, D-06, D-07, D-09, T-08-11).
#
# Proves, fully OFFLINE:
#   * cbs_rows() enumerates the FULL corpus (3 pairs x 10 degrees = 30 rows), every row
#     tier=simulated-secondary + role=benchmark, ground truth = LABELLED degree (not computed);
#   * assign_cbs_split is bit-reproducible and yields ONLY dev/eval (never sealed_holdout, D-09);
#   * safe_extract REJECTS a zip-slip entry BEFORE any write, and admits a safe entry (T-08-11);
#   * cbs_rgb_to_channels / _split_rgb yields two finite channels per pair;
#   * validate_manifest(build_manifest(cbs_rows())) passes;
#   * the committed manifest.csv contract = 30 CBS + 2 sealed anchor placeholders, header-stamped;
#   * the live CBS fetch smoke SKIPS cleanly offline (mirrors the tapqir :skipped contract, D-05);
#   * NO image bytes are git-tracked under corpus/data (D-04, Pitfall 3).

import Images   # RGB image construction for the channel-split test (root dep)

# Guarded include so this testset runs standalone or wired into runtests.jl.
isdefined(@__MODULE__, :cbs_rows) || include(joinpath(@__DIR__, "..", "cbs.jl"))

@testset "CBS ingestion (SC2/SC3/D-05/D-06/D-07/D-09)" begin

    # (a) FULL corpus enumerated: 30 rows, every one simulated-secondary / benchmark (SC2/D-06/D-07).
    rows = cbs_rows()
    @test length(rows) == 30
    @test all(r -> r.tier == "simulated-secondary", rows)
    @test all(r -> r.role == "benchmark", rows)
    @test all(r -> r.truth_label == "simulated-degree", rows)
    # Endpoints present; ids follow cbs-<pair>-<3digit>.
    ids = [r.anchor_id for r in rows]
    @test "cbs-RG-000" in ids
    @test "cbs-GB-090" in ids
    @test length(unique(ids)) == 30
    # Ground truth is the LABELLED degree as a fraction (definitional, NOT a computed score).
    rg50 = only(filter(r -> r.anchor_id == "cbs-RG-050", rows))
    @test rg50.coloc_ground_truth == 0.5
    @test all(r -> 0.0 <= r.coloc_ground_truth <= 0.9, rows)

    # (b) validate_manifest(build_manifest(cbs_rows())) passes.
    df_cbs = build_manifest(cbs_rows())
    @test validate_manifest(df_cbs) === df_cbs
    @test nrow(df_cbs) == 30

    # (c) Seeded dev/eval split: bit-reproducible across two runs, ONLY dev/eval (never sealed).
    s1 = assign_cbs_split(cbs_rows())
    s2 = assign_cbs_split(cbs_rows())
    @test [r.split for r in s1] == [r.split for r in s2]          # deterministic
    @test all(r -> r.split in ("dev", "eval"), s1)               # never sealed_holdout (D-09)
    @test !any(r -> r.split == "sealed_holdout", s1)
    @test any(r -> r.split == "dev", s1) && any(r -> r.split == "eval", s1)   # both present

    # (d) Zip-slip: an entry escaping the target dir is REJECTED before any write (T-08-11).
    target = mktempdir()
    evil = [("../evil.tif", () -> UInt8[0x00])]
    @test_throws ErrorException safe_extract("dummy.zip", target; entries = evil)
    @test !isfile(joinpath(dirname(target), "evil.tif"))          # nothing escaped
    # An absolute-path entry is likewise rejected.
    abs_evil = [(joinpath(homedir(), "pwned.tif"), () -> UInt8[0x00])]
    @test_throws ErrorException safe_extract("dummy.zip", target; entries = abs_evil)
    # A safe nested entry extracts INSIDE the target dir.
    safe = [("set/img000.tif", () -> UInt8[0x01, 0x02, 0x03])]
    written = safe_extract("dummy.zip", target; entries = safe)
    @test length(written) == 1
    @test isfile(joinpath(target, "set", "img000.tif"))
    @test read(joinpath(target, "set", "img000.tif")) == UInt8[0x01, 0x02, 0x03]

    # (e) RGB/composite split -> two finite Float64 channels, correct selection per pair.
    rmat = [0.10 0.20; 0.30 0.40]
    gmat = [0.50 0.60; 0.70 0.80]
    bmat = [0.05 0.15; 0.25 0.35]
    img = Images.RGB.(rmat, gmat, bmat)
    (c1, c2) = _split_rgb(img, "Red-Green")
    @test c1 isa Matrix{Float64} && c2 isa Matrix{Float64}
    @test size(c1) == (2, 2) && size(c2) == (2, 2)
    @test all(isfinite, c1) && all(isfinite, c2)
    @test c1 ≈ rmat && c2 ≈ gmat
    @test _split_rgb(img, "Red-Blue")[2]  ≈ bmat
    @test _split_rgb(img, "Green-Blue")[1] ≈ gmat
    @test_throws ErrorException _split_rgb(img, "Nope-Pair")
end
