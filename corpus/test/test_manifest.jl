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

# corpus/test/test_manifest.jl --- manifest data-contract testset (SC2, SC3, D-03, D-06, D-08, D-09).
#
# Proves the data contract is enforced BY CONSTRUCTION, mirroring the Phase-3
# "no standardize_all symbol exists" leakage-impossible discipline:
#   * schema is validated (columns / tier / split / role enums, D-09 physical=>sealed invariant);
#   * the default corpus iterator NEVER yields a sealed_holdout row (D-09 anti-snooping);
#   * physical anchors are reachable ONLY via open_sealed_holdout(;reason) (D-09);
#   * tiers are NEVER silently blended -- mixed_tier is the sole both-tier path (SC2/D-06);
#   * an anchor row is rejected unless it carries a biological truth_label + open license
#     + accession + an asserted (not computed) 1.0/0.0 ground truth (D-03, Pitfall 1);
#   * the committed manifest is content-addressed so any contract edit is detectable (D-08).
# Runs fully OFFLINE.

# Guarded include so this testset runs standalone or wired into runtests.jl.
isdefined(@__MODULE__, :build_manifest) || include(joinpath(@__DIR__, "..", "manifest.jl"))

# --- Row helper: a full 14-field manifest row as a NamedTuple in MANIFEST_COLUMNS order. -------
function _mrow(; anchor_id, tier, role, truth_label, coloc_ground_truth, accession,
                source_url, citation, license, format, channels, sha256, split, bytes)
    return (; anchor_id, tier, role, truth_label, coloc_ground_truth, accession,
             source_url, citation, license, format, channels, sha256, split, bytes)
end

# Canonical fixture rows: two SEALED physical anchors + two CBS benchmark rows.
_pos_row() = _mrow(; anchor_id = "pos-tandem-01", tier = "physical-primary", role = "positive",
    truth_label = "coloc", coloc_ground_truth = 1.0, accession = "10.0001/pos",
    source_url = "https://example.org/pos.tif", citation = "Doe et al. 2020",
    license = "CC-BY-4.0", format = "tiff", channels = "2:ch1,ch2",
    sha256 = repeat("a", 64), split = "sealed_holdout", bytes = 1234)

_neg_row() = _mrow(; anchor_id = "neg-nucmem-01", tier = "physical-primary", role = "negative",
    truth_label = "segregated", coloc_ground_truth = 0.0, accession = "10.0001/neg",
    source_url = "https://example.org/neg.tif", citation = "Doe et al. 2020",
    license = "CC-BY-4.0", format = "tiff", channels = "2:ch1,ch2",
    sha256 = repeat("b", 64), split = "sealed_holdout", bytes = 5678)

_cbs_row(id, gt, split) = _mrow(; anchor_id = id, tier = "simulated-secondary", role = "benchmark",
    truth_label = "simulated-degree", coloc_ground_truth = gt, accession = "CBS-$(id)",
    source_url = "https://cbs.example.org/$(id).tif", citation = "CBS benchmark",
    license = "CC0", format = "tiff", channels = "2:red,green",
    sha256 = repeat("c", 64), split = split, bytes = 9012)

_valid_rows() = [_pos_row(), _neg_row(), _cbs_row("cbs-050", 0.5, "dev"), _cbs_row("cbs-080", 0.8, "eval")]

@testset "manifest contract (SC2/SC3/D-03/D-06/D-08/D-09)" begin

    # (a) build_manifest yields a DataFrame with EXACTLY the MANIFEST_COLUMNS, in order.
    df = build_manifest(_valid_rows())
    @test propertynames(df) == collect(MANIFEST_COLUMNS)
    @test nrow(df) == 4

    # (b) A schema-valid manifest passes validate_manifest.
    @test validate_manifest(df) === df

    # (c) A bad tier / split / role each throw.
    bad_tier = build_manifest([_cbs_row("x", 0.5, "dev")])
    bad_tier.tier[1] = "bogus-tier"
    @test_throws ErrorException validate_manifest(bad_tier)

    bad_split = build_manifest([_cbs_row("x", 0.5, "dev")])
    bad_split.split[1] = "bogus-split"
    @test_throws ErrorException validate_manifest(bad_split)

    bad_role = build_manifest([_cbs_row("x", 0.5, "dev")])
    bad_role.role[1] = "bogus-role"
    @test_throws ErrorException validate_manifest(bad_role)

    # (d) D-09 invariant: a physical-primary row whose split != sealed_holdout throws.
    leaky = build_manifest([_pos_row()])
    leaky.split[1] = "dev"
    @test_throws ErrorException validate_manifest(leaky)

    # (e) A missing required column throws.
    @test_throws ErrorException validate_manifest(select(df, Not(:bytes)))

    # (f) corpus_default NEVER yields a sealed_holdout row (D-09 anti-snooping).
    d0 = corpus_default(df)
    @test nrow(d0) == 2
    @test !any(==("sealed_holdout"), d0.split)

    # (g) physical_anchors returns NO rows by default (all physical anchors are sealed).
    @test nrow(physical_anchors(df)) == 0

    # (h) cbs_benchmark returns ONLY simulated-secondary rows.
    cb = cbs_benchmark(df)
    @test nrow(cb) == 2
    @test all(==("simulated-secondary"), cb.tier)

    # (i) open_sealed_holdout is the SOLE path to the sealed anchors and requires a reason.
    @test_throws AssertionError open_sealed_holdout(df; reason = "")
    sealed = open_sealed_holdout(df; reason = "phase16-evaluation")
    @test nrow(sealed) == 2
    @test all(==("physical-primary"), sealed.tier)
    @test all(==("sealed_holdout"), sealed.split)

    # (j) mixed_tier is the ONLY function returning both tiers together.
    mt = mixed_tier(df)
    @test nrow(mt) == 4
    @test Set(mt.tier) == Set(("physical-primary", "simulated-secondary"))

    # (k) validate_anchor_row (D-03 acceptance predicate).
    @test validate_anchor_row(_pos_row()) === true          # biology label + open license + accession + 1.0
    @test validate_anchor_row(_neg_row()) === true          # segregated + 0.0
    # reject: missing license.
    @test validate_anchor_row(merge(_pos_row(), (; license = ""))) === false
    # reject: missing accession.
    @test validate_anchor_row(merge(_pos_row(), (; accession = ""))) === false
    # reject: a computed (non-biological) truth label.
    @test validate_anchor_row(merge(_pos_row(), (; truth_label = "0.73"))) === false
    # reject: a computed (non-asserted) ground-truth score on an anchor (Pitfall 1).
    @test validate_anchor_row(merge(_pos_row(), (; coloc_ground_truth = 0.73))) === false

    # (l) write_manifest content-addresses into a hash-named dir with a pre-registration
    #     CSV header (schema version + CORPUS_MASTER_SEED) plus a JLD2 twin, atomically.
    outdir = mktempdir()
    dir = write_manifest(df; outdir = outdir)
    @test isdir(dir)
    csv = joinpath(dir, "manifest.csv")
    jld = joinpath(dir, "manifest.jld2")
    @test isfile(csv) && isfile(jld)
    @test !isfile(csv * ".tmp") && !isfile(jld * ".tmp")     # atomic: no .tmp left behind
    txt = read(csv, String)
    @test occursin("MANIFEST_SCHEMA_VERSION", txt)
    @test occursin(string(MANIFEST_SCHEMA_VERSION), txt)
    @test occursin("CORPUS_MASTER_SEED", txt)
    @test occursin(repr(CORPUS_MASTER_SEED), txt)
    reloaded = JLD2.load(jld, "manifest")
    @test nrow(reloaded) == 4
    # Deterministic content address: same (df, sources) => same dir.
    @test write_manifest(df; outdir = outdir) == dir
end
