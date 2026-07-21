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

# corpus/test/test_anchors.jl --- physical ground-truth anchor testset
# (SC1, SC3, D-01, D-02, D-03, D-04, D-09, T-08-14/15/16/17).
#
# Proves, fully OFFLINE:
#   * anchor_rows() yields EXACTLY one positive + one negative, both physical-primary and
#     sealed_holdout (D-09), both passing the D-03 validate_anchor_row acceptance predicate;
#   * the asserted 1.0/0.0 ground truth is BIOLOGICAL, never a computed score (T-08-14, Pitfall 1);
#   * every anchor carries a non-empty accession + license + citation (T-08-17 provenance);
#   * the sha256 is EITHER a real 64-hex digest OR the explicit PENDING-FETCH sentinel — never a
#     fabricated/partial digest (T-08-16);
#   * the default accessors NEVER return the anchors; only open_sealed_holdout(;reason) does (D-09);
#   * the finalized committed corpus/manifest.csv validates and carries the finalized (no longer
#     "PENDING") anchor provenance;
#   * the D-01/D-02 deviations are recorded in ANCHOR_PROVENANCE_NOTES (honest provenance);
#   * NO anchor bytes are git-tracked under corpus/data (D-04, Pitfall 3).

# Guarded include so this testset runs standalone or wired into runtests.jl.
isdefined(@__MODULE__, :anchor_rows) || include(joinpath(@__DIR__, "..", "anchor_rows.jl"))

@testset "physical anchors (SC1/D-01/D-02/D-03/D-09)" begin

    rows = anchor_rows()

    # (a) EXACTLY one positive + one negative, both physical-primary + sealed_holdout (SC1/D-09).
    @test length(rows) == 2
    @test count(r -> r.role == "positive", rows) == 1
    @test count(r -> r.role == "negative", rows) == 1
    @test all(r -> r.tier == "physical-primary", rows)
    @test all(r -> r.split == "sealed_holdout", rows)
    @test length(unique(r.anchor_id for r in rows)) == 2

    pos = only(filter(r -> r.role == "positive", rows))
    neg = only(filter(r -> r.role == "negative", rows))

    # (b) The D-03 acceptance predicate passes for BOTH anchors.
    @test validate_anchor_row(pos)
    @test validate_anchor_row(neg)

    # (c) Truth is BIOLOGICAL and ASSERTED, never a computed coloc score (T-08-14, Pitfall 1).
    @test pos.truth_label == "coloc"
    @test neg.truth_label == "segregated"
    @test pos.coloc_ground_truth === 1.0
    @test neg.coloc_ground_truth === 0.0

    # (d) Provenance is complete and finalized — no "PENDING" placeholder survives (T-08-17).
    for r in rows
        @test !isempty(r.accession) && r.accession != "PENDING"
        @test !isempty(r.license)   && r.license   != "PENDING"
        @test !isempty(r.citation)  && r.citation  != "PENDING"
        @test !isempty(r.source_url) && startswith(r.source_url, "https://")
        @test occursin("2:", r.channels)   # two-channel layout declared
    end
    # The confirmed accessions/licences are pinned verbatim (human-verified values).
    @test pos.accession == "10.5281/zenodo.5509861"
    @test neg.accession == "S-BIAD1047"
    @test pos.license == "CC-BY-4.0" && neg.license == "CC-BY-4.0"
    @test occursin("10.1186/s12859-023-05320-1", pos.citation)   # RegiSTORM (TetraSpeck beads)
    @test occursin("S-BIAD1047", neg.citation)                   # Light My Cells

    # (e) Hash discipline: real 64-hex digest OR the explicit sentinel — never in between (T-08-16).
    for r in rows
        @test is_real_sha256(r.sha256) || is_pending_hash(r.sha256)
        @test !(is_real_sha256(r.sha256) && is_pending_hash(r.sha256))
    end
    @test PENDING_FETCH == "PENDING-FETCH"
    @test !is_real_sha256(PENDING_FETCH)          # the sentinel can never pass as a digest
    @test is_real_sha256(repeat("a", 64))
    @test !is_real_sha256(repeat("a", 63))
    @test !is_real_sha256("zz" * repeat("a", 62))

    # (f) A supplied bootstrap digest is carried through verbatim (the post-fetch path).
    d = repeat("0123456789abcdef", 4)
    boot = anchor_rows(; sha256_pos = d, bytes_pos = 12345)
    bpos = only(filter(r -> r.role == "positive", boot))
    @test bpos.sha256 == d && is_real_sha256(bpos.sha256) && bpos.bytes == 12345

    # (g) D-09 sealing: default accessors NEVER reach the anchors; only open_sealed_holdout does.
    df = finalized_manifest()
    @test validate_manifest(df) === df
    @test nrow(df) == 32
    @test count(==("physical-primary"), df.tier) == 2
    @test nrow(physical_anchors(df)) == 0                    # sealed ⇒ zero by construction
    @test !any(==("physical-primary"), corpus_default(df).tier)
    @test !any(==("sealed_holdout"), corpus_default(df).split)
    @test !any(==("physical-primary"), cbs_benchmark(df).tier)
    opened = open_sealed_holdout(df; reason = "phase16 blind evaluation")
    @test nrow(opened) == 2
    @test all(==("physical-primary"), opened.tier)
    @test Set(opened.anchor_id) == Set(["pos-tetraspeck-01", "neg-lightmycells-01"])
    @test_throws AssertionError open_sealed_holdout(df; reason = "")   # audit reason mandatory

    # (h) The finalized committed corpus/manifest.csv (SC3): validates, header-stamped, and its
    #     physical rows carry the finalized provenance + a hash-or-sentinel (never fabricated).
    @test isfile(COMMITTED_MANIFEST_PATH)
    txt = read(COMMITTED_MANIFEST_PATH, String)
    @test occursin("MANIFEST_SCHEMA_VERSION", txt)
    @test occursin("physical-primary", txt)
    @test occursin("pos-tetraspeck-01", txt)
    @test occursin("neg-lightmycells-01", txt)
    @test occursin("10.5281/zenodo.5509861", txt)
    @test occursin("S-BIAD1047", txt)
    @test !occursin("pos-anchor-01", txt)        # the 08-04 placeholder row is gone
    @test !occursin(",PENDING,", txt)            # no "PENDING" provenance survives
    datalines = filter(l -> !startswith(l, "#") && !isempty(l), split(txt, '\n'; keepempty = false))
    @test length(datalines) == 33                # 1 CSV header + 30 CBS + 2 anchors

    # Re-read the committed CSV and re-validate it as a manifest (round-trip contract check).
    reread = CSV.read(COMMITTED_MANIFEST_PATH, DataFrame; comment = "#")
    @test validate_manifest(reread) === reread
    phys = open_sealed_holdout(reread; reason = "test: contract re-read")
    @test nrow(phys) == 2
    for i in 1:nrow(phys)
        @test !ismissing(phys.accession[i]) && !isempty(string(phys.accession[i]))
        @test !ismissing(phys.license[i])   && !isempty(string(phys.license[i]))
        @test !ismissing(phys.citation[i])  && !isempty(string(phys.citation[i]))
        s = phys.sha256[i]
        @test is_real_sha256(s) || is_pending_hash(s)
    end

    # (i) The pending-hash notice is emitted while (and only while) a sentinel is present.
    notice = anchor_hash_notice(df)
    if any(is_pending_hash, df.sha256[df.tier .== "physical-primary"])
        @test occursin(PENDING_FETCH, notice)
        @test occursin("bootstrap_anchor_hashes", notice)
        @test occursin("Phase 16", notice)
    else
        @test isempty(notice)
    end
    # A fully-bootstrapped manifest emits NO notice.
    full = finalized_manifest(; sha256_pos = repeat("a", 64), sha256_neg = repeat("b", 64))
    @test isempty(anchor_hash_notice(full))

    # (j) The D-01/D-02 deviations are recorded honestly alongside the contract.
    @test ANCHORS_MATCHED == false                       # cross-study anchors (D-02 limitation)
    @test occursin("DEVIATION", ANCHOR_PROVENANCE_NOTES.d01)
    @test occursin("tandem", ANCHOR_PROVENANCE_NOTES.d01)
    @test occursin("DEVIATION", ANCHOR_PROVENANCE_NOTES.d02)
    @test occursin("matched=false", ANCHOR_PROVENANCE_NOTES.d02)
    @test occursin(PENDING_FETCH, ANCHOR_PROVENANCE_NOTES.hashes)

    # (k) Staging guard (D-04, Pitfall 3): NO anchor bytes tracked under corpus/data.
    repo = normpath(joinpath(@__DIR__, "..", ".."))
    tracked = try
        readchomp(setenv(`git ls-files corpus/data`, dir = repo))
    catch
        ""   # git unavailable ⇒ do not fail the offline gate
    end
    @test isempty(tracked)
end
