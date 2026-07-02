<#
.SYNOPSIS
    Build + REGRESSION-gate the ProteinCoLoc v2.0 (AmortizedColoc) manuscript with the bundled
    Libertinus fonts.

.DESCRIPTION
    Primary Windows regression gate (D-12, D-17). Compiles main.typ to main.pdf, passing the
    project-local fonts/ directory so Libertinus Serif / Sans / Math are picked up regardless of
    what is installed system-wide, and setting --root to the manuscript directory so root-absolute
    lookups of /colours.yaml and /figures/*.yaml (the per-figure YAML convention, D-14/D-15)
    resolve. This keeps the build reproducible.

    Beyond compiling, this gate ASSERTS the Wave-2/3 content is present so later phases
    regression-check the skeleton (D-12): the claim table renders (SC1); the related-work matrix +
    four Tapqir axes are present (SC2); >=7 figure specs enumerate their producing phases (SC3);
    every figure-spec 'Supports: Cx' token resolves to a claim id (D-07); and the per-figure YAML
    machinery holds — lib/figstyle.typ exists, f1_speedup.yaml exists, and f1_speedup.typ has NO
    hard-coded geometry/hex (D-14/D-16). Assertion comparisons filter comment lines so a token in a
    header comment cannot self-satisfy (or falsely trip) a check.

    Typst version: the template this manuscript is founded on (BayesInteractomics_Method_paper)
    pinned Typst 0.14; this machine runs the installed 0.15.0. We default to the installed
    `typst` on PATH and keep the $env:TYPST override hook (D-17).

    Exit-code correctness (D-17): the compile line has NO pipe, so `typst`'s own exit status is
    preserved. A broken build cannot silently pass. A non-fatal "did not converge / did not
    stabilize" warning from lib/tiered-bib.typ is EXPECTED and does NOT fail the gate — we read
    the compiler's exit code, we do NOT treat stderr warnings as failures.

.PARAMETER Watch
    Recompile automatically on file changes (typst watch). Content assertions are skipped in
    watch mode (the compiler runs indefinitely).

.PARAMETER Source
    Source file to compile (default: main.typ).

.PARAMETER Output
    Output PDF path (default: main.pdf).

.EXAMPLE
    ./build.ps1
    Compile once and run the regression assertions.

.EXAMPLE
    ./build.ps1 -Watch
    Watch and recompile on changes (no assertions).
#>
[CmdletBinding()]
param(
    [switch]$Watch,
    [string]$Source = "main.typ",
    [string]$Output = "main.pdf"
)

$ErrorActionPreference = "Stop"
$root = $PSScriptRoot
$fontPath = Join-Path $root "fonts"
$typst = if ($env:TYPST) { $env:TYPST } else { "typst" }

function Fail($msg) { Write-Host "GATE FAIL: $msg" -ForegroundColor Red; exit 1 }

if (-not (Get-Command $typst -ErrorAction SilentlyContinue)) {
    Fail "Typst binary '$typst' not found on PATH (default: installed 'typst' 0.15.0). Override with the TYPST env var."
}

if (-not (Test-Path $fontPath)) {
    Fail "Font directory not found: $fontPath (the bundled Libertinus OTFs are required)."
}

$src = Join-Path $root $Source
$out = Join-Path $root $Output

$cmd = if ($Watch) { "watch" } else { "compile" }
Write-Host "$typst $cmd --root `"$root`" --font-path `"$fontPath`" `"$src`" `"$out`"" -ForegroundColor Cyan
# NO pipe on the compile line — preserve typst's exit code (D-17, Pitfall 1).
& $typst $cmd --root $root --font-path $fontPath $src $out
if ($Watch) { exit $LASTEXITCODE }
if ($LASTEXITCODE -ne 0) { Fail "typst compile exited $LASTEXITCODE" }

# --- Content assertions (each exits non-zero on failure) --------------------------------------
# Read a file with full-line Typst comments stripped (// only — in Typst `#` is CODE, e.g.
# `#let`, `#table(`, not a comment) so a token that appears only inside a header comment cannot
# satisfy (or, for the negative F1 checks, falsely trip) an assertion.
function Get-NonComment($path) {
    Get-Content $path | Where-Object { $_ -notmatch '^\s*//' }
}

$claimsPath   = Join-Path $root "claims.typ"
$rwPath       = Join-Path $root "sections/related_work.typ"
$specsPath    = Join-Path $root "figures/specs.typ"
$figstylePath = Join-Path $root "lib/figstyle.typ"
$f1TypPath    = Join-Path $root "figures/f1_speedup.typ"
$f1YamlPath   = Join-Path $root "figures/f1_speedup.yaml"

# SC1 — compiles + the claim table renders.
if (-not (Test-Path $out)) { Fail "SC1: $out was not produced by the compile step" }
if (-not ((Get-NonComment $claimsPath) -match 'claim_table')) { Fail "SC1: claim_table render function missing from claims.typ" }

# SC2 — related-work matrix + all four Tapqir differentiation axes.
$rw = Get-NonComment $rwPath
if (-not ($rw -match 'amortiz'))      { Fail "SC2: 'amortization' axis missing from related_work.typ" }
if (-not ($rw -match 'sbc|calibrat')) { Fail "SC2: 'SBC/calibration' axis missing from related_work.typ" }
if (-not ($rw -match 'registration')) { Fail "SC2: 'registration' axis missing from related_work.typ" }
if (-not ($rw -match 'spatial'))      { Fail "SC2: 'spatial' axis missing from related_work.typ" }
if (-not ($rw -match 'table\('))      { Fail "SC2: related-work comparison matrix (table(...)) missing" }

# SC3 — >=7 figure specs enumerating their producing phases.
$specNonComment = Get-NonComment $specsPath
$specCount = ($specNonComment | Select-String -Pattern 'SPEC F[0-9]').Count
if ($specCount -lt 7) { Fail "SC3: only $specCount figure specs found (<7)" }

# D-07 — every 'Supports: Cx' token in specs.typ resolves to a claim id in claims.typ.
$claimsRaw = (Get-Content $claimsPath -Raw)
$specRaw   = ($specNonComment -join "`n")
$ids = New-Object System.Collections.Generic.HashSet[string]
foreach ($m in [regex]::Matches($specRaw, 'Supports:\s*([C0-9,\s]+)')) {
    foreach ($c in [regex]::Matches($m.Groups[1].Value, 'C[0-9]+')) { [void]$ids.Add($c.Value) }
}
foreach ($cid in $ids) {
    if ($claimsRaw -notmatch ('id:\s*"' + $cid + '"')) {
        Fail "D-07: figure spec references $cid but no matching 'id: `"$cid`"' in claims.typ"
    }
}

# D-14/D-16 — the per-figure YAML machinery holds.
if (-not (Test-Path $figstylePath)) { Fail "D-14: lib/figstyle.typ missing" }
if (-not (Select-String -Path $figstylePath -Pattern 'figstyle' -Quiet)) { Fail "D-14: figstyle loader not defined in lib/figstyle.typ" }
if (-not (Test-Path $f1YamlPath)) { Fail "D-14: figures/f1_speedup.yaml missing" }
if (-not (Select-String -Path $f1TypPath -Pattern 'figstyle\("f1_speedup"\)' -Quiet)) { Fail "D-16: f1_speedup.typ does not load its YAML via figstyle(`"f1_speedup`")" }
# F1 must be fully YAML-driven: NO hard-coded graphical length in ANY unit (cm/mm/in/pt/em) and
# NO inline hex colour. The ONLY sanctioned way a length literal appears is a unit-conversion idiom
# applied to a YAML value: `* 10mm` (via _cm), `* 1pt`, or `* 1em`. Strip those, then flag any
# REMAINING numeric length literal — so a future `width: 5mm` / `width: 200pt` (units the old
# cm-only check missed, WR-04) trips the gate. Nothing graphical may be hard-coded (D-14).
$f1 = Get-NonComment $f1TypPath
$f1stripped = $f1 -replace '\*\s*10mm', '' -replace '\*\s*1pt', '' -replace '\*\s*1em', ''
if ($f1stripped -match '[0-9]+(\.[0-9]+)?(cm|mm|in|pt|em)') { Fail "D-14: f1_speedup.typ hard-codes a graphical length (cm/mm/in/pt/em) outside the sanctioned '* 10mm' / '* 1pt' / '* 1em' idioms — geometry must come from the YAML" }
# Positive assertion: the panel width/height are explicitly YAML-driven (fs.panel.* via _cm).
if (-not ($f1 -match 'width:\s*_cm\(fs\.panel\.'))  { Fail "D-14: f1_speedup.typ panel width is not YAML-driven (expected _cm(fs.panel.width))" }
if (-not ($f1 -match 'height:\s*_cm\(fs\.panel\.')) { Fail "D-14: f1_speedup.typ panel height is not YAML-driven (expected _cm(fs.panel.height))" }
if ($f1 -match 'rgb\("#')            { Fail "D-14: f1_speedup.typ contains an inline hex colour (must come from the palette role)" }

Write-Host "OK: DoD met — compiles + claim table + matrix + >=7 specs ($specCount) + YAML-driven F1" -ForegroundColor Green
exit 0
