<#
.SYNOPSIS
    Build the ProteinCoLoc v2.0 (AmortizedColoc) manuscript with the bundled Libertinus fonts.

.DESCRIPTION
    Primary Windows compile gate. Compiles main.typ to main.pdf, passing the project-local
    fonts/ directory so Libertinus Serif / Sans / Math are picked up regardless of what is
    installed system-wide, and setting --root to the manuscript directory so root-absolute
    lookups of /colours.yaml and /figures/*.yaml (the per-figure YAML convention, D-14/D-15)
    resolve. This keeps the build reproducible.

    Typst version: the template this manuscript is founded on (BayesInteractomics_Method_paper)
    pinned Typst 0.14; this machine runs the installed 0.15.0. We default to the installed
    `typst` on PATH and keep the $env:TYPST override hook (D-17).

    Exit-code correctness (D-17): the compile line has NO pipe, so `typst`'s own exit status
    is preserved and surfaced via `exit $LASTEXITCODE`. A broken build cannot silently pass.

.PARAMETER Watch
    Recompile automatically on file changes (typst watch).

.PARAMETER Source
    Source file to compile (default: main.typ).

.PARAMETER Output
    Output PDF path (default: main.pdf).

.EXAMPLE
    ./build.ps1
    Compile once.

.EXAMPLE
    ./build.ps1 -Watch
    Watch and recompile on changes.
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

if (-not (Get-Command $typst -ErrorAction SilentlyContinue)) {
    Write-Error "Typst binary '$typst' not found on PATH (default: installed 'typst' 0.15.0). Override with the TYPST env var."
    exit 1
}

if (-not (Test-Path $fontPath)) {
    Write-Error "Font directory not found: $fontPath (the bundled Libertinus OTFs are required)."
    exit 1
}

$src = Join-Path $root $Source
$out = Join-Path $root $Output

$cmd = if ($Watch) { "watch" } else { "compile" }
Write-Host "$typst $cmd --root `"$root`" --font-path `"$fontPath`" `"$src`" `"$out`"" -ForegroundColor Cyan
# NO pipe on the compile line — preserve typst's exit code (D-17, Pitfall 1).
& $typst $cmd --root $root --font-path $fontPath $src $out
exit $LASTEXITCODE
