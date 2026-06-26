# Detailplan — Teil 3: AmortizedColoc-Spike (amortisierte SBI für ProteinCoLoc)

> **Zusammenfassung.** Der de-riskte, sofort startbare Solo-On-Ramp. Ein **kleinster lauffähiger Spike** beweist das Kernprinzip: 2-Kanal-2D-Physik-Simulator + NPE auf der *vorhandenen* Summary-Statistik + 1 SBC-Check + 1 amortisierter Bayes-Faktor + OOD-Flag. Vollständig entkoppelt (eigenes `spike/`-Env, kein src-Eingriff), null Wet-lab/AF3, CPU genügt. **Julia-nativer Stack: NeuralEstimators.jl.** Realistisch ~4 Wochen bis Go/No-Go.

**Publizierbarer Claim:** Kalibrierte, amortisierte Bayes-Kolokalisation mit Simulation-Based-Calibration-Nachweis — Inferenz in Millisekunden statt ADVI-Minuten, mit nachgewiesener nomineller Abdeckung (erstes Coloc-Tool überhaupt mit SBC/Coverage + amortisiertem BF + OOD-Meldung).

## Scope-Guard (bewusst NICHT im Spike)
- Volle 4-Ebenen-Hierarchie (Pixel→Segment→Zelltyp→Nachbarschaft), 3D/Z-Stacks, chromatische Aberration höherer Ordnung, GUI, Multi→2-Kanal-Generalisierung, Produktiv-Refactor von `colocalization()`.
- **Strikt entkoppelt** von den zwei Manuskripten: alles in `ProteinCoLoc/spike/` mit **eigenem Project.toml**; kein Eingriff in `src/bayes.jl`, `src/colocalization.jl` oder die main-Pipeline.
- Kein Wet-lab, keine AF3/Boltz/Chai-Abhängigkeit, voll solo, CPU-only.

## Inputs (vorhandene Assets)
| Asset | Pfad | Liefert |
|---|---|---|
| `colocalization()` + Turing `@model` (hierarch. Student-t, ADVI) | `ProteinCoLoc/src/colocalization.jl` | Ground-Truth-Generativstruktur (Prior-Ränge μ/ν/σ/τ) + Referenz-Posterior-Baseline |
| `correlation(x,y;method)` + `patch()` + `_prepare_data()` | `ProteinCoLoc/src/colocalization.jl` | **fertige Summary-Statistik** (Patch-Korr-Vektor) — exakt der NPE-Input, kein Neuschreiben |
| `compute_BayesFactor()` (KDE auf Δρ + quadgk) | `ProteinCoLoc/src/bayes.jl` | Baseline-BF, gegen den der amortisierte Evidence-Network-BF validiert wird; definiert Δρ |
| `MultiChannelImage(...)`, `_apply_mask!`, `otsu_thresholds` | `ProteinCoLoc/src/LoadImages.jl` | Zielformat des Simulators (damit `correlation()`/`patch()` unverändert greifen) |
| `_bin_calibration()` + `CalibrationResult` | `BayesInteractomics/src/diagnostics/calibration.jl` | ECE/MCE-Gerüst + Traffic-Light-Schwellen für die SBC/Coverage-Auswertung |
| Parametric Simulation Engine (Sweep/Coverage, Platt, JLD2-Cache) | `BayesInteractomics/src/simulation/simulation.jl` | Vorlage für Szenario-Grid-Sweep, Replikate, Konfidenzbänder, Cache |

## Externe Tools
- **NeuralEstimators.jl** — primäre SBI-Engine: `PosteriorEstimator` (NPE via NormalizingFlow) + `RatioEstimator` (NRE/Bayes-Faktor) + `train`/`sampleposterior`/`assess`. **Deckt NPE+NRE in einem Julia-nativen Paket.** CPU reicht für den 2D-Spike (Training 50–200 k Paare in Minuten).
- **Flux.jl** — NN-Backend (Summary-Netz + Flow-Conditioner); reines Julia, Windows-tauglich. (Lux.jl als Alternative.)
- **NormalizingFlows.jl** (TuringLang) — optionale Flow/Summary-Dichte (OOD-Score), Bijectors-kompatibel.
- **InvertibleNetworks.jl** — echtes BayesFlow-Äquivalent in Julia, GPU-fähig; *nur* vorgemerkt für späteren 3D/Hierarchie-Vollausbau.
- **BayesFlow (Python) via PythonCall/CondaPkg** — **dokumentierter Fallback** falls Julia-Flow nicht konvergiert (CondaPkg-Muster aus BayesInteractomics vorhanden); bewusst aus dem Kernpfad gehalten.
- **Images.jl/ImageFiltering.jl** — PSF-Faltung, Registrierungs-Shift, Rauschen (bereits Dependency).
- **HypothesisTests.jl** — KS/χ²-Uniformitätstest der SBC-Ränge.

## Arbeitspakete
| ID | Titel | Kern-Tasks | Output | Aufwand |
|---|---|---|---|---|
| **AP1** | Stack + isoliertes Env + Smoke | `ProteinCoLoc/spike/` mit **eigenem Project.toml** (kein Touch am Haupt-Manifest); NeuralEstimators+Flux+Deps; **20-Zeilen-Smoke**: `PosteriorEstimator` mit NormalizingFlow auf 1-Param-Gauß `train`+`sampleposterior`; Stack-Decision-Notiz (Julia-nativ default, BayesFlow nur Fallback) | `spike/00_smoke.jl` grün; Stack-Notiz | 1 T |
| **AP2** | Physik-Vorwärtssimulator g(θ)→Bildpaar | θ = (ρ_true, Spillover s, Autofluoreszenz a, Labeleffizienz p, Shift dx/dy, Rausch-Level); Prozess: bivariate Dichten (Korr ρ_true) → Bernoulli-Ausdünnung → **PSF-Faltung** (Gauß, `imfilter`) → 2×2-Spillover-Mischung → Autofluoreszenz → Sub-Pixel-Shift → Poisson/Gauß-Rauschen; als `MultiChannelImage` verpacken; Prior π(θ) konsistent zu Turing-Priors | `spike/simulator.jl: simulate_pair(θ)` + Plausibilitätsplots (ρ_true↑→Korr↑) | 3–4 T |
| **AP3** | Trainingsdaten-Generator | θ~π → `simulate_pair` → `_prepare_data`/`correlation` → Summary-Vektor (Patch-Korr, optional + Median/IQR); 50–200 k Paare, JLD2-Cache (Muster aus simulation.jl); fixe Summary-Dim (8×8), Standardisierung; Held-out (5 k) | `spike/training_data.jld2` + Loader (Random123-Seed) | 2 T |
| **AP4** | NPE-Training + Baseline vs. ADVI | `PosteriorEstimator`/NormalizingFlow, Summary-Netz = MLP (DeepSet optional, Patches permutationsinvariant); `train`/`sampleposterior` (Fokus ρ_true/Δρ); **Baseline: für 20–30 Stacks echte `colocalization()`-ADVI** vs. NPE (RMSE, Intervallbreite, **Laufzeit ADVI-Minuten vs. NPE-ms**) | `spike/npe_model.jld2` + Vergleichstabelle (Genauigkeit, Speedup) | 4–5 T |
| **AP5** | **SBC / Coverage** (Alleinstellung) | M≈2000× θ*~π → simulieren → L Posterior-Draws → Rang von θ* → **Rang-Histogramm** je Parameter; Uniformität (KS/χ²) + Coverage-Kurve (nominal vs. empirisch); ECE/MCE via `_bin_calibration`; Reliability + SBC-Plot (Konfidenzbänder aus simulation.jl) | `spike/sbc.jl` + SBC-Histogramme/Coverage + Traffic-Light-Verdict | 3 T |
| **AP6** | Amortisierter BF (NRE) + OOD | `RatioEstimator`/Evidence-Network (l-POP-Loss, Jeffrey&Wandelt 2023): Coloc-Modell vs. Null → **amortisierter log-BF in einem Forward-Pass**; Validierung gegen `compute_BayesFactor()` (kein quadgk/KDE/Shuffle nötig); **OOD/Misspez** via Summary-Dichte (Mahalanobis/Flow-Likelihood) + Posterior-Predictive-Mismatch; Test mit absichtlich fehlspezifizierten Bildern (doppelter Spillover, fremde PSF) → Flag muss auslösen | `spike/evidence_net.jl`: log-BF + OOD-Flag; Vergleichsplot; korrekte OOD-Auslösung | 4 T |
| **AP7** | Spike-Report + Go/No-Go | `spike/demo.jl` verkettet alles; Erfolgskriterien tabellieren; **Entscheidungs-Memo** (lohnt Vollausbau: Hierarchie/3D/Multi-Kanal?); sauberer Schnitt (was in src/ übernommen würde, ohne zu mergen) | reproduzierbares Demo + 2–3-S.-Go/No-Go-Memo | 2 T |

## Meilensteine
- **M0 Stack steht** — NeuralEstimators-Toy-NPE trainiert+sampled (`00_smoke.jl` grün); Haupt-Repo unberührt · Ende Tag 1
- **M1 Simulator plausibel** — `simulate_pair`→MultiChannelImage; monotone ρ_true→Korr; Spillover/Shift sichtbar · Ende Wo 1
- **M2 NPE trifft ADVI** — RMSE vergleichbar bei **>100× Wall-Clock-Speedup** · Ende Wo 2
- **M3 Kalibrierung nachgewiesen** — SBC-Ränge uniform (KS p>0.05) *oder* Fehlkalibrierung sauber diagnostiziert+gemeldet; Coverage + Traffic-Light · Mitte Wo 3
- **M4 Amortisierter BF + OOD** — log-BF stimmt im gut-spez. Regime mit KDE-BF überein; OOD-Flag löst bei Fehlspezifikation aus · Ende Wo 3
- **M5 Go/No-Go** — Demo reproduzierbar; Memo mit Kennzahlen + Vollausbau-Entscheidung · Ende Wo 4

## Risiken & Mitigation
- **NeuralEstimators-Flow konvergiert nicht / API instabil** → AP1-Smoke als Gate vor Investment; Fallback NormalizingFlows.jl bzw. BayesFlow-via-PythonCall.
- **Summary-Statistik insuffizient** → mit der *exakt vorhandenen* Statistik starten (Vergleichbarkeit zu `colocalization()`); bei Bedarf Momente/Manders + DeepSet.
- **Simulator-Realismus-Gap** → *genau der Spike-Wert*: OOD/Misspez-Meldung (AP6) macht den Gap **messbar** statt ihn zu verstecken.
- **Fokus-Konflikt mit Einreichungen** → vollständig entkoppeltes Env, pausierbar, niedrigste Priorität.
- **SBC zeigt Fehlkalibrierung** → **kein Misserfolg, sondern publizierbarer Befund** + Trigger für Flow-Tuning; Traffic-Light meldet ehrlich (Seefelder-Signatur).
- **Windows-Flux/CUDA-Probleme** → Spike CPU-only (2D klein genug); GPU erst im Vollausbau.

## Definition of Done
- Eine trainierte NPE liefert für ≥20 Held-out-Stacks Posteriors (ρ_true/Δρ) in einem Forward-Pass, **>100× schneller** als pro-Datensatz-ADVI bei vergleichbarem RMSE.
- SBC-Rang-Histogramm + Coverage-Kurve quantitativ bewertet (KS/χ², ECE/MCE) + Traffic-Light — Kalibrierung nominal *oder* explizit als fehlend gemeldet.
- Amortisierter NRE-log-BF reproduziert im gut-spez. Regime `compute_BayesFactor()` — ohne quadgk/KDE/Shuffle.
- OOD-Flag löst auf fehlspezifizierten Bildern aus, bleibt auf In-Distribution ruhig.
- Alles reproduzierbar aus `spike/demo.jl` (fixer Seed); Haupt-Repo + beide Manuskript-Pipelines nachweislich unberührt.
- Go/No-Go-Memo mit Kennzahlen + konkretem Vollausbau-Plan.

## Tag-1-Schritte
1. `mkdir ProteinCoLoc/spike/`, dort eigenes Env (`Pkg.activate("spike")`) — Haupt-Project.toml/Manifest **nicht** anfassen.
2. `] add NeuralEstimators Flux Distributions JLD2 Images ImageFiltering KernelDensity QuadGK StatsBase Random123 HypothesisTests`.
3. `spike/00_smoke.jl`: trivialer 1-Param-Gauß, `PosteriorEstimator` mit NormalizingFlow, `train`+`sampleposterior` (<30 Zeilen) — bestätigt API+Flux-Backend.
4. `colocalization.jl`/`bayes.jl`/`LoadImages.jl` gezielt querlesen → genutzte Signaturen + Prior-Ränge in `spike/NOTES.md` kopieren.
5. `spike/simulator.jl`-Skeleton: θ-NamedTuple + `sample_prior()` (Ränge aus Turing-`@model`) + Stub `simulate_pair(θ)::MultiChannelImage`.
6. Stack-Decision-Notiz: NeuralEstimators.jl default, BayesFlow nur Fallback.

## Offene Entscheidungen (für dich)
- Summary-Umfang: nur Patch-Korr-Vektor vs. + Manders/Momente (Empfehlung: erst minimal, identisch zur ADVI-Pipeline).
- `num_patches` fixieren (z.B. 8×8) — feste Summary-Dim nötig; variabel erst im Vollausbau.
- Inferenz-Ziel: latentes ρ_true *oder* Δρ (Empfehlung: beide — Δρ für BF-Vergleich, ρ_true für SBC).
- NPE-only vs. + NRE/Evidence-Network (Empfehlung: beide; amortisierter BF ist Kern-Claim).
- Trainingsbudget N (50 k vs. 200 k) — an M2-Genauigkeit kalibrieren.
- Bei SBC-Fehlkalibrierung: Flow vergrößern/mehr Daten vs. als ehrlichen Negativbefund berichten — Schwelle vorab definieren.
