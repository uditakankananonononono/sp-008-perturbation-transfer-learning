# SP-008 / DOC-2-037 Perturbation Transfer Learning - R1 Report

## Verdict
**R1 STOP.** Required gates failed: G4 (top-50 discrete recovery, all 5 tasks) and G5 (uncertainty calibration, both Replogle directions). G3 aggregate baseline-superiority passed all 5 tasks with thin margins; under the locked protocol this is reported as an exploratory signal only and carries no product claim. No destination-response leakage occurred at any point (G2 pass); the discrete-recovery product arm is closed on G4 under routing point 8, and the aggregate signal is preserved as exploratory.

## Locked references
- Protocol R1 sha256: `984cff11fae70f3cfce7a6189b2e1c7bf424490caed557e1cdf33dce090be285` (locked 2026-09-22 03:11:47 IST, before any outcome inspection)
- Task-freeze manifest sha256: `7d23dd661e55057b4b8258df71b385f92fd4b444b1e3c65a333c3c8af995fcf5`
- Response-matrix construction frozen in protocol (pooled-count rate vs matched NTC control rate per context); no SP-007 thresholds reused.

## Data and provenance
- Replogle K562 + RPE1 (figshare.plus 20029387 v1, CC BY 4.0; md5 d8cba17576d1a8afc0f7d71b79cad0f7 / cc7f1ec50aeb3a3e1b4a6cfa713d80fa), reused from SP-007 with provenance ledger.
- Frangieh 2021 melanoma (scPerturb Zenodo 10044268, CC-BY-4.0; md5 dc438f53c476a3dd562898654d48ed15; 218,331 cells x 23,712 genes; condition axis Control / IFNg / Co-culture). Used as an independent condition axis only - not lab replication of the Replogle cell-line task.
- DatlingerBock2021 Jurkat (md5 5cfbcb4770e9c859a07b2cc44fd12066) - evaluated only at the pre-outcome overlap gate (G7).

## Tasks (frozen)
- Replogle: 115 shared strict targets x 7,226 shared genes; K562->RPE1 and RPE1->K562.
- Frangieh: 229 targets x 11,790 genes; three leave-one-condition-out rotations.
- 5-fold target-held-out CV, seed 20260922. Descriptor = source-context response vector only; unseen-target fallback = context mean.
- Models: context_mean, xcontext_mean (baselines), ridge, kernel_ridge, lowrank, knn. Dual-form solvers throughout (1GB RAM ceiling; torch-class frameworks documented as infeasible in environment ledger).

## Gate results

### G3 baseline superiority - PASS (thin margins)
Candidate vs best baseline (context_mean), target-bootstrap 95% CI:
| task | candidate | Pearson diff | CI | R2 diff | CI |
|---|---|---|---|---|---|
| K562->RPE1 | kernel_ridge | +0.0050 | [+0.0036, +0.0065] | +0.0211 | [+0.0149, +0.0273] |
| RPE1->K562 | lowrank | +0.0446 | [+0.0309, +0.0593] | +0.0218 | [+0.0054, +0.0379] |
| Frangieh IFNg+Co-cult->Control | kernel_ridge | +0.0033 | [+0.0022, +0.0044] | +0.0018 | [+0.0009, +0.0028] |
| Frangieh Control+Co-cult->IFNg | kernel_ridge | +0.0087 | CI excl 0 | CI excl 0 | |
| Frangieh Control+IFNg->Co-cult | kernel_ridge | +0.0113 | CI excl 0 | CI excl 0 | |

Per-target Pearson deciles (candidate): K562->RPE1 0.019/0.195/0.420/0.511/0.671; RPE1->K562 0.047/0.097/0.168/0.253/0.431; Frangieh->Control 0.161/0.246/0.297/0.349/0.396. Margins are small; median absolute Pearson remains moderate (0.17-0.42).

### G4 top-50 recovery - FAIL (all 5 tasks)
Candidate top-50 overlap with observed top-50 is indistinguishable from the target-permuted null:
| task | recovery median | permuted null median | null p95 |
|---|---|---|---|
| K562->RPE1 | 0.08 | 0.08 | 0.26 |
| RPE1->K562 | 0.04 | 0.06 | 0.16 |
| Frangieh ->IFNg | 0.20 | 0.18 | 0.28 |
| Frangieh ->Co-culture | 0.14 | 0.14 | 0.24 |
| Frangieh ->Control | 0.20 | 0.20 | 0.34 |
The expression-matched null recovers nothing (0.0), so beating it is trivial and not evidence of signal. The gate required beating both nulls; the permuted null is not beaten in any task.

### G5 uncertainty calibration - FAIL Replogle / PASS Frangieh
q90-of-residuals intervals, per-fold coverage tolerance [0.80, 0.97]: both Replogle directions have at least one fold outside tolerance for every model (fold means 0.878-0.913). All three Frangieh rotations pass for all models (0.886-0.899).

### G2 leakage - PASS
Folds disjoint and exactly covering (asserted); descriptors source-context-only (code-asserted at fit time); OOF/y shape assertions passed for all 30 model-task files.

### G6 reproducibility - PASS
Independent rerun reproduced run-2 gate metrics exactly (kernel_ridge K562->RPE1 median 0.4199==0.4199; lowrank RPE1->K562 0.1682==0.1682; Frangieh medians exact).

### G7 Jurkat bonus - EXCLUDED at pre-outcome gate
Jurkat panel = 20 TCR-pathway genes; zero overlap with the 115 Replogle shared strict targets. Recorded, not evaluated (results/r1_jurkat_gate.json).

## Abstention classes
1. Targets with no valid pre-intervention descriptor fell back to context mean by design (frozen rule) and are included in all metrics - not scored as model wins.
2. Jurkat arm abstained entirely (G7).
3. GEARS-class graph deep-learning methods abstained: torch 2.14.0 (554MB wheel) exceeds the 1GB RAM ceiling; documented in environment ledger, not silently skipped.

## Relation to SP-007 negative evidence
SP-007 (DOC-2-033) R2 showed matched-null edge support does not transport across cell lines (supported-set precision 0.0295 vs null p95 0.0325). SP-008 R1 refines that boundary: aggregate response-vector prediction beats baselines weakly but consistently across 5 tasks, yet discrete top-50 recovery does not beat a target-permuted null anywhere. The negative boundary now sits between aggregate and discrete estimands. Both results are preserved as delivered.

## Product framing (honest)
No product claim graduates. The exploratory aggregate signal is suitable only for measurement-prioritisation hypotheses with uncertainty attached, and even there the thin margins (+0.003 to +0.045 Pearson) and Replogle calibration failure argue against any use-level claim. Frangieh results are condition-axis evidence (melanoma), not replication of the cell-line task.
