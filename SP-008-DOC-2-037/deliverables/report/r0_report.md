# SP-008 (DOC-2-037) Perturbation Transfer Learning - R0 Feasibility

**Verdict: GO to R1 locked benchmark (4/4 gates PASS).** Protocol sha256 bd7e29df45f71f910ca169add726720d8fc421cd8659415c9cbafc57b90b93a8 (locked before inspection).

Lineage: SP-007 showed discrete supported-edge transport fails across cell lines. This project targets a different estimand: predicting aggregate perturbation-response vectors across contexts (held-out target primary, held-out context secondary). SP-007 tau/support thresholds are NOT reused.

- G1: Replogle K562 vs RPE1 (scPerturb Zenodo 10044268, CC-BY-4.0, md5-verified in SP-007): 2,055 shared targets, 10,691+11,485 NTC. PASS.
- G2: Frangieh 2021 melanoma (md5 dc438f53c476a3dd562898654d48ed15 verified): 218,331 cells x 23,712 genes; condition axis Control/IFN-gamma/Co-culture with 15,361/23,910/18,334 control cells; 229 targets with >=30 cells in all three conditions. Plus DatlingerBock2021 Jurkat (md5 5cfbcb4770e9c859a07b2cc44fd12066). PASS.
- G3: estimand, splits (held-out target primary, held-out context secondary, never random cells), metrics (Pearson, R^2 vs context-mean, top-50 recovery, calibration), baselines (context-mean, cross-context mean, ridge), nulls (target-permuted, expression-matched, shuffled) frozen pre-outcome. PASS.
- G4: pinned stack (anndata 0.11.4, h5py 3.16.0, sklearn 1.7.2, numpy 2.2.6, pandas 2.3.3, scipy 1.15.3, wheel sha256 in ledger); torch 2.14.0 downloads (554MB) but exceeds the 1GB RAM ceiling for import+training - GEARS-class GNN methods documented as compute ceiling; R1 method family: algebraic transfer (ridge/kernel-ridge, low-rank factorization, k-NN). PASS.

Abstentions carried: no causal graph claims, no discrete-edge transport claims, cell-state claims abstained, no SP-007 threshold reuse.
