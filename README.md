# SP-008 / DOC-2-037 - Perturbation Transfer Learning

## Summary (from `SP-008-DOC-2-037/README.md`)


## Status
Closed useful negative after R1. Weak aggregate response-vector transfer exists, but discrete top-50 recovery and Replogle calibration failed.

## Useful result
Across five held-target/context tasks, leakage-safe aggregate prediction beat simple baselines by thin CI-positive margins. Yet top-50 response recovery failed every target-permuted-null comparison, and uncertainty calibration failed in both K562/RPE1 directions while passing the three Frangieh condition rotations. The boundary is sharp: weak aggregate transfer does not support discrete actionable recovery.

## What is new
The project separated aggregate-vector transfer from discrete-edge recovery under one frozen, multi-context framework. It links two apparent contradictions: some low-resolution signal transports, while the discrete features a scientist would act on do not.

## Why it matters
It prevents a statistically detectable aggregate gain from being inflated into a perturbation-selection product. The useful application is a transfer-readiness audit that blocks deployment when discrete recovery or calibration fails.

## Working application angle
A measurement-prioritization auditor can accept source/destination response matrices and emit baseline deltas, target-bootstrap uncertainty, top-k null comparisons, calibration, leakage assertions and an abstention decision. It must not recommend treatments or claim causality.

## Top-lab next question
Which prospective target descriptors, measured before destination response, can improve discrete recovery on an independent lab/context without leaking source response identity?

## Integrity boundary
Large OOF arrays and response matrices are hash-declared in `results/gate_evaluation_r1.json` but were not included in the small core archive. The core package is reproducible and auditable, but those external bytes require separate preservation before full graduation.

## Drive artifacts
- R0: https://drive.google.com/file/d/1LIYnzDC5gf1tV48xrNmfmvNAOedyXPJ1/view?usp=drivesdk&authuser=uditakankana%40gmail.com
- R1 core: https://drive.google.com/file/d/1u2lKHcBKiTvuNMA5vzrYAMiKGd2VUwj9/view?usp=drivesdk&authuser=uditakankana%40gmail.com

## Contents

- `SP-008-DOC-2-037/` - migrated unchanged from `science-program/projects/SP-008-DOC-2-037` (28 files)

## Provenance

Split out of the `science-program` repository (source commit `028a7141ed5f951a7b6e6517d4e72768d414a560`) on 2026-09-23. Every file is byte-identical to the source; `MIGRATION_MANIFEST.tsv` lists sha256, original path and new path for each of the 28 files.

Part of Udita Phookan's computational science program: every experiment locks its question, validation design, success gate and failure policy before outcome analysis, and negative results are preserved. Program-wide ledgers and standards live in the `science-program-ledger` repository.
