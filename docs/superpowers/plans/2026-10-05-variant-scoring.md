# ChromBPNet Variant Scoring Implementation Plan

> Implement in this session with superpowers:executing-plans. Track steps below.

**Goal:** Score chromatin variant effects across models and merge cell-type results.
**Architecture:** Validate typed manifest, scatter GPU scores, merge CPU outputs.
**Tech Stack:** WDL 1.0, Python, pinned Kundaje variant-scorer, micromamba, GitHub Actions.
**Spec:** ../specs/2026-10-05-variant-scoring-design.md

## Global Constraints

Use Terra File localization. No workflow-scope file-writing functions.
Use named CLI arguments and task-local newline lists. Quote all task inputs.
Require g2-standard-16 and nvidia-l4. Never submit Terra jobs or build Docker locally.

## Review Focus

- Unresolved cloud URIs must fail before files are opened.
- Duplicate model IDs and empty manifests must fail before GPU tasks.
- Missing variants, duplicate scores, and allele changes must fail before merging.
- Apostrophes and shell expansion characters must remain literal in task paths.
- REF mismatches and indels near contig boundaries must fail before scoring.

### Task 1: Input and merge scripts

Files: scripts/io_utils.py, scripts/validate_manifest.py, scripts/merge_scores.py,
scripts/score_variants.py, tests/test_pipeline.py.
Interfaces: named CLI arguments; per-model labelled TSV; long and wide TSV outputs.
- [x] Write failing input/localization/merge tests and run unittest.
- [x] Implement validation, scoring wrapper, and merge. Run the complete unit suite.

### Task 2: WDL, image, and CI

Files: workflows/score_variants.wdl, containers/Dockerfile,
containers/environment.yml, tests/check_wdl.py, tests/smoke_image.py,
.github/workflows/test.yml, .github/workflows/image.yml, examples/inputs.json, README.md.
Interfaces: typed ModelSpec manifest; localized score File array for merge.
- [x] Add regression tests for workflow-scope file writing and command quoting.
- [x] Implement WDL, pinned image, CPU image smoke test, and Actions.
- [x] Run unit tests, miniwdl validation, static checks, and command rendering tests.
- [ ] Review all changes, commit, push, create and attach the pull request.

## Verification record

Fifteen unit tests passed. WDL syntax and workflow-scope write checks passed.
An independent review found two input-handling issues. Both were fixed with
regression tests: Java-compatible shell quoting, and preservation of reserved
or numeric-looking variant IDs. The real scorer CPU smoke test passed in an
existing local TensorFlow environment. The Linux image will be checked in Actions.
Terra and L4 execution have not been tested.

A task-local FASTA index copy also protects localized inputs when parallel
downloads make the index timestamp older than the FASTA. A regression test
checks that updates cannot change the input index.
