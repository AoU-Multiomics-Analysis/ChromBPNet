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
- [ ] Write failing input/localization/merge tests and run unittest.
- [ ] Implement validation, scoring wrapper, and merge. Run the complete unit suite.

### Task 2: WDL, image, and CI

Files: workflows/score_variants.wdl, containers/Dockerfile,
containers/environment.yml, tests/check_wdl.py, tests/smoke_image.py,
.github/workflows/test.yml, .github/workflows/image.yml, examples/inputs.json, README.md.
Interfaces: typed ModelSpec manifest; localized score File array for merge.
- [ ] Add regression tests for workflow-scope file writing and command quoting.
- [ ] Implement WDL, pinned image, CPU image smoke test, and Actions.
- [ ] Run unit tests, miniwdl validation, static checks, and command rendering tests.
- [ ] Review all changes, commit, push, create and attach the pull request.
