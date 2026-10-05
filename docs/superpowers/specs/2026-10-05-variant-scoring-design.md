# Terra ChromBPNet variant scoring

Use WDL 1.0 in Terra. Accept a headerless, five-column ChromBPNet variant TSV,
an indexed reference FASTA, chromosome sizes, and an array of manifest structs.
Each struct has a unique model ID, cell type, model File, and peaks File.
The manifest is structured workflow input data. Files stay typed until each task runs.

Validate the manifest before the scatter. Score all variants once per model with
the pinned official Kundaje variant-scorer. Use bias-corrected HDF5 models.
Peaks provide active allele quantiles; they do not filter variants.
Reject REF mismatches, invalid sequence windows, unreadable files, and lost variants.
Require a GPU in the scoring task. Request g2-standard-16 and nvidia-l4.
Keep per-model scores and logs. Merge into a long TSV and a wide TSV with model
IDs in score column names. Keep separate models, including folds, separate.

Build a pinned micromamba image with conda-forge and bioconda packages.
Build and smoke-test it in GitHub Actions, never locally.
Image rebuilds require changes to scripts contained in the image.
Check WDL syntax, workflow-scope file writes, task localization, and safe quoting.
A CPU synthetic-model smoke test checks the container; it does not validate L4 or Terra.
Open a pull request after local checks. Do not submit cloud jobs.
