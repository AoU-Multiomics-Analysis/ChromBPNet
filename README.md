# ChromBPNet variant scoring in Terra

[workflows/score_variants.wdl](workflows/score_variants.wdl) scores predicted chromatin
effects for a variant list across cell-type models. The workflow uses WDL 1.0.
It runs one GPU task per model and merges the score files in a CPU task.

## Inputs

1. Supply a **headerless TSV** with these five columns, in this order:
   `chr`, `pos`, `allele1`, `allele2`, `variant_id`.
   Use 1-based positions, uppercase alleles, and unique variant IDs.
   `allele1` must match the reference. Use `-` for an empty allele in an indel.
   Use chromosome names with the `chr` prefix. Names must match the FASTA,
   chromosome sizes, and peaks. The scorer does not change genome builds.
2. Supply an uncompressed reference FASTA, its matching `.fai` index, and a
   headerless, two-column chromosome-size TSV. All three are explicit `File` inputs.
3. Supply a JSON model manifest. See [examples/models.json](examples/models.json).
   Each row has `model_id`, `cell_type`, `model`, and `peaks`. Model IDs must be
   unique and contain only letters, digits, underscores, or hyphens.
   Use a separate row and model ID for each fold or replicate.
4. Use standard, bias-corrected `chrombpnet_nobias.h5` models with one sequence
   input. SavedModel directories, model archive downloads, and ChromBPNet-lite
   models are not supported by this workflow.
5. Supply BED3 through narrowPeak files in the same reference build. For BED
   files with fewer than ten columns, the scorer uses the peak midpoint.
   Ten-column narrowPeak files must have valid summit offsets; `-1` is rejected.
6. Set `docker_image` to the digest of the tested Linux image.

The JSON manifest is structured model data. It is not a script argument wrapper.
The workflow converts each row to a `ModelSpec`. The `model` and `peaks` members
have WDL `File` types. Cromwell can thus localize them before scoring. A metadata
validation task checks the manifest without opening the model or peak URIs.
Each scoring task receives files as explicit `File` inputs and named CLI arguments.
The merge file list is created during command rendering, after localization.

Copy [examples/inputs.json](examples/inputs.json) and replace the example URIs
and image digest. Register the WDL in Terra, then use the input JSON for your
workspace submission. Ensure that the workspace can read all inputs and the image.
The example variant is a format example; check its REF against your selected genome.

## GPU and storage

The scoring runtime requests:

```wdl
predefinedMachineType: "g2-standard-16"
gpuType: "nvidia-l4"
gpuCount: 1
cpu: 16
zones: ["us-central1-a", "us-central1-b", "us-central1-c", "us-central1-f"]
```

The default memory is 64 GB. The default task disk is 100 GB SSD. The boot disk
is 50 GB. `num_preempt` defaults to zero. The g2-standard-16 machine has one L4
GPU, so the GPU count is fixed at one. The script fails if TensorFlow cannot
find a GPU. Select enough task disk for the reference, localized model, peaks,
and outputs. The manifest and merge tasks do not request GPUs.

## Scores and outputs

The image pins [Kundaje variant-scorer](https://github.com/kundajelab/variant-scorer)
to commit `0e1e34199e63112aa618748bb79a206fc491300a`. Its source and license are
retained in `/opt/variant-scorer` in the image. Forward and reverse-complement
predictions are averaged. All input variants are scored, including variants
outside peaks. Peaks supply the distribution for active allele quantiles.
Peaks without a full model sequence window are excluded by the upstream scorer.
The workflow rejects a variant with an invalid window or a REF mismatch.

The output includes these upstream metrics:

- `logfc`: log2 predicted count ratio, allele2 / allele1. Positive values predict
  higher accessibility for allele2.
- `jsd`: Jensen–Shannon profile distance, with the upstream adjustment for indels.
- `active_allele_quantile`: predicted coverage percentile of the stronger allele
  relative to the supplied peaks.
- The upstream products of these metrics, including effect and prioritization scores.

The wrapper uses temporary internal IDs during scoring. It restores original IDs,
including `NA` and numeric-looking IDs, in main and shuffled score files.
Each per-model TSV adds `model_id` and `cell_type`. The workflow also returns peak
scores, run metadata, and logs. `merged_effects` has one row per variant per model.
`wide_effects` has one row per variant, with columns such as `CD4_fold0.logfc`
and `CD4_fold0.cell_type`. Folds and cell types remain separate. No averaging is
performed. Merging rejects missing variants, changed alleles, duplicate models,
and nonfinite core scores.

`num_shuf` defaults to zero, so no shuffle-derived P values are calculated.
Set it above zero to enable the upstream shuffled-sequence null distribution.
A small shuffle count is suitable for a smoke test, not stable P-value estimates.
Shuffle output files are returned when present. `max_peaks` is optional. If set,
it limits the number of valid peaks sampled for the quantile distribution.
`random_seed` defaults to 1234. `batch_size` defaults to 128.

The upstream scorer stores predictions in memory. Increase task memory or split
large variant lists if needed. The merge task also stores score rows in memory.
Predicted accessibility effects do not prove a causal expression effect.

## Container and checks

The Dockerfile uses `mambaorg/micromamba:2.3.2`. Conda packages use pinned versions
from conda-forge and bioconda. TensorFlow 2.15.1 and its CUDA dependencies use
TensorFlow's `and-cuda` wheel extra. TensorFlow Probability is pinned to 0.23.0.
See [containers/environment.yml](containers/environment.yml).

GitHub Actions builds and smoke-tests the image on Linux. The image workflow runs
only when a script under `scripts/` changes. WDL, documentation, tests, or container
recipe changes alone do not rebuild the image. When changing the image recipe or
its packages, include the related change to the script that uses them.
Pull requests test the image. After a merge to `main`, the image is published to
`ghcr.io/aou-multiomics-analysis/chrombpnet/variant-scoring:<commit>` only if the
smoke test passes. The Actions summary gives its digest. Use that digest in Terra.
Make the package readable by the Terra runtime, or copy the tested image to an
accessible registry and use that registry's digest.

The image smoke test creates a synthetic HDF5 model. It runs the actual scorer
on SNPs, an indel, peaks, shuffled sequences, and an identical-allele control.
It checks both merge outputs and the default GPU requirement. This test runs
on CPU because the GitHub runner has no GPU.

Run the other checks with Python 3.11 and `miniwdl==1.14.2`:

```sh
python -m unittest discover -s tests -v
miniwdl check workflows/score_variants.wdl
python tests/check_wdl.py workflows/score_variants.wdl
```

Tests check File types in the parsed manifest, command-time localization,
shell quoting, optional CLI arguments, and workflow-scope file writes.
Miniwdl can warn that `predefinedMachineType` is unknown; Cromwell uses that field.

**The complete workflow has not been tested on Terra.** Syntax checks and the
CPU image smoke test do not validate Terra localization or L4 GPU execution.
No Terra or other cloud jobs have been submitted.
