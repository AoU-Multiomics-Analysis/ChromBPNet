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
2. Supply a plain or gzip-compressed reference FASTA and a headerless,
   two-column chromosome-size TSV. The `genome_index` input is optional.
   If supplied, it must index the **uncompressed FASTA**, even when the FASTA
   input is compressed. Otherwise, each scoring task creates its own index.
   Do not supply a compressed-file index or a `.gzi` file.
3. Supply a TSV model manifest. See [examples/models.tsv](examples/models.tsv).
   The first line must contain `model_id`, `cell_type`, `model`, and `peaks`,
   in that order, separated by tabs. Model IDs must be
   unique and contain only letters, digits, underscores, or hyphens.
   Use a separate row and model ID for each fold or replicate.
4. Use standard, bias-corrected `chrombpnet_nobias.h5` models with one sequence
   input. SavedModel directories, model archive downloads, and ChromBPNet-lite
   models are not supported by this workflow.
5. Supply BED3 through narrowPeak files in the same reference build. For BED
   files with fewer than ten columns, the scorer uses the peak midpoint.
   Ten-column narrowPeak files must have valid summit offsets; `-1` is rejected.
6. Set `docker_image` to the digest of the tested Linux image.

Manifest fields are literal text. Do not add CSV-style quotes around them.
Use full cloud URIs for model and peak files in Terra. For example:

```tsv
model_id	cell_type	model	peaks
CD4_fold0	CD4 T cell	gs://bucket/cd4/model.h5	gs://bucket/cd4/peaks.bed
NK_fold0	NK	gs://bucket/nk/model.h5	gs://bucket/nk/peaks.bed
```

A metadata validation task checks the manifest without opening the model or
peak URIs. It returns a headerless TSV. The workflow reads that file and converts
each row to a `ModelSpec`. The `model` and `peaks` members have WDL `File` types.
Cromwell can thus localize them before scoring.
Each scoring task receives files as explicit `File` inputs and named CLI arguments.
The merge command writes its file list in the task execution directory, using
the localized score-file paths. It does not pass an engine-generated file URI
to the merge script.

Copy [examples/inputs.json](examples/inputs.json) and replace the example URIs
and image digest. Register the WDL in Terra, then use the input JSON for your
workspace submission. Ensure that the workspace can read all inputs and the image.
The example variant is a format example; check its REF against your selected genome.

For a pipeline test, use
[the 100-variant GM12878 ATAC TSV](examples/GM12878_ATAC.test_100_variants.tsv).
It contains synthetic SNPs in 100 distinct GM12878 peak intervals, with reference
alleles checked against UCSC hg38. See the
[test-set notes](examples/GM12878_ATAC.test_100_variants.README.md) for sources,
checks, and an upload command. Set `max_peaks` to 100 and keep `num_shuf` at 0
for a short execution test.

## Download GM12878 ATAC inputs

[tools/prepare_gm12878.py](tools/prepare_gm12878.py) prepares the five-fold ATAC
model set [ENCSR389HIH](https://www.encodeproject.org/annotations/ENCSR389HIH/),
trained on experiment ENCSR637XSC in GRCh38. It uses Python 3.9 or later with no
additional Python packages. Install and authenticate `gsutil` first. Your account
must have write access to the destination bucket.

```sh
python3 tools/prepare_gm12878.py \
  --output-dir ./GM12878_ATAC \
  --gcs-prefix gs://your-bucket/chrombpnet/GM12878_ATAC
```

The script downloads the ENCODE model archive `ENCFF142IOR` and region archive
`ENCFF971WEQ`. It verifies their published MD5 checksums. It selects only the five
bias-corrected HDF5 models and decompresses the full input peak set from the
region archive. The peak file keeps its ten narrowPeak columns and summit offsets.

The script also downloads the UCSC `hg38.fa.gz`, verifies its published checksum,
and creates matching chromosome sizes from the compressed FASTA. To use your
existing hg38 reference instead, add `--skip-reference`.

Files remain in the local directory. The script uploads each prepared data file
with `gsutil cp`, then creates and uploads `models.tsv`. The manifest has the
workflow's four columns and one row per fold. All five rows use `GM12878` as the
cell type and the full input peak set. Model IDs include the assay, experiment,
and fold. The script retains source URLs and checksums in local `sources.json`.
It reuses downloaded archives and FASTA files when their checksums match.
Uploads replace files at the selected destination paths.

Use these workflow inputs after the script succeeds:

| WDL input | GCS path |
| --- | --- |
| `model_manifest` | `gs://your-bucket/chrombpnet/GM12878_ATAC/models.tsv` |
| `genome` | `gs://your-bucket/chrombpnet/GM12878_ATAC/reference/hg38.fa.gz` |
| `chrom_sizes` | `gs://your-bucket/chrombpnet/GM12878_ATAC/reference/hg38.chrom.sizes` |

Leave `genome_index` unset. Supply your own variant list and tested scoring image
digest. If you use `--skip-reference`, also supply your own genome and chromosome
sizes. These uploads do not submit the WDL or start a Terra job.

## Download TenK10K ATAC inputs

[tools/prepare_tenk10k.py](tools/prepare_tenk10k.py) uses the same manifest format
as the GM12878 tool. It prepares **130 bias-corrected models: folds 0–4 for 26
cell types**. It also downloads the combined annotated peaks and the peak archive
for all 28 cell types. See [the TenK10K instructions](examples/TenK10K_ATAC.README.md)
for source details, file formats, and workflow inputs.

Keep the data on T7. Install and authenticate `gsutil`, then replace the bucket
path in this command:

```sh
python3 tools/prepare_tenk10k.py \
  --output-dir /Volumes/T7/AoU/ChromBPNet/TenK10K_ATAC \
  --gcs-prefix gs://your-bucket/chrombpnet/TenK10K_ATAC \
  --skip-reference
```

The tool verifies source checksums, converts peak CSV files to narrowPeak TSV,
and retains the two released BED6 files. Each model row uses the peaks for its
cell type. The ready `models.tsv` is uploaded after all data files. Set
`model_manifest` to `gs://your-bucket/chrombpnet/TenK10K_ATAC/models.tsv` in the
existing scoring WDL. Use your existing hg38 genome and chromosome sizes with
`--skip-reference`; omit that flag to prepare and upload the reference as well.

Add `--plan-only` to inspect the full manifest without downloads or uploads.
The resulting `models.planned.tsv` has future GCS paths and is not ready for Terra.
The bundled source inventory pins an immutable Hugging Face revision.

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
find a GPU. Compressed FASTA files reduce the data transferred during localization.
Each scoring task decompresses its localized FASTA before scoring. Select enough
task disk for both the compressed and uncompressed reference, its index, the
localized model, peaks, and outputs. The manifest and merge tasks do not request GPUs.

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

The `fold_summary` output is `variant_effects.fold_summary.tsv`. It has one row
per variant per model group, with arithmetic means of all numeric score columns
named `<metric>.mean`. This follows the upstream method for averaging fold scores.
For example, it provides `logfc.mean`, `jsd.mean`,
`active_allele_quantile.mean`, and `abs_logfc_x_jsd_x_active_allele_quantile.mean`.
Products and absolute values are averaged from the per-fold scores. They are
not calculated again from the averaged component scores. P-value columns are
excluded from this summary; the per-fold files retain them.

To group folds, use model IDs with a shared prefix and a final `_foldN` suffix,
where `N` is a nonnegative integer. The prepared GM12878 IDs already follow this
rule: `GM12878_ATAC_ENCSR637XSC_fold0` through `fold4` become the model group
`GM12878_ATAC_ENCSR637XSC`. Use a different prefix for each cell type, assay,
model set, or configuration that must remain separate. IDs without this suffix
produce separate one-model summaries. The script rejects duplicate fold indices,
groups with different cell types, and groups that mix suffixed and unsuffixed IDs.

The summary also provides:

- `n_folds` and `model_ids`: the actual number and IDs of the supplied models.
  Missing folds are not filled in. All supplied models must score the same variants.
- `logfc.sd`: sample standard deviation across folds, with denominator `n_folds - 1`.
  It is blank when only one model is supplied.
- `percent_change`: `100 * (2^logfc.mean - 1)`, for allele2 relative to allele1.
- `n_positive`, `n_negative`, and `n_zero`: fold counts by effect direction.
- `direction_agreement`: the larger of the positive and negative counts divided
  by `n_folds`. Zero effects remain in the denominator; all-zero effects give zero.
  This measures model agreement, not statistical significance.

To summarize an existing merged file without running GPU scoring again:

```sh
python scripts/summarize_folds.py \
  --scores variant_effects.all_models.tsv \
  --output variant_effects.fold_summary.tsv
```

The updated WDL requires an image that includes `scripts/summarize_folds.py`.
GitHub Actions builds and smoke-tests that image for this script change. After
the change is merged, use the newly published image digest in Terra.

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
It checks the TSV manifest, both merge outputs, and the default GPU requirement.
It compares scores from a plain FASTA with a supplied index against a compressed
FASTA with a task-created index. This test runs on CPU because the GitHub runner
has no GPU.

Run the other checks with Python 3.11 and `miniwdl==1.14.2`:

```sh
python -m unittest discover -s tests -v
miniwdl check workflows/score_variants.wdl
python tests/check_wdl.py workflows/score_variants.wdl
```

Tests check File types in the parsed manifest, command-time localization,
shell quoting, optional CLI arguments, and workflow-scope file writes.
The merge regression test starts with cloud File inputs, applies their local
path mapping, and executes the rendered command with the actual merge script.
It also models the cloud URI that Terra can return for `write_lines`. The static
check rejects WDL file-writing functions inside task command expressions;
commands must create required local files directly.
It also rejects functions outside the WDL 1.0 standard library. Miniwdl accepts
some newer functions in a 1.0 document, including `sep()`. Terra rejects that
function. The merge command uses the WDL 1.0 placeholder option
`~{sep="\n" score_files}` inside a quoted shell file block. GitHub Actions also
validates the workflow with Cromwell's `womtool` 85 and Java 17.
Miniwdl can warn that `predefinedMachineType` is unknown; Cromwell uses that field.

**The complete workflow has not been verified successfully on Terra.** A reported
Terra run reached the merge task but failed because its generated file-list path
remained a cloud URI. The first fix then failed Terra's parser because it used
`sep()`. The revised command creates the list locally with WDL 1.0 syntax. It
passed the command regression test, but has not been rerun on Terra. Syntax checks
and the CPU image smoke test do not validate L4 GPU execution. No Terra or other
analysis jobs were submitted for this fix.

## Dockstore

[.dockstore.yml](.dockstore.yml) uses Dockstore schema 1.2. It registers
`chrombpnet-variant-scoring` from the WDL, README, and example input files.
The configuration requests public publication with `publish: true` and limits
registration to the `main` branch.

Enable the [Dockstore GitHub App](https://docs.dockstore.org/en/stable/getting-started/github-apps/github-apps.html)
for this repository, then merge the configuration into `main`. Dockstore uses
the app connection to register and publish the workflow. Adding this file alone
does not confirm publication. Replace the example input URIs and image digest
before a Terra run.
