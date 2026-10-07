# ChromBPNet variant scoring in Terra

[workflows/score_variants.wdl](workflows/score_variants.wdl) scores predicted chromatin
effects for a variant list across cell-type models. The workflow uses WDL 1.0.
It runs one GPU scoring task per fold. It then runs one CPU summary task per
model group and a final CPU task to combine the outputs.

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

The scoring task defaults to 64 GB memory and 100 GB SSD. The boot disk
is 50 GB. `num_preempt` defaults to zero. The g2-standard-16 machine has one L4
GPU, so the GPU count is fixed at one. The script fails if TensorFlow cannot
find a GPU. Compressed FASTA files reduce the data transferred during localization.
Each scoring task decompresses its localized FASTA before scoring. Select enough
task disk for both the compressed and uncompressed reference, its index, the
localized model, peaks, and outputs. The manifest and merge tasks do not request GPUs.

Each model summary task has separate inputs: `summary_memory_gb` defaults to 64,
`summary_disk_gb` defaults to 500 (SSD), and `summary_max_retries` defaults to 2.
The final combine task has separate inputs: `merge_memory_gb` defaults to 64,
`merge_disk_gb` defaults to 500 (SSD), and `merge_max_retries` defaults to 2.
Two retries permit up to three attempts. Retries use the same memory request;
increase the memory input for the task that needs more memory.
These WDL resource inputs use the existing Docker image. Changing them does
not trigger an image rebuild. The image workflow runs only for `scripts/**` changes.

## Recover a failed merge

Use [workflows/merge_variants.wdl](workflows/merge_variants.wdl) to rerun only the
model summaries and final combination from saved score files. Register the
`chrombpnet-merge-variant-effects` workflow in Terra through Dockstore, or import
the WDL with its `summarize_variants.wdl` dependency. It has no scoring calls.
Supply the original TSV model manifest so each summary can use the correct
cell type and peaks. The workflow localizes peak files, but does not open model
files or need the genome or a separate variant input.

Copy [examples/merge_inputs.json](examples/merge_inputs.json), then replace its
`score_files` array with the full `gs://` paths for all completed
`call-ScoreVariants/shard-*/effects/variant_effects.tsv` files from the failed run.
Use the exact object paths shown in the call outputs; do not supply a wildcard,
a file-list TSV, or the already merged table. Include each model/fold file once,
**in the same order as the data rows in `model_manifest`**. A summary task checks
each score file's model ID and cell type against the manifest. An incorrect order
fails with a clear error instead of assigning a file to the wrong group.
The typed `Array[File]` input lets Cromwell localize every score file.
Set `docker_image` to the same image digest used by the current scoring workflow.
The summary and combine defaults are 64 GB memory, 500 GB SSD, and two retries.
The outputs are the long merge, wide merge, combined fold summary, per-model fold
summaries, and task logs. The peak annotation is included in each fold summary.

You can also resubmit the full workflow with call caching enabled, the same
scoring inputs and image digest, and only the merge inputs changed. Terra can
reuse successful scoring calls when their inputs and outputs match and remain
accessible. Changing scoring resource inputs can prevent these cache hits.
See [Terra call caching](https://support.terra.bio/hc/en-us/articles/360047664872-Call-caching-How-it-works-and-when-to-use-it).
The merge-only workflow avoids scoring calls even when a cache hit is unavailable.

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

The shared [workflows/summarize_variants.wdl](workflows/summarize_variants.wdl)
first creates a model-group plan from the manifest. That plan contains group names,
cell types, and row indices. Its JSON format represents structured output data;
it is not a task argument wrapper and contains no input file paths. The workflow
selects typed score and peak `File` inputs before task localization.

Each `SummarizeModel` task receives only the folds for one model group and its
peak file. It averages those folds and adds the peak annotation. `CombineModelEffects`
then combines these completed outputs. It does not calculate new means across
model groups. It streams the long and summary tables and stores the wide rows in
a temporary SQLite database on the task disk. Variant identities must match across
groups; variant order can differ. Existing long and wide outputs keep all folds.

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
All folds in a group must also use the same peak URI in the manifest.

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
- `in_peak`: `true` when the variant's reference interval overlaps at least one
  supplied peak for that model group; otherwise `false`. The input position is
  1-based. For a reference allele of length `L`, the tested BED interval is
  `[pos - 1, pos - 1 + L)`. An empty reference allele (`-`) uses one base at that
  position. Peak intervals are 0-based with an excluded end position. SNPs test
  one base; deletions test their full reference span. Peak starts count as overlap;
  peak ends do not. Chromosome names must match exactly. The annotation uses the
  full supplied peak file, regardless of the scoring task's `max_peaks` setting.

`per_model_fold_summaries` returns each model group's annotated summary before
combination. `summary_logs` and `grouping_log` report the summary and grouping steps.

To summarize an existing merged file without running GPU scoring again:

```sh
python scripts/summarize_folds.py \
  --scores variant_effects.all_models.tsv \
  --output variant_effects.fold_summary.tsv
```

The workflows reuse the existing image's `merge_scores.py`, `summarize_folds.py`,
and validation functions. New grouping, peak annotation, and combination code runs
inside WDL commands using Python's standard library. **This change does not require
an image rebuild.** Use the same digest as your current fold-summary pipeline.
The direct Python command above averages folds only; it does not add `in_peak`.

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
large variant lists if needed. Each model summary stores that group's fold scores
in memory; the final combine task stores rows on disk.
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
miniwdl check workflows/merge_variants.wdl
miniwdl check workflows/summarize_variants.wdl
python tests/check_wdl.py workflows/score_variants.wdl
python tests/check_wdl.py workflows/merge_variants.wdl
python tests/check_wdl.py workflows/summarize_variants.wdl
```

Tests check File types in the parsed manifest, command-time localization,
shell quoting, optional CLI arguments, and workflow-scope file writes.
The grouping, summary, and combination tests start with cloud File inputs, apply
their local path mapping, and execute the actual rendered task commands.
They check group selection, wrong model assignments, peak boundaries, indels,
missing or changed variants, and different variant orders across groups.
It also models the cloud URI that Terra can return for `write_lines`. The static
check rejects WDL file-writing functions inside task command expressions;
commands must create required local files directly.
It also rejects functions outside the WDL 1.0 standard library. Miniwdl accepts
some newer functions in a 1.0 document, including `sep()`. Terra rejects that
function. The merge command uses the WDL 1.0 placeholder option
`~{sep="\n" score_files}` inside a quoted shell file block. GitHub Actions also
validates the workflow with Cromwell's `womtool` 85 and Java 17.
Miniwdl can warn that `predefinedMachineType` is unknown; Cromwell uses that field.

**The per-model summary, peak annotation, and final combination have not been run
on Terra.** Local command tests exercise cloud-to-local File mapping and actual outputs.
Syntax checks and the CPU image smoke test do not validate L4 GPU execution.
No Terra or other analysis jobs were submitted for this change.

## Dockstore

[.dockstore.yml](.dockstore.yml) uses Dockstore schema 1.2. It registers
`chrombpnet-variant-scoring` and `chrombpnet-merge-variant-effects` from their WDLs,
README, and example input files.
The configuration requests public publication with `publish: true` and limits
registration to the `main` branch.

Enable the [Dockstore GitHub App](https://docs.dockstore.org/en/stable/getting-started/github-apps/github-apps.html)
for this repository, then merge the configuration into `main`. Dockstore uses
the app connection to register and publish the workflow. Adding this file alone
does not confirm publication. Replace the example input URIs and image digest
before a Terra run.
