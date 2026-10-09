# TenK10K models and peaks for the scoring workflow

Run [the preparation tool](../tools/prepare_tenk10k.py) with Python 3.9 or later.
Keep it with [the GM12878 helper](../tools/prepare_gm12878.py) and
[the source inventory](../tools/tenk10k_sources.json) in the `tools` directory.
No extra Python packages are required. Install and authenticate `gsutil` before
uploading. Your account must have write access to the destination bucket.

From this repository directory, use:

```sh
python3 tools/prepare_tenk10k.py \
  --output-dir /Volumes/T7/AoU/ChromBPNet/TenK10K_ATAC \
  --gcs-prefix gs://your-bucket/chrombpnet/TenK10K_ATAC \
  --skip-reference
```

Replace `your-bucket` with your bucket name. T7 must be mounted before the tool
writes this directory. All downloads, temporary download files, converted peaks,
and local output manifests stay in that output directory. Source code stays local.

This command uses your existing hg38 reference in Terra. Omit `--skip-reference`
to download and upload the same UCSC reference used by the GM12878 tool. The tool
derives matching chromosome sizes from the FASTA. Leave `genome_index` unset if
you want the scoring task to create its index.

To inspect the planned files first, add `--plan-only`. This creates
`models.planned.tsv` and `download_plan.json` without downloading or uploading
data. It does not create a ready `models.tsv`.

## Released inputs

The source is the [official TenK10K Hugging Face dataset](https://huggingface.co/datasets/anglixue/TenK10K_multiome).
The inventory pins commit `ba186a3790eb1a8298a08e31c17f16df5b6625b0`.
Each download uses that commit in its URL. Its byte count and SHA-256 checksum
must match the inventory. A repeated run reuses verified files.

The tool downloads:

- All 130 `chrombpnet_nobias.h5` files: folds 0–4 for each of 26 cell types.
- The original `TenK10K_ATAC_ct_specific_peaks.zip`, with peaks for 28 cell types.
- `TenK10K_ATAC_MACS3_Combined_Peaks_Annotated.bed`, retained without changes.

These source files total about 3.50 GB. Converted peaks need additional space.
The optional reference needs additional space as well. The tool selects the
bias-corrected models required by the scoring WDL. The full models with enzyme
bias and the separate bias models are not scoring inputs for this workflow.

The peak archive has mixed formats. Most files are CSV tables with an index
column and ten narrowPeak fields. The tool removes the CSV index and header,
writes tab-delimited narrowPeak, and retains the supplied summit offset. The
`CD14_Mono.csv` and `CD4_TCM.csv` members are headerless BED6 files. The tool
preserves those six columns and uses `.bed` output names. The existing scorer
centers BED6 peaks at their midpoint. It uses the summit for narrowPeak files.
The original archive remains available for reference.

All 28 peak files are converted and uploaded. `ILC` and `CD8_Proliferating` have
released peaks but no released models at this revision. They have no manifest
rows. The tool reports this difference in its log and `sources.json`.
The combined annotated peak file is an annotation resource. The manifest pairs
each model with the matching cell type peak file.

## Manifest and Terra inputs

The tool creates and uploads `models.tsv` after the data uploads succeed.
It has 130 data rows. The four columns are:

```text
model_id    cell_type    model    peaks
```

The separator is a tab. The five folds have separate model IDs. Each row uses
the same cell type label as the released model folder. The `model` and `peaks`
values are complete `gs://` URIs. See [the full example manifest](models.tenk10k.example.tsv).
Its bucket paths are examples; use the manifest produced by your successful run.

Use the existing [WDL 1.0 workflow](../workflows/score_variants.wdl).
Its `ModelSpec` keeps model and peak inputs as `File` values before task
localization. This preparation tool does not require a WDL or image change.

| Workflow input | Value |
| --- | --- |
| `model_manifest` | `gs://your-bucket/chrombpnet/TenK10K_ATAC/models.tsv` |
| `variants` | Your headerless, five-column variant TSV in hg38 |
| `genome` | Your existing hg38 FASTA; or the uploaded `reference/hg38.fa.gz` |
| `chrom_sizes` | Sizes matched to that FASTA; or the uploaded `reference/hg38.chrom.sizes` |
| `docker_image` | The tested scoring image digest used for your GM12878 run |

All other inputs can use the existing settings. The full manifest starts one
GPU scoring task for each of the 130 models. The workflow keeps folds separate
in the long and wide merged effects. A separate summary task averages the five
folds within each cell type's model group and adds `in_peak` from its supplied
peak file. The final task combines those summaries without averaging cell types.

`sources.json` records the source revision, source checksums, converted peak
counts and checksums, and cell types without models. It is source data provenance,
not a JSON wrapper for arguments passed to a WDL task. WDL continues to use typed
inputs and named CLI arguments.

Uploads replace objects at the destination paths. If a download or upload fails,
run the same command again. Use a new GCS prefix when preparing a different source
revision. Do not use `models.planned.tsv` to start scoring before the upload ends.

## Validation

The tests check file verification, both peak formats, missing and duplicate folds,
missing cell types, unsafe archives, manifest parsing, GCS file pairing, upload
failure, and the T7 mount check. The existing GitHub Actions test workflow runs
these tests with the other script and WDL checks. This tool is not part of the
scoring image, so its changes do not require an image rebuild.

All 28 actual peak files were converted and validated: 4,546,975 peak records
in total. The original peak archive, combined peak file, and one CD14 monocyte
fold were downloaded and checked against their source SHA-256 checksums.
The full 130-model download and GCS upload have not been run in this task.
This TenK10K manifest has not been tested in a complete Terra workflow run.
No Terra job was submitted.
