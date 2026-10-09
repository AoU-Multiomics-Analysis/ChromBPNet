# ChromBPNet TF-MoDISco in Terra

[../workflows/discover_motifs.wdl](../workflows/discover_motifs.wdl) discovers
motifs for one cell type using five trained ChromBPNet folds. It uses WDL 1.0.
The workflow does not train models or annotate variant alleles. Its motif file
can be used in a later Fi-NeMo analysis of peak or REF/ALT contributions.

## What happens to the five folds

1. A CPU task checks the five model files and the reference. It prepares the
   common peak file once, excludes windows that cross chromosome boundaries,
   and shuffles peak order with a fixed seed.
2. Five GPU tasks compute hypothetical per-base DeepSHAP contributions for
   the same sequences. Each task uses one trained fold. The dinucleotide
   background count and shuffle seed are the same across folds.
3. A CPU task checks exact sequence identity, coordinates, row order, array
   dimensions, fold IDs, cell type, and output head. It computes the arithmetic
   mean of the **signed hypothetical contributions** in row chunks. It copies
   the common sequences and computes projected contributions from the mean.
4. TF-MoDISco runs once on the averaged contributions.
5. A report task generates sequence motif matrices, matches them to the
   supplied reference database with MEME Tomtom, and generates the HTML report.

The [ChromBPNet FAQ](https://github.com/kundajelab/chrombpnet/wiki/FAQ) recommends
averaging contribution HDF5 files before TF-MoDISco. The
[HDMA workflow](https://greenleaflab.github.io/HDMA/code/03-chrombpnet/) also
averages fold contributions separately for counts and profile before discovery.

The TF-MoDISco step takes contribution arrays, not five model files directly.
Do not average absolute scores, scalar variant log2FC scores, or motif matrices
from separate discovery runs as substitutes for this step.

## Required inputs

Copy [examples/inputs.cd4.json](examples/inputs.cd4.json), replace all example
paths, and supply these inputs to Terra:

| Input | Use |
| --- | --- |
| `fold_models` | Array of five distinct bias-corrected Keras HDF5 model files, ordered folds 0–4. |
| `peaks` | Full peak set for the selected cell type, in headerless BED3–BED10/narrowPeak TSV. |
| `genome` | Matching reference FASTA, plain or gzip/BGZF compressed. |
| `genome_index` | Optional index of the **uncompressed** FASTA. An index is built if absent. |
| `chrom_sizes` | Headerless chromosome name and length TSV matching the FASTA. |
| `motif_database` | Known human motif database in MEME format, such as JASPAR or a named MotifCompendium release. |
| `cell_type` | Label for the model and peak set, such as `CD4_Naive`. |
| `docker_image` | Dedicated tested image digest; see the GitHub image job. |

The example JSON is a Terra submission input file. It is not passed to a task
script. Every script receives safely quoted, named CLI arguments.

Use files from the same genome build. The workflow does not convert builds.
It checks chromosome lengths and names. It uses narrowPeak summit offsets
when ten columns are supplied. A summit of `-1` is rejected. For BED files
with fewer columns, it uses the interval midpoint. Duplicate chromosome/center
windows are rejected. Other narrowPeak scores are not used for filtering.

Use the peaks supplied with the model training set when available. For
TenK10K, supply the peaks for the exact cell subtype, not an aggregate CD4
peak set with a different model. The existing model manifest has one row per
fold; copy the five selected model URIs into `fold_models` and their shared
peak URI into `peaks`. Cromwell localizes these explicit WDL File inputs.
Paths stored inside a TSV or JSON file are not task inputs in this workflow.

Models must have one sequence input and the standard profile/counts outputs.
SavedModel directories, model download archives, bias-only models, and
ChromBPNet-lite models are not supported. `input_length` defaults to 2114;
change it only when all five models have a different input length. Model files
with duplicate contents are rejected.

## Settings

| Setting | Default | Meaning |
| --- | --- | --- |
| `head` | `counts` | Explain total accessibility. Use `profile` in a separate submission to explain profile shape. |
| `num_backgrounds` | 20 | Dinucleotide-shuffled reference sequences per interpreted sequence. |
| `batch_size` | 64 | Peak rows interpreted at once; also used for averaging chunks. |
| `random_seed` | 1234 | Peak-order shuffle and dinucleotide-background seed. |
| `max_peaks` | unset | Use all valid peaks. Set a limit only for a trial run. |
| `discovery_window` | 400 | Central bases searched by TF-MoDISco; the full model window is scored and retained. |
| `max_seqlets` | 1000000 | Maximum seqlets per positive/negative group. |
| `n_leiden` | 2 | Leiden clustering runs used by TF-MoDISco. |
| `match_qvalue` | 0.05 | Tomtom q-value threshold for the exported candidate TF table. |
| `n_matches` | 5 | Top candidate matches displayed in the HTML report. |

Counts and profile contributions remain separate. Submit the workflow twice
if both are required. Do not combine cell types in the fold array. A broader
motif library across cell types can be assembled after individual discovery
runs; that is a separate analysis.

The seqlet cap follows input order, so peaks are shuffled once before scoring.
The mean contribution file retains that order. Keep `interpreted_regions`
with the HDF5 files to map sequence rows and seqlets to genomic coordinates.

## Runtime and storage

The five contribution tasks request `g2-standard-16`, one NVIDIA L4 GPU,
16 CPUs, and 64 GB memory. They fail if TensorFlow cannot detect a GPU.
The CPU-only flag is available in the Python script for the CI smoke test;
the WDL never enables it. No other task requests a GPU.

TF-MoDISco defaults to 16 CPUs and 128 GB memory. These limits are configurable.
Motif discovery can use substantially more memory as the seqlet cap increases.
The averaging task processes row chunks, but TF-MoDISco loads its input arrays
and constructs similarity matrices in memory. Review the logs and adjust
resources for a full peak set. The default discovery and contribution disks
are 100 GB; averaging disk is 150 GB because it localizes five contribution
files and writes a sixth. The boot disk is 50 GB. Inputs, reference staging,
intermediate files, and final outputs must fit on each task disk.

## Outputs and interpretation

| Output | Use |
| --- | --- |
| `modisco_motifs` | TF-MoDISco v2 HDF5 motif patterns/CWMs for Fi-NeMo. |
| `averaged_contributions` | Mean `shap/seq`, common `raw/seq`, and `projected_shap/seq`, with arrays `(regions, 4, input_length)` in ACGT order. |
| `interpreted_regions` | Ordered narrowPeak regions corresponding to the HDF5 rows. |
| `report_bundle` | ZIP containing the report and all report assets. Extract and open `report/report.html`. |
| `report_html` | HTML alone; use the bundle when local image assets are required. |
| `discovered_motifs_meme` | Sequence probability matrices trimmed at 30% of maximum absolute CWM importance. These are for motif matching, not the Fi-NeMo CWM input. |
| `candidate_tf_matches` | All Tomtom matches with q-values at or below `match_qvalue`. |
| `tomtom_results` | Raw MEME Tomtom table for the trimmed motif queries. |
| `motif_inventory` | Pattern names, HDF5 group paths, seqlet counts, and trim coordinates. |
| Metadata and logs | Model/peak checksums, versions, settings, row counts, and task logs. |
| Per-fold contributions | Original fold HDF5 files for inspection or sensitivity analyses. |

Pattern IDs use names such as `pos_patterns.pattern_0`; the inventory maps
these IDs to HDF5 paths such as `pos_patterns/pattern_0`. The HDF5 patterns
keep their original IDs. They are not renamed after the top TF match.

The HTML report shows its top similarity matches, which can include weak
matches. Use the exported table and its q-values for thresholded annotations.
A motif match gives candidate TFs or TF families. Similar motifs can belong
to multiple TFs; cell-type expression and binding data can help resolve them.
Positive/negative patterns refer to predicted accessibility contributions,
not the direction of methylation change. A match does not prove TF occupancy
or a causal effect on methylation.

## Dedicated GitHub image and checks

The separate [Dockerfile](containers/Dockerfile) uses
`mambaorg/micromamba:2.3.2`, conda-forge/bioconda, and pinned packages.
It pins TF-MoDISco 2.5.2, TensorFlow 2.15.1, MEME 5.5.7, and the ChromBPNet
attribution source at commit `09938fdb4397ec0006510e5251e48920a505d4de`.
NumPy 1.23.5 is retained for the older Kundaje DeepSHAP fork. The image imports
the pinned ChromBPNet attribution helpers without installing its legacy
training dependency stack. The upstream source and license remain in the image.
Pinned igraph and Leiden wheels avoid the conflicting ICU requirements of
the MEME and conda graph packages. Discovery metadata records their versions.

[The image workflow](../.github/workflows/motif-image.yml) runs only when
`motif_pipeline/scripts/**` changes. It builds on Linux in GitHub Actions and
runs actual counts and profile DeepSHAP, five-fold averaging, TF-MoDISco,
Tomtom, and reporting on synthetic data. The read-only BuildKit test mount
does not add test files to the image. Tests run before image export and
publication. They also check that CPU use is rejected by default.
No local Docker build is required. Recipe-only, test-only, and WDL-only
changes do not trigger a rebuild; include the related script change when an
image recipe or package update is needed.

PRs build and test the image. A qualifying push to `main` publishes the tested
image to `ghcr.io/aou-multiomics-analysis/chrombpnet/motif-discovery:<commit>`.
The GitHub job summary supplies the digest. Use that digest in Terra. Ensure
that Terra can read the package, or copy the tested image to an accessible
registry. The existing `variant-scoring` image is not changed by this workflow.

The script/WDL workflow runs unit tests, miniwdl validation, static checks
against workflow-scope file writes, and Cromwell womtool 85 validation.
Regression tests map cloud File inputs to local paths before rendering every
task command. They run the real preparation and averaging scripts, and check
readable CLI paths for the remaining tasks. Preparation covers the optional
FASTA index. File lists are written within task commands from individually
quoted paths. A regression test checks that paths with line breaks cannot
run shell commands or enter these lists.

**The complete workflow has not been tested on Terra.** Local syntax checks
and the GitHub CPU smoke test do not validate managed Cromwell localization,
NVIDIA L4 execution, compatibility with each released TenK10K model, or full
peak-set resource needs. No Terra analysis job has been submitted.
