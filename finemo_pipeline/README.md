# Fi-NeMo motif calls and variant annotation in Terra

This pipeline uses a cell-type TF-MoDISco motif library to call motif instances
from ChromBPNet contributions. It has two WDL 1.0 entry points:

- [call_motifs.wdl](../workflows/call_motifs.wdl) uses existing five-fold averaged
  peak contributions from [motif discovery](../motif_pipeline/README.md).
- [annotate_variant_motifs.wdl](../workflows/annotate_variant_motifs.wdl) generates
  five-fold REF/ALT contributions, calls motifs, and compares allele calls.

Run one cell type and one output head per submission. Counts is the default.
Use profile in a separate submission. The variant workflow does not require
ATAC peak overlap, so variants outside peaks remain eligible.

## Inputs

Copy [inputs.variants.cd4.json](examples/inputs.variants.cd4.json) or
[inputs.peaks.cd4.json](examples/inputs.peaks.cd4.json) into Terra. Replace all
example paths and image placeholders. JSON is used only for Terra submission;
tasks receive typed inputs and safely quoted named CLI arguments.

| Input | Use |
| --- | --- |
| `motifs` | The `modisco_motifs` HDF5 output from the matching cell-type discovery run. Use the CWM file, not the exported MEME sequence matrices. |
| `motif_provenance` | The `preparation_metadata` JSON output from that same discovery run. This is provenance data, not an argument wrapper. |
| `tf_matches` | Optional `candidate_tf_matches` TSV output from that run. All candidate matches are retained. |
| `cell_type`, `head` | Must match discovery and contribution provenance. |
| `contributions`, `regions` | Peak workflow only: exact `averaged_contributions` HDF5 and `interpreted_regions` BED outputs from discovery. Do not sort or filter one without the other. |
| `fold_models` | Variant workflow only: five distinct bias-corrected HDF5 models, ordered folds 0–4. Their hashes must match the discovery provenance. |
| `variants` | Variant workflow only: headerless TSV with chromosome, 1-based position, REF, ALT, and unique variant ID. |
| `genome`, `genome_index`, `chrom_sizes` | Variant workflow only: matching FASTA, optional index of its uncompressed form, and chromosome sizes. Plain and gzip/BGZF FASTA are accepted. |
| `contribution_image` | Variant workflow only: the tested updated motif-discovery image containing the variant attribution scripts. |
| `finemo_image` / `docker_image` | Tested dedicated Fi-NeMo image digest for variant / peak workflow. |

Use nonempty uppercase ACGT alleles in the REF/ALT columns, including the VCF
anchor for insertions and deletions. `-`, symbolic alleles, multiallelic ALT
strings, and identical REF/ALT are rejected. Convert these inputs before
submission. The workflow checks REF against the genome and requires a complete
model window for both alleles. It does not convert genome builds.

The averaged contribution HDF5 must contain the exact region-file hash.
Older mean files from the first discovery image do not have this attribute.
Regenerate those means from the saved per-fold contribution files and exact
region BED with the updated averaging script/image. No GPU scoring is needed
for this regeneration.

The variant window starts half a model-input length before the first REF base.
An ALT insertion replaces REF and truncates the right edge to the model length.
A deletion takes additional reference bases at the right edge. Coordinate maps
match unchanged allele prefixes and suffixes; replacement bases pair from the
left. Extra inserted bases are marked `-1` in the map. Large alleles at least
half the model-input length are rejected.

## Five folds and motif calling

The preparation task writes interleaved REF/ALT sequences once. Five GPU tasks
interpret the same rows with the five trained models. They use the same
dinucleotide shuffle settings, with a shared seed for each REF/ALT pair.
Shuffled sequences can differ between alleles because their sequences differ.
Scores are averaged as signed hypothetical contributions, separately by allele.
Absolute scores and scalar log2FC values are not used for this average.

Fi-NeMo then uses the averaged contributions and the fixed discovery motifs.
Both alleles use the same library and settings. The defaults are:

| Setting | Default | Meaning |
| --- | --- | --- |
| `input_length` | 2114 | Variant model input length; must match all five models. |
| `num_backgrounds` | 20 | Shuffled DeepSHAP backgrounds per allele. |
| `random_seed` | 1234 | Background seed base. |
| `contribution_batch_size` | 64 | Rows written per attribution batch. Each allele is explained with its paired seed. |
| `shard_size` | 500 | Variants per shard. Each shard has twice as many allele rows. |
| `region_width` | 0 | Use the complete input window. An explicit width must be even and no larger than the input. |
| `mode` | `pp` | Projected motif CWMs fitted to projected contributions. `ph`, `hp`, and `hh` are also exposed. |
| `global_lambda` | 0.7 | Fi-NeMo sparsity and similarity setting; must be between 0 and 1. |
| `trim_threshold` | 0.3 | Fraction of maximum CWM importance used for motif trimming. |
| `finemo_batch_size` / `batch_size` | 256 | Variant / peak optimization batch size. |
| `max_steps` | 10000 | Maximum optimization steps. Check `dual_gap` in QC. |
| `pairing_tolerance` | 3 | Maximum difference in aligned start **and** end coordinates when pairing REF/ALT calls of the same motif and strand. |

The [Fi-NeMo documentation](https://github.com/kundajelab/Fi-NeMo) gives 0.7 as
the default lambda for accessibility data. Larger values call fewer motifs;
smaller values can call weaker motifs. These values are settings, not P-values.

Zero projected-contribution rows are not sent to the optimizer in `pp`/`hp`
mode. They have no calls, remain in the variant summary, and appear in a
separate zero-row table. For hypothetical-contribution modes, zero hypothetical
rows receive the same handling. Original row IDs are restored after subsetting.

## Variant outputs

| Output | Contents |
| --- | --- |
| `annotated_hits` | One row per allele call, with variant ID, motif ID, strand, full allele-relative coordinates, mapped reference coordinates, insertion flags, edit overlap, hit scores, and candidate TFs. |
| `motif_changes` | Paired retained calls and unpaired gained/lost calls, with REF and ALT scores and ALT-minus-REF differences. |
| `variant_summary` | One row per input variant, including those with no calls. Counts of REF/ALT calls, gains, losses, retained calls, and edit-overlapping gains/losses. |
| `allele_sequences`, `allele_rows` | Per-shard prepared HDF5 sequence/coordinate maps and exact ordered row TSVs. |
| `averaged_contributions` | Per-shard signed five-fold means. |
| `raw_call_bundles`, `optimizer_qc`, `zero_contribution_rows` | Original Fi-NeMo calls, motif metadata, optimizer results, and skipped zero rows per shard. |
| `report_bundles` | Per-shard HTML summary and TSVs. Extract the ZIP and open `report.html`. |
| Metadata and logs | Versions, hashes, settings, row counts, and logs for each stage. |

All reference hit coordinates are zero-based, with an exclusive end. Allele
coordinates refer to the full model window even when Fi-NeMo uses a central
crop. A hit with inserted bases carries their count and an insertion anchor.
A hit made entirely of inserted bases has blank reference start/end fields;
its allele coordinates and insertion anchor remain available.
Keep the per-base map when inspecting indel-spanning calls; a reference interval
can span deleted bases and is not a contiguous allele alignment.

Calls are paired independently for each variant, motif ID, and strand. The
comparison maximizes the number of valid pairs, then minimizes total aligned
boundary distance. Fixed REF/ALT coordinate order makes tied assignments
stable. Inserted-only calls cannot be paired to REF loci. A call that shifts
beyond the selected tolerance can appear as a loss
and a gain. Change the tolerance for a sensitivity analysis when needed.

`n_changes` counts gains plus losses. `n_retained` includes calls whose strength
changes and calls whose strength stays the same. Use the retained rows and
score differences in `motif_changes` to study strength changes.

For an unpaired call, the uncalled allele score is blank. Its contribution to
the call-level difference is zero. This is a bookkeeping convention, not
evidence of a zero biological effect below the threshold. Fi-NeMo scores are
not probabilities. `hit_coefficient_global` is a scale-adjusted regression
coefficient; `hit_importance` measures absolute contribution magnitude.
Neither supplies a methylation-effect direction. Positive/negative motif
groups describe accessibility contributions.

The TF match join preserves original pattern IDs, such as
`pos_patterns.pattern_0`. Multiple candidate TF names appear in one field and
do not multiply variant rows. Related TFs can share a motif. These annotations
do not prove TF occupancy or establish a methylation mechanism.

The HTML report gives motif counts and convergence totals. It does not compute
TF-MoDISco seqlet recall. Keep the original motif-discovery report for seqlet
inspection. Variant windows cannot use discovery-region recall statistics.

## Peak outputs

The peak workflow exports annotated raw hits, deduplicated hits, a genomic BED,
motif counts, optimizer QC, zero rows, an HTML bundle, and provenance. Use raw
hits for region-level analysis: the same site can occur in overlapping input
windows. Deduplicated hits select one instance per site/motif/strand and are
suited to genome-browser display.

## Resources and shard limit

Attribution and Fi-NeMo calling request one NVIDIA L4 on `g2-standard-16`,
16 CPUs, and 64 GB by default. They fail if the required framework cannot
detect a GPU. CPU flags are used only by the GitHub smoke tests. No WDL enables
them. Disk and memory settings are configurable. The other tasks use CPUs.

The final merge supports at most **128 shards**. Preparation rejects larger
requests before any GPU tasks start. Increase `shard_size` for larger inputs,
then check attribution/calling memory and disk needs. The limit permits safe
individual File quoting in WDL 1.0, which has no array map expression. It does
not limit the number of variants when a suitable shard size is selected.

Fi-NeMo conversion and optimization load a shard into memory. Peak conversion
loads the full peak contribution file. Use enough memory for these arrays;
the complete model input is retained by default. Raw call/report ZIPs can add
substantial disk use. Limit concurrent shards in Terra when planning GPU cost.

## Images and validation

The dedicated [Dockerfile](containers/Dockerfile) uses micromamba 2.3.2,
conda-forge/bioconda, Fi-NeMo 0.41, PyTorch 2.5.1, and pinned scientific packages.
The existing motif image retains TensorFlow 2.15.1 and its NumPy version for
the pinned DeepSHAP implementation. These environments are kept separate.

The Fi-NeMo image workflow runs only when `finemo_pipeline/scripts/**` changes.
The motif image workflow runs only when `motif_pipeline/scripts/**` changes.
GitHub builds run real algorithm tests before export: five-fold REF/ALT
DeepSHAP, SNV/indel preparation, separate profile scoring, motif gains/losses,
strength changes, empty rows, candidate TF joins, peak coordinates, reports,
and GPU requirements. No local Docker build is required.

A qualifying push to main publishes the tested Fi-NeMo image as
`ghcr.io/aou-multiomics-analysis/chrombpnet/finemo:<commit>` and the updated
contribution image as `ghcr.io/aou-multiomics-analysis/chrombpnet/motif-discovery:<commit>`.
Use the digest recorded in each GitHub job summary. Terra must have registry
access; copy a tested digest to an accessible registry if required.
Recipe-only, test-only, and WDL-only changes do not rebuild images; include a
related contained-script change when changing the recipe or dependencies.

Local tests cover allele maps, deterministic comparisons, missing/invalid
inputs, TF ambiguity, static workflow checks, and cloud-to-local command
rendering for every new task with optional inputs and hostile quoted paths.
GitHub also validates both workflows with miniwdl and Cromwell womtool 85.

**The complete workflows have not run on Terra.** GitHub CPU tests do not
validate managed localization, GPU execution, all released production models,
or full-run resources. No Terra analysis job has been submitted.
