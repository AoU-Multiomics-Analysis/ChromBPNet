# Fi-NeMo workflow design for Terra

Status: proposed design for user review. Implementation has not started.

## Purpose

Use the cell-type motif library from the five-fold TF-MoDISco workflow to
identify motif instances and candidate TF effects of methylation-QTL variants.
The recommended first release includes REF/ALT variant annotation. A smaller
release can start with motif calls from existing contribution files.

## Two entry points

1. `call_motifs.wdl`: accept the averaged ChromBPNet contribution HDF5,
   its exact ordered region BED, the TF-MoDISco motif HDF5, and an optional
   candidate-TF match table. Validate row identity, dimensions, finite scores,
   cell type, and output head. Convert contributions to Fi-NeMo input, call
   motif instances, and export annotated hits, per-region QC, and a report.
2. `annotate_variant_motifs.wdl`: accept five ordered ChromBPNet model Files,
   the reference FASTA and optional index, chromosome sizes, the existing
   headerless variant TSV, and the corresponding motif library. Prepare
   REF/ALT sequences once. Compute contributions for each allele in each
   fold. Average signed hypothetical scores separately for REF and ALT.
   Call Fi-NeMo on each allele, then compare motif instances per variant.

Run one cell type and one output head per submission. Use `counts` by default;
run `profile` separately. Use matching cell-type models and motif libraries.
Do not require a variant to overlap an ATAC peak. This permits candidate
motif gains outside the existing peak set.

## Alleles and coordinates

Reuse the current TSV fields: chromosome, 1-based position, REF, ALT, and
variant ID. Check REF against the selected genome. Reject duplicate IDs,
invalid alleles, and incomplete windows with a clear error.

Handle substitutions and indels with a per-base coordinate map for each
allele. Inserted bases must retain an allele-relative position and an insertion
anchor; they must not receive fabricated reference-genome coordinates.
Keep all motif calls linked to their variant and allele. Do not deduplicate
hits across variant windows before the REF/ALT comparison.

Export all allele calls. Pair calls using motif ID, strand, and aligned locus,
with a documented, configurable positional tolerance. Record unpaired REF
calls as candidate losses and unpaired ALT calls as candidate gains. For paired
calls, report ALT minus REF global hit coefficient, importance, and similarity.
Separate calls that overlap the edited allele from other calls in the context
window. Retain variants with no calls in the variant summary.

## Fi-NeMo settings and annotations

Pin Fi-NeMo 0.41. Start with upstream defaults: projected motif/contribution
mode (`pp`), CWM trim threshold 0.3, and global lambda 0.7. Expose these as
typed settings. Do not use a hit coefficient as a probability or a P-value.
Motif gain/loss depends on the call threshold and is a model prediction.

Keep original pattern IDs such as `pos_patterns.pattern_0`. Join the existing
Tomtom candidate-TF table without replacing pattern IDs or discarding multiple
TF matches. TF similarity does not establish occupancy or methylation causality.

For existing discovery regions, recall statistics require the exact discovery
row order and discovery window. Disable recall for REF/ALT variant sequences.
Record checksums, package versions, parameters, row counts, and task logs.

## Images and WDL

Use WDL 1.0 and Terra managed Cromwell. Supply every opened input as File,
File?, or Array[File]. Preserve File types until command rendering. Use safely
quoted named CLI arguments and task-local newline lists for file arrays.
Do not use workflow-scope file writers or JSON argument wrappers. Check local
input readability and report unresolved cloud URIs before computation.

Reuse the existing contribution image for DeepSHAP. Add a separate Fi-NeMo
image based on pinned micromamba, with conda-forge/bioconda packages and pinned
PyTorch CUDA support. Only scoring and motif-call tasks request GPUs; CPU
preparation, averaging, and summary tasks remain separate. Make memory, disk,
batch size, and variant shard size configurable.

Build and test images in GitHub Actions. Rebuild only when scripts contained
in that image change. Run real synthetic motif calling and REF/ALT comparison
tests during the image build. Publish tested images on qualifying main-branch
pushes and record the digest. Do not build Docker locally for the smoke test.

## Required checks

- Known planted motif retained, gained, and lost across REF/ALT sequences.
- Real CPU Fi-NeMo smoke test; production tasks require a detected GPU.
- Signed averaging, exact sequence/row identity, and missing-fold rejection.
- Indel coordinate maps and nearby calls that shift after an insertion.
- Empty hit sets, zero-effect variants, reverse-strand calls, and ambiguous TFs.
- Cloud-to-local path rendering for every task, including optional inputs.
- Safe quoting of file paths and rejection of line breaks in file lists.
- Static workflow-scope file-write checks, miniwdl, and Cromwell validation.

Deliver the code, example Terra inputs, documentation, and passed GitHub checks
in a draft PR. Do not merge or submit a Terra analysis job without user approval.
State explicitly that full Terra execution remains untested until it runs.

## Primary references

- https://github.com/kundajelab/Fi-NeMo
- https://kundajelab.github.io/Fi-NeMo/finemo.html
- https://pypi.org/project/finemo/0.41/
- Existing `motif_pipeline/README.md` and `workflows/discover_motifs.wdl`.
