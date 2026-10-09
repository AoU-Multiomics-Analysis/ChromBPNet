# TF-MoDISco workflow implementation plan

Build the workflow approved in chat. Target WDL 1.0 and Terra managed Cromwell.
Use one cell type, one output (`counts` by default), and five bias-corrected
HDF5 models per submission. Do not submit analysis jobs.

1. Add tests for common peak preparation, signed fold averaging, rejected
   sequence/coordinate mismatches, finite scores, and missing folds.
2. Implement small named-argument scripts under `motif_pipeline/scripts`.
   Prepare and shuffle the common peak set once. Score it in five GPU tasks,
   using the pinned ChromBPNet contribution functions and the same shuffle
   seed/settings. Average hypothetical contributions in bounded row chunks.
3. Add a WDL with explicit File inputs, task-local newline file lists, safe
   shell quoting, logs, CPU motif discovery, and a separate report task.
   Keep counts and profile in separate submissions. Export Fi-NeMo-compatible
   TF-MoDISco HDF5, averaged contribution HDF5, region coordinates, MEME
   motifs, Tomtom results, report assets, and provenance.
4. Add a dedicated micromamba image with pinned packages/source. Add a
   GitHub Actions build triggered only by changes to scripts in that image.
   Run real contribution scoring, motif discovery, and motif matching on
   synthetic data in the Linux image; do not build Docker locally.
5. Run the whole repository test suite, miniwdl/static Terra checks and
   command-localization regression tests. Push a branch, create a reviewable
   PR, and inspect GitHub checks. State which checks ran and that a complete
   Terra run remains untested.

Review focus: an unresolved cloud URI must fail before computation; File
arrays must remain typed until command rendering; averaging must preserve
negative scores and exact peak order; no incompatible model width or missing
fold may silently change the denominator; a discovery result with no motifs
must report that condition clearly.
