# GM12878 ATAC pipeline test variants

`GM12878_ATAC.test_100_variants.tsv` contains 100 synthetic single-nucleotide
substitutions for a pipeline test. These are test inputs, not observed GM12878
genotypes or validated chromatin effects.

The file is headerless and has the workflow's five columns:

1. Chromosome, with the `chr` prefix.
2. Position, using 1-based coordinates.
3. Reference allele, uppercase.
4. Alternate allele, uppercase.
5. Unique variant ID, starting with `GM12878_ATAC_TEST_`.

There are ten variants on each of chromosomes 1 through 10. Each variant is at
the summit of a distinct peak from the full input peak set used for GM12878
ATAC experiment ENCSR637XSC and model annotation ENCSR389HIH. The selected
positions are at least one million bases from the start of their chromosome.

Sources:

- [ENCODE model annotation ENCSR389HIH](https://www.encodeproject.org/annotations/ENCSR389HIH/).
- [ENCODE region archive ENCFF971WEQ](https://www.encodeproject.org/files/ENCFF971WEQ/).
  Archive member: `peaks.all_input_regions.ENCSR637XSC.bed.gz`.
  Archive MD5: `6b7ccb352e0fc2f76853144ed1cd846d`.
- [UCSC hg38 sequence API](https://genome.ucsc.edu/goldenPath/help/api.html).
  The provenance TSV gives the exact sequence URL for each site.

Generation used the first 30 peaks after position 1,000,000 on each selected
chromosome, sorted by start, end, and name. Ten unique summit sites in distinct peak intervals per
chromosome passed the sequence checks. Reference bases came from UCSC hg38;
alternate bases use the fixed mapping A to C, C to G, G to T, and T to A.
A 2,401-base window around every site contains only A, C, G, and T. This window
covers the published model input length of 2,114 bases.

Checks passed on 2026-10-06: 100 rows, five fields per row, unique variant IDs
and sites, reference allele matches, peak containment, and model window bounds.
The production workflow input parser and reference checks accepted the file.
The full scoring workflow has not been run on this test file or on Terra.

SHA256 of the variant TSV:
`7a204433c8cea23c6a018f33b4794855f8ae1a3077fb5694d67e30c39a6bb793`.

For a short pipeline test, set `max_peaks` to 100 and keep `num_shuf` at 0.
This limits the peak prediction step. These settings test execution; they do
not give a calibrated biological benchmark.

From the repository directory, upload the variant file with:

```sh
gsutil cp examples/GM12878_ATAC.test_100_variants.tsv \
  gs://your-bucket/chrombpnet/GM12878_ATAC/test_100_variants.tsv
```

Then set the workflow's `variants` input to that GCS URI. Continue to use the
full matching peak file and the hg38 reference from the download script.
