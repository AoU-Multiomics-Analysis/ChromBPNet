# Model fold summary and peak annotation

The user requested one summary task per model group, followed by a task that
combines outputs. Keep the existing image digest and saved scoring files usable.

## Design

- A task validates the TSV manifest and creates a group plan containing model
  names and row indices. JSON is structured output data, not a CLI argument
  wrapper. It contains no input file paths.
- A shared WDL 1.0 workflow selects typed score and peak Files using those
  indices. Score files must follow manifest row order; tasks check their labels.
- Each group task reuses the image's merge and fold summary functions for only
  that group. Inline Python adds `in_peak` using reference-allele overlap with
  BED intervals; an empty reference allele uses one base at the variant position.
  Each group must use one cell type and the same peak URI for all folds.
- A final task streams long and summary tables and uses a temporary SQLite
  database to combine wide tables without holding every model in memory.
  It checks variant identity across groups, including files with different orders.
- Preserve scoring task definitions, existing output names, and merge resource
  settings. Add separate summary memory, disk, and retry inputs. Recovery takes
  saved score files and the original model manifest; it has no scoring calls.
- Only WDL, docs, examples, Dockstore configuration, and tests change. Do not
  modify image scripts, the Dockerfile, or the image build workflow.

## Implementation and validation

1. Add failing command tests for group isolation, peak boundaries and indels,
   final combination, invalid inputs, and cloud-to-local File handling.
2. Implement the shared tasks and wire full and recovery workflows to them.
3. Update input examples, Dockstore imports, CI validation, and user instructions.
4. Run the complete unit suite, static Terra checks, miniwdl and Cromwell
   validation. Review the diff and open a PR. Confirm image build is not triggered.

No cloud jobs are authorized for validation. Report that Terra execution has not
been tested for this change.
