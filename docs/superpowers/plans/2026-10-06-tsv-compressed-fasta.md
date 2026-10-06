# TSV manifest and compressed FASTA update

Implement in this session with the existing test and review workflow.

- Accept a TSV manifest with the exact header model_id, cell_type, model, peaks.
- Validate metadata in a CPU task. Return a headerless TSV for WDL read_tsv.
- Construct typed ModelSpec values before passing model and peaks to GPU tasks.
- Accept gzip/BGZF FASTA or plain FASTA. Decompress only inside the scoring task.
- Make genome_index optional. Use a writable copy when supplied, or create it locally.
- Keep safe command quoting, GPU settings, logging, and merge behavior.
- Add .dockstore.yml for the WDL and example input JSON using current Dockstore syntax.
- Write failing tests, implement, run unit/static checks, and run the image smoke test in Actions.
- Update examples and README. Open a follow-up PR because PR #1 has been merged.
- Do not build Docker locally or submit Terra jobs. Record Terra execution as untested.
