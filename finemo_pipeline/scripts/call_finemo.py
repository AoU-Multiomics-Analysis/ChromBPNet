"""Run pinned Fi-NeMo with a GPU requirement and explicit empty-row handling."""
import argparse
import importlib.metadata
import json
import shutil
from pathlib import Path

import numpy as np

from finemo_io import readable, sha256, read_table, write_table, write_json


def main():
    parser = argparse.ArgumentParser()
    for name in ['regions', 'motifs', 'preparation', 'output-dir']:
        parser.add_argument('--' + name, required=True)
    parser.add_argument('--mode', choices=['pp', 'ph', 'hp', 'hh'], default='pp')
    parser.add_argument('--global-lambda', type=float, default=0.7)
    parser.add_argument('--trim-threshold', type=float, default=0.3)
    parser.add_argument('--batch-size', type=int, default=256)
    parser.add_argument('--max-steps', type=int, default=10000)
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--allow-cpu', action='store_true', help='GitHub smoke test only')
    args = parser.parse_args()
    args.regions, args.motifs, args.preparation = map(readable, [args.regions, args.motifs, args.preparation])
    meta = json.loads(args.preparation.read_text())
    if sha256(args.regions) != meta['npz_sha256'] or sha256(args.motifs) != meta['motifs_sha256']:
        raise ValueError('Fi-NeMo input or motif identity differs from preparation')
    if not 0 < args.global_lambda < 1 or not 0 <= args.trim_threshold < 1 or min(args.batch_size, args.max_steps, args.threads) < 1:
        raise ValueError('Invalid Fi-NeMo sensitivity, trimming, batch, steps, or threads setting')
    import torch
    torch.set_num_threads(args.threads)
    if not torch.cuda.is_available() and not args.allow_cpu:
        raise ValueError('GPU is required, but PyTorch cannot detect one')
    from finemo import data_io, main as finemo
    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    with np.load(args.regions, allow_pickle=False) as data:
        arrays = {key: data[key] for key in data.files}
    # A zero contribution row has no motifs and cannot be normalized by Fi-NeMo.
    signal = arrays['contributions']
    if args.mode in ['pp', 'hp']:
        signal = signal * arrays['sequences']
    active = np.flatnonzero(np.any(signal != 0, axis=(1, 2)))
    if len(active):
        np.savez_compressed(output / 'active.npz', **{key: value[active] for key, value in arrays.items()})
        finemo.call_hits(str(output / 'active.npz'), None, str(args.motifs), None, None, None, None,
                        str(output), cwm_trim_threshold_default=args.trim_threshold,
                        lambda_default=args.global_lambda, batch_size=args.batch_size, max_steps=args.max_steps,
                        mode=args.mode, device='cpu' if args.allow_cpu and not torch.cuda.is_available() else None)
        # Without absolute coordinates Fi-NeMo assigns row IDs after subsetting.
        # Restore the original REF/ALT row index before comparison.
        if 'chr' not in arrays:
            for name in ['hits.tsv', 'peaks_qc.tsv']:
                rows = read_table(output / name)
                for row in rows:
                    row['peak_id'] = str(active[int(row['peak_id'])])
                with (output / name).open() as stream:
                    fields = stream.readline().rstrip('\n').split('\t')
                write_table(output / name, rows, fields)
        (output / 'active.npz').unlink()
    else:
        write_table(output / 'hits.tsv', [], list(data_io.HITS_DTYPES))
        write_table(output / 'hits_unique.tsv', [], list(data_io.HITS_DTYPES))
        (output / 'hits.bed').write_text('')
        write_table(output / 'peaks_qc.tsv', [], ['peak_id', 'nll', 'dual_gap', 'num_steps', 'step_size', 'global_scale'])
    active_set = set(active.tolist())
    write_table(output / 'zero_contribution_rows.tsv', [{'region_index': i} for i in range(len(arrays['sequences'])) if i not in active_set], ['region_index'])
    write_json(output / 'metadata.json', dict(meta, finemo=importlib.metadata.version('finemo'),
               torch=torch.__version__, mode=args.mode, global_lambda=args.global_lambda,
               trim_threshold=args.trim_threshold, max_steps=args.max_steps,
               batch_size=args.batch_size, gpu_available=torch.cuda.is_available(),
               allow_cpu=args.allow_cpu, zero_contribution_rows=len(arrays['sequences']) - len(active)))
    archive = shutil.make_archive(str(output.parent / 'raw_calls'), 'zip', output)
    Path(archive).replace(output / 'raw_calls.zip')
    print(f'[Fi-NeMo] complete active_rows={len(active)} zero_rows={len(arrays["sequences"])-len(active)}', flush=True)


if __name__ == '__main__':
    main()
