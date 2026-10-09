"""Average signed hypothetical contributions without loading all peaks at once."""
import argparse
from contextlib import ExitStack
from pathlib import Path
import shutil

import h5py
import numpy as np

from motif_io import (contribution_shape, create_contribution_datasets,
                      read_file_list, readable, sha256, write_json)


def average_contributions(files, regions, output, output_regions, expected_folds=5, chunk_rows=128):
    if len(files) != expected_folds or len(regions) != expected_folds or chunk_rows < 1:
        raise ValueError(f'Expected {expected_folds} fold score and region files; chunk size must be positive')
    if len({str(readable(p)) for p in files}) != expected_folds:
        raise ValueError('Duplicate fold score file')
    reference_regions = readable(regions[0]).read_bytes()
    if any(readable(p).read_bytes() != reference_regions for p in regions[1:]):
        raise ValueError('Fold region coordinates or order differ')
    temporary = Path(str(output) + '.partial')
    try:
        with ExitStack() as stack:
            handles = [stack.enter_context(h5py.File(readable(p), 'r')) for p in files]
            shape = contribution_shape(handles[0])
            if any(contribution_shape(h) != shape for h in handles[1:]):
                raise ValueError('Fold contribution dimensions differ')
            if len(reference_regions.splitlines()) != shape[0]:
                raise ValueError('Region count and contribution array row count differ')
            for key in ['head', 'cell_type']:
                if key not in handles[0].attrs or any(h.attrs.get(key) != handles[0].attrs[key] for h in handles[1:]):
                    raise ValueError(f'Fold {key} values differ or are missing')
            if handles[0].attrs['head'] not in ['counts', 'profile']:
                raise ValueError('Unknown contribution head')
            fold_indices = [h.attrs.get('fold_index', -1) for h in handles]
            if sorted(fold_indices) != list(range(expected_folds)):
                raise ValueError('Missing or duplicate fold indices')
            target = stack.enter_context(h5py.File(temporary, 'w'))
            create_contribution_datasets(target, shape, chunk_rows)
            target.attrs.update(head=handles[0].attrs['head'], cell_type=handles[0].attrs['cell_type'],
                                fold_count=expected_folds, base_order='ACGT', averaging='arithmetic_mean_signed')
            for start in range(0, shape[0], chunk_rows):
                stop = min(start + chunk_rows, shape[0])
                raw = handles[0]['raw/seq'][start:stop]
                total = np.zeros(raw.shape, dtype=np.float64)
                for handle in handles:
                    other = handle['raw/seq'][start:stop]
                    if not np.all((other == 0) | (other == 1)) or not np.all(other.sum(axis=1) <= 1):
                        raise ValueError('Invalid one-hot sequence')
                    if not np.array_equal(raw, other):
                        raise ValueError('Fold sequence identity or order differs')
                    scores = handle['shap/seq'][start:stop]
                    if not np.isfinite(scores).all():
                        raise ValueError('Fold contribution scores must be finite')
                    total += scores
                mean = (total / expected_folds).astype(np.float32)
                target['raw/seq'][start:stop] = raw
                target['shap/seq'][start:stop] = mean
                target['projected_shap/seq'][start:stop] = mean * raw
                print(f'[average] Processed {stop}/{shape[0]} regions', flush=True)
            result = dict(fold_count=expected_folds, region_count=shape[0], input_length=shape[2],
                          head=str(handles[0].attrs['head']), cell_type=str(handles[0].attrs['cell_type']),
                          region_sha256=sha256(regions[0]), base_order='ACGT',
                          averaging='arithmetic mean of signed hypothetical contributions')
        temporary.replace(output)
        shutil.copyfile(regions[0], output_regions)
        return result
    except Exception:
        temporary.unlink(missing_ok=True)
        raise


def main():
    parser = argparse.ArgumentParser()
    for name in ['scores-list', 'regions-list', 'output', 'output-regions', 'metadata']:
        parser.add_argument('--' + name, required=True)
    parser.add_argument('--expected-folds', type=int, default=5)
    parser.add_argument('--chunk-rows', type=int, default=128)
    args = parser.parse_args()
    result = average_contributions(read_file_list(args.scores_list), read_file_list(args.regions_list),
                                   args.output, args.output_regions, args.expected_folds, args.chunk_rows)
    write_json(args.metadata, result)
    print('[average] Fold averaging complete', flush=True)


if __name__ == '__main__':
    main()
