"""Validate five localized models and prepare one common, shuffled peak set."""
import argparse
import csv
import json
from pathlib import Path

import h5py
import numpy as np

from motif_io import read_file_list, readable, sha256, stage_reference, write_json


def validate_models(models, input_length, expected_folds=5):
    if len(models) != expected_folds or len({str(readable(p)) for p in models}) != expected_folds:
        raise ValueError(f'Expected {expected_folds} distinct fold model files')
    hashes = []
    for model in models:
        with h5py.File(readable(model), 'r') as handle:
            if 'model_config' not in handle.attrs:
                raise ValueError('Use standard Keras HDF5 bias-corrected ChromBPNet models')
            config = json.loads(handle.attrs['model_config'])
            layers = config.get('config', {}).get('layers', [])
            inputs = [layer['config'] for layer in layers if layer.get('class_name') == 'InputLayer']
            if len(inputs) != 1:
                raise ValueError('Each model must have one sequence input')
            shape = inputs[0].get('batch_input_shape', inputs[0].get('batch_shape'))
            if shape != [None, input_length, 4]:
                raise ValueError(f'Model input length/shape differs from configured length: {shape}')
        hashes.append(sha256(model))
    if len(set(hashes)) != expected_folds:
        raise ValueError('Fold model contents are duplicated; supply five distinct trained models')
    return hashes


def prepare_peaks(peaks, sizes, input_length, output, random_seed=1234, max_peaks=None):
    if input_length < 2 or input_length % 2 or (max_peaks is not None and max_peaks < 1):
        raise ValueError('Input length must be positive and even; max peaks must be positive')
    rows, seen = [], set()
    input_count = edge_count = 0
    with readable(peaks).open() as stream:
        for values in csv.reader(stream, delimiter='\t'):
            input_count += 1
            if not 3 <= len(values) <= 10:
                raise ValueError('Peaks must be a headerless BED3 through narrowPeak TSV')
            chrom, start, end = values[0], int(values[1]), int(values[2])
            if chrom not in sizes or not 0 <= start < end <= sizes[chrom]:
                raise ValueError('Peak chromosome or coordinates do not match the reference')
            summit = int(values[9]) if len(values) == 10 else (end - start) // 2
            if not 0 <= summit < end - start:
                raise ValueError('Peak summit must be inside the peak; -1 is not supported')
            center = start + summit
            if (chrom, center) in seen:
                raise ValueError('Duplicate peak center; remove duplicate windows before discovery')
            seen.add((chrom, center))
            if center - input_length // 2 < 0 or center + input_length // 2 > sizes[chrom]:
                edge_count += 1
                continue
            label = values[3] if len(values) > 3 else f'peak_{input_count}'
            rows.append([chrom, str(start), str(end), label, '0', '.', '0', '0', '0', str(summit)])
    if not rows:
        raise ValueError('No peak has a valid full model input window')
    # MoDISco caps seqlets in input order, so shuffle once before all fold calls.
    order = np.random.default_rng(random_seed).permutation(len(rows))
    if max_peaks is not None:
        order = order[:max_peaks]
    with Path(output).open('w') as stream:
        writer = csv.writer(stream, delimiter='\t', lineterminator='\n')
        writer.writerows(rows[int(i)] for i in order)
    return dict(input_peaks=input_count, retained_peaks=len(order),
                excluded_edge_peaks=edge_count, eligible_peaks=len(rows),
                max_peaks=max_peaks, random_seed=random_seed, input_length=input_length)


def main():
    parser = argparse.ArgumentParser()
    for name in ['models-list', 'peaks', 'genome', 'chrom-sizes', 'cell-type', 'output-dir']:
        parser.add_argument('--' + name, required=True)
    parser.add_argument('--genome-index')
    parser.add_argument('--input-length', type=int, default=2114)
    parser.add_argument('--discovery-window', type=int, default=400)
    parser.add_argument('--expected-folds', type=int, default=5)
    parser.add_argument('--random-seed', type=int, default=1234)
    parser.add_argument('--max-peaks', type=int)
    parser.add_argument('--head', choices=['counts', 'profile'], default='counts')
    args = parser.parse_args()
    if not 40 <= args.discovery_window <= args.input_length:
        raise ValueError('Discovery window must be 40 bp through the model input length')
    if not args.cell_type.strip() or any(c in args.cell_type for c in '\r\n\t'):
        raise ValueError('Cell type must be a nonempty single-line label')
    models = read_file_list(args.models_list)
    print('[prepare] Check model files and input lengths', flush=True)
    model_hashes = validate_models(models, args.input_length, args.expected_folds)
    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    reference = stage_reference(args.genome, args.genome_index, output)
    sizes = {}
    with readable(args.chrom_sizes).open() as stream:
        for row in csv.reader(stream, delimiter='\t'):
            if len(row) != 2 or row[0] in sizes or int(row[1]) < 1:
                raise ValueError('Chromosome sizes must contain unique names and positive lengths')
            sizes[row[0]] = int(row[1])
    from pyfaidx import Fasta
    with Fasta(str(reference), build_index=False, rebuild=False) as genome:
        if not sizes or any(chrom not in genome or len(genome[chrom]) != length for chrom, length in sizes.items()):
            raise ValueError('FASTA and chromosome sizes differ; check the genome build')
    meta = prepare_peaks(args.peaks, sizes, args.input_length, output / 'peaks.bed',
                         args.random_seed, args.max_peaks)
    meta.update(cell_type=args.cell_type, head=args.head, fold_count=len(models),
                model_sha256=model_hashes, original_peak_sha256=sha256(args.peaks),
                prepared_peak_sha256=sha256(output / 'peaks.bed'))
    write_json(output / 'preparation.json', meta)
    print(f"[prepare] Retained {meta['retained_peaks']} peaks; excluded {meta['excluded_edge_peaks']} edge windows", flush=True)


if __name__ == '__main__':
    main()
