"""Validate five models and write variant allele shards for GPU scoring."""
import argparse
import csv
import json
from pathlib import Path

import h5py

from motif_io import readable, read_file_list, sha256, stage_reference, write_json
from prepare_motif_inputs import validate_models
from variant_sequences import allele_windows, edit_bounds, one_hot

FIELDS = ['chr', 'pos', 'allele1', 'allele2', 'variant_id']
ROW_FIELDS = ['region_index', 'variant_id', 'allele', 'chr', 'pos', 'ref', 'alt',
              'edit_start', 'edit_end', 'window_start', 'insertion_anchor']


def read_variants(path):
    rows, seen = [], set()
    with readable(path).open() as stream:
        for values in csv.reader(stream, delimiter='\t'):
            if len(values) != 5 or any(not v for v in values):
                raise ValueError('Variants require five nonempty headerless TSV fields')
            row = dict(zip(FIELDS, values))
            if not row['pos'].isdigit() or int(row['pos']) < 1 or row['variant_id'] in seen:
                raise ValueError('Variant positions must be positive and IDs must be unique')
            seen.add(row['variant_id'])
            rows.append(row)
    if not rows:
        raise ValueError('Variant list is empty')
    return rows


def prepare_variants(path, genome, width, shard_size, output, cell_type, head, model_hashes):
    if shard_size < 1:
        raise ValueError('Shard size must be positive')
    rows = read_variants(path)
    for shard, begin in enumerate(range(0, len(rows), shard_size)):
        sequences, maps, records = [], [], []
        for row in rows[begin:begin + shard_size]:
            seqs, positions, _ = allele_windows(row, genome, width)
            prefix, suffix = edit_bounds(row['allele1'], row['allele2'])
            for allele, seq, coord in zip(['REF', 'ALT'], seqs, positions):
                edited = row['allele1'] if allele == 'REF' else row['allele2']
                records.append([len(records), row['variant_id'], allele, row['chr'], int(row['pos']),
                                row['allele1'], row['allele2'], width // 2 + prefix,
                                width // 2 + len(edited) - suffix,
                                int(row['pos']) - 1 - width // 2,
                                max(int(row['pos']) - 1, int(row['pos']) + prefix - 2)])
                sequences.append(seq)
                maps.append(coord)
        stem = Path(output) / f'shard_{shard:06d}'
        with stem.with_suffix('.tsv').open('w') as stream:
            csv.writer(stream, delimiter='\t', lineterminator='\n').writerows(records)
        with h5py.File(stem.with_suffix('.h5'), 'w') as handle:
            handle.create_dataset('raw/seq', data=one_hot(sequences), compression='gzip')
            handle.create_dataset('reference_positions', data=maps, compression='gzip')
            handle.attrs.update(cell_type=cell_type, head=head, row_count=len(records), base_order='ACGT',
                                region_sha256=sha256(stem.with_suffix('.tsv')))
        write_json(stem.with_suffix('.json'), dict(cell_type=cell_type, head=head, input_length=width,
                   fold_count=5, model_sha256=model_hashes, row_count=len(records),
                   region_sha256=sha256(stem.with_suffix('.tsv')),
                   sequences_sha256=sha256(stem.with_suffix('.h5'))))
        print(f'[prepare variants] shard={shard} variants={len(records)//2}', flush=True)
    return len(rows)


def main():
    parser = argparse.ArgumentParser()
    for name in ['models-list', 'variants', 'genome', 'chrom-sizes', 'cell-type', 'motif-provenance', 'output-dir']:
        parser.add_argument('--' + name, required=True)
    parser.add_argument('--genome-index')
    parser.add_argument('--head', choices=['counts', 'profile'], default='counts')
    parser.add_argument('--input-length', type=int, default=2114)
    parser.add_argument('--shard-size', type=int, default=500)
    args = parser.parse_args()
    if not args.cell_type.strip() or any(c in args.cell_type for c in '\t\r\n'):
        raise ValueError('Cell type must be a nonempty single-line label')
    hashes = validate_models(read_file_list(args.models_list), args.input_length)
    provenance = json.loads(readable(args.motif_provenance).read_text())
    if provenance.get('cell_type') != args.cell_type or provenance.get('head') != args.head or provenance.get('model_sha256') != hashes:
        raise ValueError('Variant models, cell type, or head differ from motif discovery preparation')
    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    reference = stage_reference(args.genome, args.genome_index, output)
    from pyfaidx import Fasta
    with readable(args.chrom_sizes).open() as stream:
        sizes = list(csv.reader(stream, delimiter='\t'))
    if not sizes or any(len(r) != 2 or not r[1].isdigit() or int(r[1]) < 1 for r in sizes) or len({r[0] for r in sizes}) != len(sizes):
        raise ValueError('Chromosome sizes require unique names and positive lengths')
    with Fasta(str(reference), build_index=False, rebuild=False) as genome:
        if any(chrom not in genome or len(genome[chrom]) != int(length) for chrom, length in sizes):
            raise ValueError('Chromosome sizes and reference FASTA differ')
        rows = read_variants(args.variants)
        if any(row['chr'] not in {r[0] for r in sizes} for row in rows):
            raise ValueError('Variant chromosome is absent from chromosome sizes')
        count = prepare_variants(args.variants, genome, args.input_length, args.shard_size, output,
                                 args.cell_type, args.head, hashes)
    write_json(output / 'metadata.json', dict(variant_count=count, cell_type=args.cell_type,
               head=args.head, model_sha256=hashes, variants_sha256=sha256(args.variants),
               input_length=args.input_length, shard_size=args.shard_size,
               coordinate_convention='zero-based reference positions; -1 marks inserted bases'))


if __name__ == '__main__':
    main()
