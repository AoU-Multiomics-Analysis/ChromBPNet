"""Validate signed mean contributions and convert them to Fi-NeMo NPZ."""
import argparse
import csv
import json
from pathlib import Path

import h5py
import numpy as np

from finemo_io import readable, read_rows, sha256, write_json


def validate_arrays(raw, hyp):
    if raw.shape != hyp.shape or raw.ndim != 3 or raw.shape[1] != 4 or not raw.shape[0]:
        raise ValueError('Sequences and hypothetical scores require equal (N,4,L) dimensions')
    if not np.all((raw == 0) | (raw == 1)) or not np.all(raw.sum(axis=1) <= 1):
        raise ValueError('Invalid one-hot sequences')
    if not np.isfinite(hyp).all():
        raise ValueError('Contribution scores must be finite')


def validate_motifs(path, width):
    names = []
    with h5py.File(readable(path)) as handle:
        for group in ['pos_patterns', 'neg_patterns']:
            if group not in handle:
                continue
            for name, pattern in handle[group].items():
                cwm = pattern['contrib_scores'][:]
                if cwm.ndim != 2 or cwm.shape[1] != 4 or not 1 <= cwm.shape[0] <= width or not np.isfinite(cwm).all() or np.linalg.norm(cwm) == 0:
                    raise ValueError('Motifs must have finite nonzero (length,4) CWMs shorter than the region')
                names.append(group + '.' + name)
    if not names:
        raise ValueError('Motif library contains no patterns')
    return names


def prepare(contributions, regions, motifs, provenance, output, cell_type, head, kind='peaks', width=0, coordinates=None):
    provenance = json.loads(readable(provenance).read_text())
    if provenance.get('cell_type') != cell_type or provenance.get('head') != head:
        raise ValueError('Motif discovery cell type or output head differs')
    with h5py.File(readable(contributions)) as handle:
        if handle.attrs.get('cell_type') != cell_type or handle.attrs.get('head') != head:
            raise ValueError('Contribution cell type or head differs')
        if handle.attrs.get('fold_count') != 5 or handle.attrs.get('averaging') != 'arithmetic_mean_signed':
            raise ValueError('Use the signed five-fold averaged contribution file')
        raw = handle['raw/seq'][:]
        hyp = handle['shap/seq'][:]
        validate_arrays(raw, hyp)
    length = raw.shape[2]
    width = width or length
    if width < 2 or width % 2 or width > length or length % 2:
        raise ValueError('Fi-NeMo width must be positive, even, and no larger than the model input')
    start = (length - width) // 2
    names = validate_motifs(motifs, width)
    kwargs = dict(sequences=raw[:, :, start:start+width], contributions=hyp[:, :, start:start+width].astype(np.float32))
    if kind == 'variants':
        rows = read_rows(regions)
        if len(rows) != len(raw) or coordinates is None:
            raise ValueError('Allele row count differs or coordinate map is missing')
        with h5py.File(readable(coordinates)) as coords:
            if coords.attrs.get('region_sha256') != sha256(regions) or not np.array_equal(coords['raw/seq'][:], raw):
                raise ValueError('Allele coordinate map sequence/row identity differs')
            if coords['reference_positions'].shape != (len(raw), length):
                raise ValueError('Allele coordinate map dimensions differ')
    else:
        with readable(regions).open() as stream:
            peaks = list(csv.reader(stream, delimiter='\t'))
        if len(peaks) != len(raw) or any(len(row) != 10 for row in peaks):
            raise ValueError('Exact ordered narrowPeak rows must match contribution rows')
        if provenance.get('prepared_peak_sha256') != sha256(regions):
            raise ValueError('Peak order differs from motif discovery preparation')
        chroms = list(dict.fromkeys(row[0] for row in peaks))
        starts = []
        for row in peaks:
            a, b, summit = int(row[1]), int(row[2]), int(row[9])
            region_start = a + summit - width // 2
            if not 0 <= a < b or not 0 <= summit < b-a or region_start < 0:
                raise ValueError('Invalid peak coordinates or summit')
            starts.append(region_start)
        kwargs.update(chr=np.array([r[0] for r in peaks]), chr_id=np.array([chroms.index(r[0]) for r in peaks], dtype=np.uint32),
                      start=np.array(starts, dtype=np.int64), peak_id=np.arange(len(peaks), dtype=np.uint32),
                      peak_name=np.array([r[3] for r in peaks]))
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output / 'regions.npz', **kwargs)
    meta = dict(kind=kind, cell_type=cell_type, head=head, region_count=len(raw), region_width=width,
                crop_start=start, input_length=length, motif_count=len(names), motif_names=names,
                regions_sha256=sha256(regions), contributions_sha256=sha256(contributions),
                motifs_sha256=sha256(motifs), npz_sha256=sha256(output / 'regions.npz'),
                motif_provenance_sha256=sha256_file_provenance(provenance))
    write_json(output / 'metadata.json', meta)
    print(f'[prepare Fi-NeMo] rows={len(raw)} motifs={len(names)} width={width}', flush=True)
    return meta


def sha256_file_provenance(value):
    import hashlib
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    for name in ['contributions', 'regions', 'motifs', 'motif-provenance', 'cell-type', 'output-dir']:
        parser.add_argument('--' + name, required=True)
    parser.add_argument('--head', choices=['counts', 'profile'], default='counts')
    parser.add_argument('--kind', choices=['peaks', 'variants'], default='peaks')
    parser.add_argument('--region-width', type=int, default=0)
    parser.add_argument('--coordinates')
    args = parser.parse_args()
    prepare(args.contributions, args.regions, args.motifs, args.motif_provenance, args.output_dir,
            args.cell_type, args.head, args.kind, args.region_width, args.coordinates)


if __name__ == '__main__':
    main()
