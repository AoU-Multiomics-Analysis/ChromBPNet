"""File localization, reference staging, and output helpers for motif tasks."""
import gzip
import hashlib
import json
import shutil
from pathlib import Path

import h5py


def readable(value):
    value = str(value)
    if '://' in value or any(c in value for c in '\r\n'):
        raise ValueError(f'File localization error: unresolved URI or invalid path: {value}')
    path = Path(value)
    if not path.is_file():
        raise ValueError(f'File localization error: input is not readable: {value}')
    with path.open('rb') as stream:
        stream.read(1)
    return path.resolve()


def read_file_list(value):
    lines = readable(value).read_text().splitlines()
    if not lines or any(not line for line in lines):
        raise ValueError('File list is empty or contains an empty path')
    return [readable(line) for line in lines]


def sha256(path):
    digest = hashlib.sha256()
    with readable(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def stage_reference(genome, index, output):
    from pyfaidx import Fasta
    genome = readable(genome)
    target = Path(output) / 'reference.fa'
    with genome.open('rb') as stream:
        compressed = stream.read(2) == b'\x1f\x8b'
    if compressed:
        print('[reference] Decompress localized reference', flush=True)
        with gzip.open(genome, 'rb') as source, target.open('wb') as destination:
            shutil.copyfileobj(source, destination, 1024 * 1024)
    else:
        shutil.copyfile(genome, target)
    if index is not None:
        shutil.copyfile(readable(index), str(target) + '.fai')
    with Fasta(str(target), rebuild=index is None) as reference:
        if not len(reference.keys()):
            raise ValueError('Reference FASTA has no sequences')
    return target


def create_contribution_datasets(handle, shape, chunk_rows=128):
    chunks = (min(chunk_rows, shape[0]), shape[1], shape[2])
    for key, dtype in [('raw/seq', 'int8'), ('shap/seq', 'float32'),
                       ('projected_shap/seq', 'float32')]:
        handle.create_dataset(key, shape=shape, dtype=dtype, chunks=chunks,
                              compression='gzip', compression_opts=4, shuffle=True)


def contribution_shape(handle):
    if any(key not in handle for key in ['raw/seq', 'shap/seq', 'projected_shap/seq']):
        raise ValueError('Contribution HDF5 must contain raw, shap, and projected_shap seq datasets')
    shape = handle['raw/seq'].shape
    if len(shape) != 3 or shape[0] < 1 or shape[1] != 4 or shape[2] < 2:
        raise ValueError('Contribution array must have shape (regions, 4, sequence_length)')
    if any(handle[key].shape != shape for key in ['shap/seq', 'projected_shap/seq']):
        raise ValueError('Contribution array dimensions differ')
    return shape
