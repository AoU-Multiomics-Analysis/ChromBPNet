"""Run the released TF-MoDISco CLI on averaged hypothetical scores."""
import argparse
import importlib.metadata
from pathlib import Path
import runpy
import shutil
import sys

import h5py
import numpy as np

from motif_io import contribution_shape, readable, write_json


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--contributions', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--metadata', required=True)
    parser.add_argument('--window', type=int, default=400)
    parser.add_argument('--max-seqlets', type=int, default=1000000)
    parser.add_argument('--n-leiden', type=int, default=2)
    parser.add_argument('--random-seed', type=int, default=1234)
    args = parser.parse_args()
    source = readable(args.contributions)
    with h5py.File(source, 'r') as handle:
        shape = contribution_shape(handle)
    if args.window < 40 or args.window > shape[2] or min(args.max_seqlets, args.n_leiden) < 1:
        raise ValueError('Discovery window must be 40 bp through input length; seqlet and Leiden limits must be positive')
    binary = shutil.which('modisco')
    if binary is None:
        raise ValueError('TF-MoDISco executable is missing from the image')
    np.random.seed(args.random_seed)
    sys.argv = [binary, 'motifs', '-i', str(source), '-n', str(args.max_seqlets),
                '-w', str(args.window), '-l', str(args.n_leiden), '-o', args.output, '-v']
    print('[modisco] Start motif discovery on fold-averaged contributions', flush=True)
    runpy.run_path(binary, run_name='__main__')
    with h5py.File(args.output, 'r') as handle:
        patterns = [f'{group}.{name}' for group in ['pos_patterns', 'neg_patterns']
                    if group in handle for name in handle[group]]
    versions = {name: importlib.metadata.version(name) for name in ['modisco', 'igraph', 'leidenalg']}
    print(f'[modisco] Algorithm versions: {versions}', flush=True)
    write_json(args.metadata, dict(versions=versions, modisco=versions['modisco'],
                                  window=args.window, max_seqlets=args.max_seqlets,
                                  n_leiden=args.n_leiden, random_seed=args.random_seed,
                                  pattern_count=len(patterns), patterns=patterns))
    if not patterns:
        raise ValueError('No motifs discovered; inspect the contributions and peak set before changing discovery settings')
    print(f'[modisco] Discovered {len(patterns)} patterns', flush=True)


if __name__ == '__main__':
    main()
