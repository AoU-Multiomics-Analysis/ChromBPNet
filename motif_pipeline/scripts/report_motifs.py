"""Generate a portable report and retain significant candidate TF matches."""
import argparse
import csv
import importlib.metadata
import os
from pathlib import Path
import shutil
import subprocess

import h5py
import numpy as np

from motif_io import readable, sha256, write_json


def export_trimmed_motifs(source, output, threshold=0.3):
    inventory = []
    with h5py.File(readable(source), 'r') as handle, Path(output).open('w') as stream:
        stream.write('MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\nBackground letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n\n')
        for group in ['pos_patterns', 'neg_patterns']:
            if group not in handle:
                continue
            for name, pattern in handle[group].items():
                cwm = pattern['contrib_scores'][:]
                pfm = pattern['sequence'][:]
                if cwm.shape != pfm.shape or pfm.ndim != 2 or pfm.shape[1] != 4 or not np.isfinite(cwm).all():
                    raise ValueError('Invalid motif matrix')
                importance = np.abs(cwm).sum(axis=1)
                positions = np.flatnonzero(importance >= threshold * importance.max())
                if importance.max() <= 0 or len(positions) == 0:
                    raise ValueError('Motif has no nonzero contribution positions')
                start, end = int(positions[0]), int(positions[-1]) + 1
                matrix = pfm[start:end].astype(np.float64) + 1e-6
                if not np.isfinite(matrix).all() or np.any(matrix < 0):
                    raise ValueError('Invalid motif sequence probabilities')
                matrix /= matrix.sum(axis=1, keepdims=True)
                tag = f'{group}.{name}'
                nseqlets = int(pattern['seqlets/n_seqlets'][()][0])
                inventory.append(dict(pattern=tag, h5_group=f'{group}/{name}', seqlet_count=nseqlets,
                                      trim_start=start, trim_end=end, width=end - start))
                stream.write(f'MOTIF {tag}\nletter-probability matrix: alength= 4 w= {len(matrix)} nsites= {nseqlets}\n')
                np.savetxt(stream, matrix, fmt='%.8f')
                stream.write('\n')
    if not inventory:
        raise ValueError('No motif patterns are available for reporting')
    return inventory


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--motifs', required=True)
    parser.add_argument('--motif-database', required=True)
    parser.add_argument('--output-dir', required=True)
    parser.add_argument('--match-qvalue', type=float, default=0.05)
    parser.add_argument('--n-matches', type=int, default=5)
    args = parser.parse_args()
    if not 0 < args.match_qvalue <= 1 or args.n_matches < 1:
        raise ValueError('Match q-value must be in (0,1]; match count must be positive')
    source, database = readable(args.motifs), readable(args.motif_database)
    output = Path(args.output_dir).resolve()
    output.mkdir(parents=True, exist_ok=True)
    # Upstream reporting uses shell commands internally. Use fixed safe names
    # inside its working directory, regardless of the localized input names.
    shutil.copyfile(source, output / 'modisco_results.h5')
    shutil.copyfile(database, output / 'known_motifs.meme')
    inventory = export_trimmed_motifs(source, output / 'discovered_motifs.meme')
    print('[report] Match trimmed sequence motifs with MEME Tomtom', flush=True)
    subprocess.run(['tomtom', '-oc', 'tomtom', '-dist', 'pearson', '-min-overlap', '5',
                    '-thresh', str(args.match_qvalue), 'discovered_motifs.meme', 'known_motifs.meme'],
                   cwd=output, check=True)
    labels = {}
    for line in database.read_text().splitlines():
        fields = line.split()
        if fields and fields[0] == 'MOTIF' and len(fields) >= 2:
            labels[fields[1]] = ' '.join(fields[2:]) or fields[1]
    matches = []
    with (output / 'tomtom/tomtom.tsv').open() as stream:
        for row in csv.DictReader((line for line in stream if line.strip() and not line.startswith('#')), delimiter='\t'):
            if float(row['q-value']) <= args.match_qvalue:
                matches.append(dict(pattern=row['Query_ID'], target_motif=row['Target_ID'],
                                    candidate_tf=labels.get(row['Target_ID'], row['Target_ID']),
                                    q_value=row['q-value'], p_value=row['p-value'],
                                    overlap=row['Overlap'], orientation=row['Orientation']))
    with (output / 'candidate_tf_matches.tsv').open('w') as stream:
        fields = ['pattern', 'target_motif', 'candidate_tf', 'q_value', 'p_value', 'overlap', 'orientation']
        writer = csv.DictWriter(stream, fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(matches)
    with (output / 'motif_inventory.tsv').open('w') as stream:
        writer = csv.DictWriter(stream, list(inventory[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(inventory)
    print('[report] Generate TF-MoDISco report and motif logos', flush=True)
    subprocess.run(['modisco', 'report', '-i', 'modisco_results.h5', '-o', 'report',
                    '-s', './', '-m', 'known_motifs.meme', '-n', str(args.n_matches)],
                   cwd=output, env=dict(os.environ, MPLBACKEND='Agg'), check=True)
    write_json(output / 'report_metadata.json', dict(
        modisco=importlib.metadata.version('modisco'), motif_database_sha256=sha256(database),
        pattern_count=len(inventory), reference_motif_count=len(labels),
        significant_match_count=len(matches), match_qvalue=args.match_qvalue,
        annotation='candidate TF motif similarity; not proof of TF occupancy',
        tomtom_distance='pearson', min_overlap=5, trim_threshold=0.3, plotting_backend='Agg'))
    shutil.make_archive(str(output / 'report_bundle'), 'zip', root_dir=output,
                        base_dir='report')
    print(f'[report] Complete: {len(matches)} candidate motif matches', flush=True)


if __name__ == '__main__':
    main()
