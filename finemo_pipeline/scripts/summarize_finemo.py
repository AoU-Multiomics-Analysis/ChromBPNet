"""Annotate motif hits and compare allele calls in aligned reference coordinates."""
import argparse
import html
import json
import math
import shutil
from collections import Counter
from pathlib import Path

import h5py

from finemo_io import readable, read_rows, read_table, write_table, tf_matches, sha256, write_json

METRICS = ['hit_coefficient_global', 'hit_importance', 'hit_similarity']
CALL_FIELDS = ['variant_id', 'allele', 'chr', 'pos', 'ref', 'alt', 'region_index', 'motif_name',
               'strand', 'allele_start', 'allele_end', 'reference_start', 'reference_end',
               'inserted_base_count', 'insertion_anchor', 'overlaps_edit', 'candidate_tfs'] + METRICS
CHANGE_FIELDS = ['variant_id', 'chr', 'pos', 'ref', 'alt', 'motif_name', 'strand', 'status',
                 'ref_start', 'alt_start', 'reference_start', 'reference_end', 'overlaps_edit', 'candidate_tfs']
for metric in METRICS:
    CHANGE_FIELDS += ['ref_' + metric, 'alt_' + metric, 'delta_' + metric]
VARIANT_FIELDS = ['variant_id', 'chr', 'pos', 'ref', 'alt', 'n_ref_hits', 'n_alt_hits', 'n_changes',
                  'n_gained', 'n_lost', 'n_retained', 'n_edit_overlapping_changes']


def compare_variants(rows, hits, maps, matches, tolerance=3, crop_start=0):
    if tolerance < 0 or maps.shape[0] != len(rows):
        raise ValueError('Invalid pairing tolerance or coordinate row count')
    variants, index = {}, {}
    for i, row in enumerate(rows):
        if int(row['region_index']) != i or row['allele'] not in ['REF', 'ALT']:
            raise ValueError('Invalid allele identity or row order')
        key = row['variant_id']
        pair = variants.setdefault(key, {})
        if row['allele'] in pair:
            raise ValueError('Duplicate allele row')
        pair[row['allele']] = row
        index[i] = row
    for pair in variants.values():
        if set(pair) != {'REF', 'ALT'} or any(pair['REF'][k] != pair['ALT'][k] for k in ['chr', 'pos', 'ref', 'alt']):
            raise ValueError('Every variant requires matching REF and ALT rows')
    calls, groups = [], {}
    for hit in hits:
        row_id = int(hit['peak_id'])
        if row_id not in index:
            raise ValueError('Hit references an unknown allele row')
        row = index[row_id]
        a, b = int(hit['start']) + crop_start, int(hit['end']) + crop_start
        if not 0 <= a < b <= maps.shape[1] or hit['strand'] not in ['+', '-']:
            raise ValueError('Invalid motif hit coordinates or strand')
        values = {k: float(hit[k]) for k in METRICS}
        if not all(math.isfinite(v) for v in values.values()):
            raise ValueError('Hit scores must be finite')
        positions = maps[row_id, a:b]
        mapped = positions[positions >= 0]
        anchor = int(row['insertion_anchor'])
        ref_start = int(mapped.min()) if len(mapped) else anchor
        ref_end = int(mapped.max()) + 1 if len(mapped) else anchor + 1
        x, y = int(row['edit_start']), int(row['edit_end'])
        overlap = (a < y and b > x) if x != y else a <= x < b
        call = dict(variant_id=row['variant_id'], allele=row['allele'], chr=row['chr'], pos=int(row['pos']),
                    ref=row['ref'], alt=row['alt'], region_index=row_id, motif_name=hit['motif_name'], strand=hit['strand'],
                    allele_start=a, allele_end=b, reference_start=ref_start, reference_end=ref_end,
                    inserted_base_count=int((positions < 0).sum()), insertion_anchor=anchor if (positions < 0).any() else '',
                    overlaps_edit=overlap, candidate_tfs=';'.join(matches.get(hit['motif_name'], [])), **values)
        calls.append(call)
        groups.setdefault((row['variant_id'], hit['motif_name'], hit['strand']), {'REF': [], 'ALT': []})[row['allele']].append(call)
    changes = []
    for (variant_id, motif, strand), pair in sorted(groups.items()):
        refs = sorted(pair['REF'], key=lambda r: (r['reference_start'], r['allele_start']))
        alts = sorted(pair['ALT'], key=lambda r: (r['reference_start'], r['allele_start']))
        # Greedy minimum boundary-distance matching, with explicit deterministic ties.
        candidates = []
        for i, ref in enumerate(refs):
            for j, alt in enumerate(alts):
                delta = max(abs(ref['reference_start'] - alt['reference_start']), abs(ref['reference_end'] - alt['reference_end']))
                if delta <= tolerance:
                    candidates.append((delta, i, j))
        used_ref, used_alt, paired = set(), set(), []
        for _, i, j in sorted(candidates):
            if i not in used_ref and j not in used_alt:
                used_ref.add(i)
                used_alt.add(j)
                paired.append((refs[i], alts[j]))
        paired += [(ref, None) for i, ref in enumerate(refs) if i not in used_ref]
        paired += [(None, alt) for j, alt in enumerate(alts) if j not in used_alt]
        for ref, alt in paired:
            source = ref or alt
            change = {k: source[k] for k in ['variant_id', 'chr', 'pos', 'ref', 'alt', 'motif_name', 'strand', 'candidate_tfs']}
            change.update(status='retained' if ref and alt else 'gained' if alt else 'lost',
                          ref_start=ref['allele_start'] if ref else '', alt_start=alt['allele_start'] if alt else '',
                          reference_start=source['reference_start'], reference_end=source['reference_end'],
                          overlaps_edit=bool((ref and ref['overlaps_edit']) or (alt and alt['overlaps_edit'])))
            for metric in METRICS:
                change['ref_' + metric] = ref[metric] if ref else ''
                change['alt_' + metric] = alt[metric] if alt else ''
                # Zero is a call-level placeholder for absence, not evidence that
                # the true biological motif effect is zero below the call threshold.
                change['delta_' + metric] = (alt[metric] if alt else 0) - (ref[metric] if ref else 0)
            changes.append(change)
    summaries = []
    for variant_id, pair in variants.items():
        row = pair['REF']
        relevant = [r for r in changes if r['variant_id'] == variant_id]
        counts = Counter(r['status'] for r in relevant)
        summaries.append(dict(variant_id=variant_id, chr=row['chr'], pos=int(row['pos']), ref=row['ref'], alt=row['alt'],
                    n_ref_hits=sum(r['variant_id'] == variant_id and r['allele'] == 'REF' for r in calls),
                    n_alt_hits=sum(r['variant_id'] == variant_id and r['allele'] == 'ALT' for r in calls),
                    n_changes=counts['gained']+counts['lost'], n_gained=counts['gained'], n_lost=counts['lost'],
                    n_retained=counts['retained'], n_edit_overlapping_changes=sum(r['overlaps_edit'] and r['status'] != 'retained' for r in relevant)))
    return calls, changes, summaries


def make_report(output, tables, metadata):
    parts = ['<!doctype html><html><head><meta charset="utf-8"><title>Fi-NeMo results</title>',
             '<style>body{font:16px system-ui;max-width:1100px;margin:40px auto;color:#253238}table{border-collapse:collapse}td,th{padding:8px 14px;text-align:left;border-bottom:1px solid #ddd}</style></head><body>',
             '<h1>Fi-NeMo results</h1><p>Motif calls are model predictions. Candidate TF matches can identify several related TFs. Scores are not probabilities or P-values.</p>',
             '<p>REF/ALT gains and losses depend on the call threshold. Uncalled scores are recorded as zero only when calculating call-level differences.</p>',
             '<p>Seqlet recall is not computed by this report. Retain the raw QC table to assess optimizer convergence.</p>',
             '<pre>' + html.escape(json.dumps(metadata, indent=2)) + '</pre>']
    for name, rows, fields in tables:
        parts.append('<h2>' + html.escape(name) + '</h2><table><tr>' + ''.join('<th>' + html.escape(f) + '</th>' for f in fields) + '</tr>')
        for row in rows[:100]:
            parts.append('<tr>' + ''.join('<td>' + html.escape(str(row[f])) + '</td>' for f in fields) + '</tr>')
        parts.append('</table>')
    parts.append('</body></html>')
    (output / 'report.html').write_text('\n'.join(parts))
    archive = shutil.make_archive(str(output.parent / 'finemo_report'), 'zip', output)
    Path(archive).replace(output / 'report.zip')


def main():
    parser = argparse.ArgumentParser()
    for name in ['hits', 'qc', 'regions', 'preparation', 'call-metadata', 'output-dir']:
        parser.add_argument('--' + name, required=True)
    parser.add_argument('--coordinates')
    parser.add_argument('--tf-matches')
    parser.add_argument('--pairing-tolerance', type=int, default=3)
    args = parser.parse_args()
    meta = json.loads(readable(args.preparation).read_text())
    call_meta = json.loads(readable(args.call_metadata).read_text())
    if meta['npz_sha256'] != call_meta['npz_sha256'] or sha256(args.regions) != meta['regions_sha256']:
        raise ValueError('Summary region or call identity differs from preparation')
    matches, hits = tf_matches(args.tf_matches), read_table(args.hits)
    unknown = {hit['motif_name'] for hit in hits} - set(meta['motif_names'])
    if unknown:
        raise ValueError('Hits contain motifs outside the prepared library')
    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    if meta['kind'] == 'variants':
        if not args.coordinates:
            raise ValueError('Variant summary requires allele coordinates')
        with h5py.File(readable(args.coordinates)) as coords:
            if coords.attrs.get('region_sha256') != meta['regions_sha256']:
                raise ValueError('Coordinate map and region identity differ')
            calls, changes, variants = compare_variants(read_rows(args.regions), hits, coords['reference_positions'][:],
                                                        matches, args.pairing_tolerance, meta['crop_start'])
        write_table(output / 'annotated_hits.tsv', calls, CALL_FIELDS)
        write_table(output / 'motif_changes.tsv', changes, CHANGE_FIELDS)
        write_table(output / 'variant_summary.tsv', variants, VARIANT_FIELDS)
    else:
        fields = list(hits[0]) if hits else list(readable(args.hits).read_text().splitlines()[0].split('\t'))
        for hit in hits:
            hit['candidate_tfs'] = ';'.join(matches.get(hit['motif_name'], []))
        write_table(output / 'annotated_hits.tsv', hits, fields + ['candidate_tfs'])
        write_table(output / 'motif_changes.tsv', [], CHANGE_FIELDS)
        write_table(output / 'variant_summary.tsv', [], VARIANT_FIELDS)
    qc = read_table(args.qc)
    counts = Counter(hit['motif_name'] for hit in hits)
    motif_rows = [dict(motif_name=name, n_hits=counts[name], candidate_tfs=';'.join(matches.get(name, []))) for name in meta['motif_names']]
    write_table(output / 'motif_summary.tsv', motif_rows, ['motif_name', 'n_hits', 'candidate_tfs'])
    result = dict(call_meta, pairing_tolerance=args.pairing_tolerance, hit_count=len(hits), qc_rows=len(qc),
                  unconverged_rows=sum(float(row['dual_gap']) > 0.0005 for row in qc),
                  tf_matches_sha256=sha256(args.tf_matches) if args.tf_matches else None,
                  seqlet_recall_computed=False)
    write_json(output / 'metadata.json', result)
    make_report(output, [('Motif counts', motif_rows, ['motif_name', 'n_hits', 'candidate_tfs'])], result)
    print(f'[summary] hits={len(hits)} unconverged_rows={result["unconverged_rows"]}', flush=True)


if __name__ == '__main__':
    main()
