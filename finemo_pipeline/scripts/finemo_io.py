"""Local input checks and simple TSV/provenance helpers."""
import csv
import hashlib
import json
from pathlib import Path

ROW_FIELDS = ['region_index', 'variant_id', 'allele', 'chr', 'pos', 'ref', 'alt',
              'edit_start', 'edit_end', 'window_start', 'insertion_anchor']


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


def sha256(path):
    digest = hashlib.sha256()
    with readable(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def write_json(path, data):
    Path(path).write_text(json.dumps(data, indent=2, allow_nan=False) + '\n')


def read_rows(path):
    with readable(path).open() as stream:
        values = list(csv.reader(stream, delimiter='\t'))
    if not values or any(len(row) != len(ROW_FIELDS) for row in values):
        raise ValueError('Allele rows require the exact headerless preparation schema')
    rows = [dict(zip(ROW_FIELDS, row)) for row in values]
    if [int(row['region_index']) for row in rows] != list(range(len(rows))):
        raise ValueError('Allele rows must retain exact consecutive region indices')
    return rows


def read_table(path):
    with readable(path).open() as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        if not reader.fieldnames or len(set(reader.fieldnames)) != len(reader.fieldnames):
            raise ValueError('TSV requires unique column names')
        rows = list(reader)
    if any(None in row or any(v is None for v in row.values()) for row in rows):
        raise ValueError('TSV row width differs from header')
    return rows


def write_table(path, rows, fields):
    with Path(path).open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def tf_matches(path):
    if path is None:
        return {}
    rows = read_table(path)
    result = {}
    for row in rows:
        key = row.get('pattern', row.get('motif_id', row.get('Query_ID')))
        tf = row.get('candidate_tf', row.get('target_name', row.get('target_label', row.get('Target_ID'))))
        if not key or not tf:
            raise ValueError('TF matches require pattern/candidate_tf (or Tomtom Query_ID/Target_ID)')
        result.setdefault(key, set()).add(tf)
    return {key: sorted(values) for key, values in result.items()}
