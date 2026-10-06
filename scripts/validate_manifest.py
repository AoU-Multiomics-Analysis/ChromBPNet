"""Validate manifest metadata without opening model or peak URIs."""
import argparse
import csv
from pathlib import Path
from io_utils import check_label, readable

FIELDS = ['model_id', 'cell_type', 'model', 'peaks']


def validate_manifest(rows):
    if not isinstance(rows, list) or not rows:
        raise ValueError('Model manifest must contain at least one TSV data row')
    seen = set()
    for row in rows:
        if not isinstance(row, dict) or set(row) != {'model_id', 'cell_type', 'model', 'peaks'}:
            raise ValueError('Each manifest row needs model_id, cell_type, model, and peaks')
        if any(not isinstance(value, str) or not value.strip() or any(c in value for c in '\t\r\n') for value in row.values()):
            raise ValueError('Manifest values must be nonempty single-line strings')
        check_label(row['model_id'], row['cell_type'])
        if row['model_id'] in seen:
            raise ValueError(f"duplicate model_id: {row['model_id']}")
        seen.add(row['model_id'])
    return len(rows)


def read_manifest(path):
    with readable(path).open(encoding='utf-8-sig', newline='') as stream:
        reader = csv.DictReader(stream, delimiter='\t', quoting=csv.QUOTE_NONE)
        if reader.fieldnames != FIELDS:
            raise ValueError('TSV header must be model_id, cell_type, model, peaks in that order')
        rows = list(reader)
    validate_manifest(rows)
    return rows


def write_manifest_rows(rows, path):
    """Write plain headerless TSV metadata for WDL read_tsv, without CSV quoting."""
    validate_manifest(rows)
    with open(path, 'w') as stream:
        for row in rows:
            stream.write('\t'.join(row[field] for field in FIELDS) + '\n')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--manifest', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--rows-output', required=True)
    args = parser.parse_args()
    print('[manifest] Check model metadata', flush=True)
    rows = read_manifest(args.manifest)
    write_manifest_rows(rows, args.rows_output)
    count = len(rows)
    Path(args.output).write_text(f'{count}\n')
    print(f'[manifest] Validated {count} models', flush=True)


if __name__ == '__main__':
    main()
