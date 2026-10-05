"""Validate manifest metadata without opening model or peak URIs."""
import argparse
import json
from pathlib import Path
from io_utils import check_label, readable


def validate_manifest(rows):
    if not isinstance(rows, list) or not rows:
        raise ValueError('Model manifest must be a nonempty JSON array')
    seen = set()
    for row in rows:
        if not isinstance(row, dict) or set(row) != {'model_id', 'cell_type', 'model', 'peaks'}:
            raise ValueError('Each manifest row needs model_id, cell_type, model, and peaks')
        if any(not isinstance(value, str) or not value.strip() for value in row.values()):
            raise ValueError('Manifest values must be nonempty strings')
        check_label(row['model_id'], row['cell_type'])
        if row['model_id'] in seen:
            raise ValueError(f"duplicate model_id: {row['model_id']}")
        seen.add(row['model_id'])
    return len(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--manifest', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    print('[manifest] Check model metadata', flush=True)
    count = validate_manifest(json.loads(readable(args.manifest).read_text()))
    Path(args.output).write_text(f'{count}\n')
    print(f'[manifest] Validated {count} models', flush=True)


if __name__ == '__main__':
    main()
