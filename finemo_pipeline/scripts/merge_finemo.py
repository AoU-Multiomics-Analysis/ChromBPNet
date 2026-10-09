"""Merge consistent shard tables without losing variants or headers."""
import argparse
from pathlib import Path

from finemo_io import readable, read_table, write_table, write_json


def merge_tables(paths, output, unique_variants=False):
    fields, rows, seen = None, [], set()
    for value in paths:
        path = readable(value)
        current = path.read_text().splitlines()[0].split('\t')
        if fields is not None and fields != current:
            raise ValueError('Shard table headers differ')
        fields = current
        for row in read_table(path):
            if unique_variants and row['variant_id'] in seen:
                raise ValueError('Duplicate variant across shard summaries')
            if unique_variants:
                seen.add(row['variant_id'])
            rows.append(row)
    if fields is None:
        raise ValueError('No shard tables')
    write_table(output, rows, fields)
    return len(rows)


def main():
    parser = argparse.ArgumentParser()
    for name in ['hits-list', 'changes-list', 'variants-list', 'output-dir']:
        parser.add_argument('--' + name, required=True)
    args = parser.parse_args()
    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    counts = {}
    for name, value in [('annotated_hits', args.hits_list), ('motif_changes', args.changes_list), ('variant_summary', args.variants_list)]:
        counts[name] = merge_tables(readable(value).read_text().splitlines(), output / (name + '.tsv'), name == 'variant_summary')
    write_json(output / 'metadata.json', counts)
    print(f'[merge] variants={counts["variant_summary"]}', flush=True)


if __name__ == '__main__':
    main()
