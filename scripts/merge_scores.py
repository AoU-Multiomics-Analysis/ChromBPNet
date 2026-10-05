"""Merge per-model scores without averaging distinct cell types or folds."""
import argparse
import csv
from io_utils import VARIANT_FIELDS, read_scores, readable


def merge_scores(paths, long_output, wide_output):
    if not paths:
        raise ValueError('Score file list is empty')
    loaded = [read_scores(path) for path in paths]
    fields, reference = loaded[0]
    identities = {row['variant_id']: tuple(row[f] for f in VARIANT_FIELDS) for row in reference}
    models = []
    indexed = []
    metrics = [field for field in fields if field not in VARIANT_FIELDS + ['model_id', 'cell_type']]
    for current_fields, rows in loaded:
        if current_fields != fields:
            raise ValueError('Score column sets differ between models')
        labels = {(row['model_id'], row['cell_type']) for row in rows}
        if len(labels) != 1:
            raise ValueError('Each score file must contain exactly one model and cell type')
        model, cell_type = next(iter(labels))
        if model in models:
            raise ValueError(f'duplicate model_id in score files: {model}')
        current = {row['variant_id']: tuple(row[f] for f in VARIANT_FIELDS) for row in rows}
        if current != identities:
            raise ValueError('Variant set or variant alleles differ between models')
        models.append(model)
        indexed.append({row['variant_id']: row for row in rows})
    with open(long_output, 'w') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for _, rows in loaded:
            writer.writerows(rows)
    wide_fields = VARIANT_FIELDS + [f'{model}.{field}' for model in models for field in ['cell_type'] + metrics]
    with open(wide_output, 'w') as stream:
        writer = csv.DictWriter(stream, fieldnames=wide_fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for variant in reference:
            row = {field: variant[field] for field in VARIANT_FIELDS}
            for model, scores in zip(models, indexed):
                row.update({f'{model}.{field}': scores[variant['variant_id']][field] for field in ['cell_type'] + metrics})
            writer.writerow(row)
    print(f'[merge] Wrote {len(reference)} variants across {len(models)} models', flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--score-files', required=True, help='Task-local newline-delimited file list')
    parser.add_argument('--long-output', required=True)
    parser.add_argument('--wide-output', required=True)
    args = parser.parse_args()
    paths = readable(args.score_files).read_text().splitlines()
    merge_scores(paths, args.long_output, args.wide_output)


if __name__ == '__main__':
    main()
