"""Average fold scores within model groups, with effect consistency measures."""
import argparse
import csv
import math
import re
import statistics
from io_utils import CORE_METRICS, VARIANT_FIELDS, check_label, readable


def summarize_folds(scores, output):
    print('[summary] Read merged variant scores', flush=True)
    with readable(scores).open() as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        fields = reader.fieldnames or []
        required = VARIANT_FIELDS + ['model_id', 'cell_type'] + CORE_METRICS
        if len(fields) != len(set(fields)) or not set(required).issubset(fields):
            raise ValueError('Score file has missing or duplicate columns')
        rows = list(reader)
    if not rows:
        raise ValueError('Score file is empty')
    # P values are retained in the per-fold output; averaging them is not a test.
    metrics = [field for field in fields if field not in VARIANT_FIELDS + ['model_id', 'cell_type']
               and not field.endswith(('.pval', '_pval', '.pvalue', '_pvalue'))]
    identities = {}
    models = {}
    groups = {}
    for row in rows:
        if None in row or any(value is None for value in row.values()):
            raise ValueError('Score row has a different number of columns')
        model, cell = row['model_id'], row['cell_type']
        check_label(model, cell)
        match = re.fullmatch(r'(.+)_fold([0-9]+)', model)
        group, fold = (match[1], int(match[2])) if match else (model, None)
        if model not in models:
            members = groups.setdefault(group, [])
            if members:
                if cell != models[members[0]]['cell_type']:
                    raise ValueError(f'Model group has more than one cell type: {group}')
                if fold is None or any(models[m]['fold'] is None for m in members):
                    raise ValueError(f'Model group mixes folds with an unsuffixed model: {group}')
                if any(models[m]['fold'] == fold for m in members):
                    raise ValueError(f'duplicate fold index in model group: {group}')
            models[model] = dict(cell_type=cell, fold=fold, rows={})
            members.append(model)
        elif models[model]['cell_type'] != cell:
            raise ValueError(f'Model has more than one cell type: {model}')
        variant = row['variant_id']
        identity = tuple(row[field] for field in VARIANT_FIELDS)
        if variant in identities and identities[variant] != identity:
            raise ValueError(f'Variant coordinates or alleles differ: {variant}')
        identities[variant] = identity
        if variant in models[model]['rows']:
            raise ValueError(f'duplicate model/variant pair: {model}, {variant}')
        for metric in metrics:
            try:
                value = float(row[metric])
            except ValueError as error:
                raise ValueError(f'Nonnumeric score: {metric}') from error
            if not math.isfinite(value):
                raise ValueError(f'Nonfinite score: {metric}')
            row[metric] = value
        models[model]['rows'][variant] = row
    for model, data in models.items():
        if set(data['rows']) != set(identities):
            raise ValueError(f'Variant set differs between models: {model}')
    summary_fields = ['model_group', 'cell_type'] + VARIANT_FIELDS + [
        'n_folds', 'model_ids'] + [metric + '.mean' for metric in metrics] + [
        'logfc.sd', 'percent_change', 'n_positive', 'n_negative', 'n_zero',
        'direction_agreement']
    with open(output, 'w') as stream:
        writer = csv.DictWriter(stream, fieldnames=summary_fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for group in sorted(groups):
            members = sorted(groups[group], key=lambda m: models[m]['fold'] if models[m]['fold'] is not None else -1)
            for variant, identity in identities.items():
                fold_rows = [models[model]['rows'][variant] for model in members]
                result = dict(zip(VARIANT_FIELDS, identity))
                result.update(model_group=group, cell_type=models[members[0]]['cell_type'],
                              n_folds=len(members), model_ids=','.join(members))
                result.update({metric + '.mean': statistics.mean(row[metric] for row in fold_rows)
                               for metric in metrics})
                logfc = [row['logfc'] for row in fold_rows]
                positive = sum(value > 0 for value in logfc)
                negative = sum(value < 0 for value in logfc)
                result.update({'logfc.sd': statistics.stdev(logfc) if len(members) > 1 else '',
                    'percent_change': 100 * (2 ** result['logfc.mean'] - 1),
                    'n_positive': positive, 'n_negative': negative,
                    'n_zero': len(members) - positive - negative,
                    'direction_agreement': max(positive, negative) / len(members)})
                writer.writerow(result)
    print(f'[summary] Wrote {len(identities)} variants across {len(groups)} model groups', flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--scores', required=True, help='Localized merged score TSV')
    parser.add_argument('--output', required=True, help='Fold summary TSV')
    args = parser.parse_args()
    summarize_folds(args.scores, args.output)


if __name__ == '__main__':
    main()
