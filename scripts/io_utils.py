"""Input checks shared by scoring and merging."""
import csv
import math
import re
from pathlib import Path

VARIANT_FIELDS = ['chr', 'pos', 'allele1', 'allele2', 'variant_id']
CORE_METRICS = ['logfc', 'jsd', 'active_allele_quantile']


def readable(value):
    value = str(value)
    if '://' in value:
        raise ValueError(f'File localization error: unresolved URI: {value}')
    path = Path(value)
    if not path.is_file():
        raise ValueError(f'File localization error: input is not a file: {value}')
    with path.open('rb') as stream:
        stream.read(1)
    return path


def check_label(model_id, cell_type):
    if not re.fullmatch(r'[A-Za-z0-9_-]+', model_id):
        raise ValueError('model_id must contain only letters, digits, underscores, or hyphens')
    if not cell_type.strip() or any(c in cell_type for c in '\t\r\n'):
        raise ValueError('cell_type must be a nonempty single-line label')


def read_variants(path):
    rows = []
    seen = set()
    with readable(path).open() as stream:
        for line, values in enumerate(csv.reader(stream, delimiter='\t'), 1):
            if values == VARIANT_FIELDS:
                raise ValueError('Variants must be a headerless ChromBPNet TSV')
            if len(values) != 5 or any(not value or any(c in value for c in '\r\n') for value in values):
                raise ValueError(f'Variant row {line} must contain five nonempty fields')
            row = dict(zip(VARIANT_FIELDS, values))
            if not row['pos'].isdigit() or int(row['pos']) < 1:
                raise ValueError(f'Variant row {line} must have a positive 1-based position')
            for field in ['allele1', 'allele2']:
                if not re.fullmatch(r'[ACGT]+|-', row[field]):
                    raise ValueError(f'Invalid {field} at row {line}; use uppercase ACGT or -')
            if row['allele1'] == row['allele2'] == '-':
                raise ValueError('Both alleles cannot be empty')
            if row['variant_id'] in seen:
                raise ValueError(f"duplicate variant_id: {row['variant_id']}")
            seen.add(row['variant_id'])
            rows.append(row)
    if not rows:
        raise ValueError('Variant list is empty')
    return rows


def read_chrom_sizes(path):
    sizes = {}
    with readable(path).open() as stream:
        for row in csv.reader(stream, delimiter='\t'):
            if len(row) != 2 or row[0] in sizes or not row[1].isdigit() or int(row[1]) < 1:
                raise ValueError('Chromosome sizes must have unique names and positive lengths')
            sizes[row[0]] = int(row[1])
    if not sizes:
        raise ValueError('Chromosome sizes are empty')
    return sizes


def validate_reference(rows, genome, sizes, input_length):
    if input_length < 2 or input_length % 2:
        raise ValueError('Model input length must be a positive even number')
    half = input_length // 2
    for row in rows:
        chrom, pos = row['chr'], int(row['pos'])
        ref = '' if row['allele1'] == '-' else row['allele1']
        alt = '' if row['allele2'] == '-' else row['allele2']
        if chrom not in sizes or chrom not in genome or len(genome[chrom]) != sizes[chrom]:
            raise ValueError(f'Chromosome is missing or its FASTA length differs: {chrom}')
        extra = max(0, len(ref) - len(alt))
        if pos - half <= 0 or pos + half + extra > sizes[chrom] or max(len(ref), len(alt)) > half:
            raise ValueError(f"Invalid sequence window: {row['variant_id']}")
        observed = str(genome[chrom][pos - 1:pos - 1 + len(ref)]).upper()
        if observed != ref:
            raise ValueError(f"REF mismatch for {row['variant_id']}: expected {ref}, found {observed}")


def read_scores(path):
    with readable(path).open() as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        fields = reader.fieldnames or []
        required = VARIANT_FIELDS + ['model_id', 'cell_type'] + CORE_METRICS
        if len(fields) != len(set(fields)) or not set(required).issubset(fields):
            raise ValueError('Score file has missing or duplicate columns')
        rows = list(reader)
    if not rows:
        raise ValueError('Score file is empty')
    seen = set()
    for row in rows:
        if None in row or any(value is None for value in row.values()):
            raise ValueError('Score row has a different number of columns')
        check_label(row['model_id'], row['cell_type'])
        if row['variant_id'] in seen:
            raise ValueError('duplicate variant_id in score file')
        seen.add(row['variant_id'])
        for field in CORE_METRICS:
            if not math.isfinite(float(row[field])):
                raise ValueError(f'Nonfinite score: {field}')
    return fields, rows
