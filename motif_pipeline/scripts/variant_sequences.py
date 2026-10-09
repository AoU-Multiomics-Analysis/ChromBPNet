"""Make fixed-length allele windows with explicit reference coordinate maps."""
import re

import numpy as np


def edit_bounds(ref, alt):
    prefix = 0
    while prefix < min(len(ref), len(alt)) and ref[prefix] == alt[prefix]:
        prefix += 1
    suffix = 0
    while suffix < min(len(ref), len(alt)) - prefix and ref[-1-suffix] == alt[-1-suffix]:
        suffix += 1
    return prefix, suffix


def allele_windows(row, genome, width):
    ref, alt = row['allele1'], row['allele2']
    if any(not re.fullmatch('[ACGT]+', a) for a in [ref, alt]):
        raise ValueError('Use uppercase ACGT alleles; indels must use VCF-style anchored alleles')
    if ref == alt:
        raise ValueError('REF and ALT must differ')
    if width < 2 or width % 2 or max(len(ref), len(alt)) >= width // 2:
        raise ValueError('Window must be positive and even; alleles must be shorter than half the window')
    chrom, pos = row['chr'], int(row['pos']) - 1
    if chrom not in genome or pos < 0:
        raise ValueError('Variant chromosome or position is invalid')
    start = pos - width // 2
    extra = max(0, len(ref) - len(alt))
    end = start + width + extra
    if start < 0 or end > len(genome[chrom]):
        raise ValueError('Variant does not have a complete model input window')
    sequence = str(genome[chrom][start:end]).upper()
    half = width // 2
    if sequence[half:half + len(ref)] != ref:
        raise ValueError(f"REF mismatch for {row['variant_id']}")
    alt_sequence = (sequence[:half] + alt + sequence[half + len(ref):])[:width]
    prefix, suffix = edit_bounds(ref, alt)
    # Match unchanged prefix/suffix. Pair replacement bases from the left.
    middle_ref = len(ref) - prefix - suffix
    middle_alt = len(alt) - prefix - suffix
    allele_map = list(range(pos, pos + prefix))
    allele_map += [pos + prefix + i if i < middle_ref else -1 for i in range(middle_alt)]
    allele_map += list(range(pos + len(ref) - suffix, pos + len(ref)))
    alt_map = (list(range(start, pos)) + allele_map + list(range(pos + len(ref), end)))[:width]
    maps = np.array([list(range(start, start + width)), alt_map], dtype=np.int64)
    return [sequence[:width], alt_sequence], maps, [(half, half + len(ref)), (half, half + len(alt))]


def one_hot(sequences):
    lookup = np.zeros((256, 4), dtype=np.int8)
    for i, base in enumerate(b'ACGT'):
        lookup[base, i] = 1
    encoded = np.array([list(s.encode('ascii')) for s in sequences], dtype=np.uint8)
    return lookup[encoded].transpose(0, 2, 1)
