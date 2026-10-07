#!/usr/bin/env python3
"""Download all released TenK10K folds and peaks; upload a Terra model manifest.

Requires Python 3.9+ and authenticated gsutil. No analysis job is submitted.
The adjacent source index pins Hugging Face paths and SHA-256 checksums.
"""
import argparse
import csv
import hashlib
import io
import itertools
import json
import math
import re
import shutil
import stat
import subprocess
import time
import urllib.error
import urllib.request
import zipfile
from pathlib import Path, PurePosixPath
from urllib.parse import quote

from prepare_gm12878 import HDF5_MAGIC, prepare_reference, validate_prefix

SOURCE_INDEX = Path(__file__).with_name('tenk10k_sources.json')
PEAK_ARCHIVE = 'Miscellaneous/TenK10K_ATAC_ct_specific_peaks.zip'
COMBINED_PEAKS = 'Miscellaneous/TenK10K_ATAC_MACS3_Combined_Peaks_Annotated.bed'
CSV_COLUMNS = ['chrom', 'start', 'end', 'name', 'score', 'strand',
               'signal_value', 'p_value', 'q_value', 'peak']
# These two released .csv files contain headerless BED6 rather than CSV.
BED_CELL_TYPES = {'CD14_Mono', 'CD4_TCM'}
FIELDS = ['model_id', 'cell_type', 'model', 'peaks']


def log(message):
    print('[tenk10k] ' + message, flush=True)


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def download(url, destination, expected_sha256, expected_size):
    """Reuse verified bytes. Never publish a partial or unverified download."""
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    if (destination.is_file() and destination.stat().st_size == expected_size
            and sha256_file(destination) == expected_sha256):
        log(f'Reuse verified download: {destination}')
        return
    temporary = destination.with_name(destination.name + '.partial')
    log(f'Download {url}')
    try:
        for attempt in range(3):
            try:
                request = urllib.request.Request(url, headers={'User-Agent': 'ChromBPNet-TenK10K/1.0'})
                with urllib.request.urlopen(request, timeout=60) as source, temporary.open('wb') as target:
                    shutil.copyfileobj(source, target, length=1024 * 1024)
                break
            except (urllib.error.URLError, TimeoutError, OSError):
                if attempt == 2:
                    raise
                log('Download interrupted; retry from the start')
                time.sleep(2)
        if temporary.stat().st_size != expected_size or sha256_file(temporary) != expected_sha256:
            raise ValueError(f'Download size or checksum mismatch: {url}')
        temporary.replace(destination)
    finally:
        temporary.unlink(missing_ok=True)


def select_assets(catalog):
    """Require one bias-corrected model for each of five folds per cell type."""
    if catalog.get('repository') != 'anglixue/TenK10K_multiome':
        raise ValueError('Source index must name the official TenK10K repository')
    if not re.fullmatch(r'[0-9a-f]{40}', catalog.get('revision', '')):
        raise ValueError('Source revision must be an immutable 40-character commit SHA')
    cells = catalog.get('peak_cell_types', [])
    if (not cells or len(cells) != len(set(cells))
            or any(not re.fullmatch(r'[A-Za-z0-9_-]+', cell) for cell in cells)):
        raise ValueError('Source index needs unique, safe peak cell types')
    sources, models, folds = {}, [], {}
    pattern = re.compile(r'ChromBPNet/([A-Za-z0-9_-]+)/([A-Za-z0-9_-]+_fold([0-4]))/models/chrombpnet_nobias\.h5')
    for item in catalog.get('files', []):
        path = item['path']
        if path in sources:
            raise ValueError(f'duplicate source path: {path}')
        if (PurePosixPath(path).is_absolute() or '..' in PurePosixPath(path).parts
                or not isinstance(item.get('size'), int) or item['size'] <= 0
                or not re.fullmatch(r'[0-9a-f]{64}', item.get('sha256', ''))):
            raise ValueError(f'Invalid source path, size, or SHA-256: {path}')
        sources[path] = item
        if path in {PEAK_ARCHIVE, COMBINED_PEAKS}:
            continue
        match = pattern.fullmatch(path)
        if match is None:
            raise ValueError(f'Invalid model path or fold: {path}')
        cell, _, fold = match.groups()
        fold = int(fold)
        if cell not in cells:
            raise ValueError(f'Model cell type has no released peaks: {cell}')
        if fold in folds.setdefault(cell, set()):
            raise ValueError(f'duplicate model fold: {cell}, {fold}')
        folds[cell].add(fold)
        models.append({'cell_type': cell, 'fold': fold, 'source': item,
                       'relative': f'models/{cell}/fold_{fold}/chrombpnet_nobias.h5'})
    if not models:
        raise ValueError('Source index has no model folds')
    for cell, found in folds.items():
        if found != set(range(5)):
            raise ValueError(f'Incomplete model folds for {cell}: {sorted(found)}')
    expected_models = catalog.get('model_cell_types', [])
    if (not expected_models or len(expected_models) != len(set(expected_models))
            or not set(expected_models).issubset(cells)):
        raise ValueError('Source index needs unique model cell types with matching peak cells')
    if set(folds) != set(expected_models):
        raise ValueError(f'Model cell types missing or unexpected: {sorted(set(folds) ^ set(expected_models))}')
    if PEAK_ARCHIVE not in sources or COMBINED_PEAKS not in sources:
        raise ValueError('Source index is missing the peak archive or combined peaks')
    return sorted(models, key=lambda item: (item['cell_type'], item['fold'])), sources


def asset_url(catalog, path):
    return ('https://huggingface.co/datasets/' + catalog['repository'] + '/resolve/'
            + catalog['revision'] + '/' + quote(path, safe='/'))


def peak_relative(cell):
    return 'peaks/' + cell + ('.bed' if cell in BED_CELL_TYPES else '.narrowPeak')


def validate_peak(fields, line_number):
    if (not 3 <= len(fields) <= 10 or not re.fullmatch(r'chr[A-Za-z0-9_.-]+', fields[0])
            or any(not field or any(ord(c) < 32 or ord(c) == 127 for c in field) for field in fields)):
        raise ValueError(f'Invalid peak columns or chromosome on row {line_number}')
    start, end = int(fields[1]), int(fields[2])
    if not 0 <= start < end:
        raise ValueError(f'Invalid peak interval on row {line_number}')
    if len(fields) == 10:
        summit = int(fields[9])
        if not 0 <= summit < end - start:
            raise ValueError(f'Invalid peak summit on row {line_number}')
        if any(not math.isfinite(float(value)) for value in [fields[4]] + fields[6:9]):
            raise ValueError(f'Nonfinite narrowPeak value on row {line_number}')


def extract_peaks(archive_path, output, expected_cells):
    """Read exact ZIP members. Convert CSV; retain released BED6 and summits."""
    output = Path(output)
    report = {}
    with zipfile.ZipFile(archive_path) as archive:
        members = {}
        for member in archive.infolist():
            path = PurePosixPath(member.filename)
            mode = member.external_attr >> 16
            if path.is_absolute() or '..' in path.parts or stat.S_ISLNK(mode):
                raise ValueError(f'Unsafe archive member: {member.filename}')
            if member.is_dir() or path.parts[0] == '__MACOSX':
                continue
            if len(path.parts) != 2 or path.parts[0] != 'TenK10K_ATAC_ct_specific_peaks' or path.suffix != '.csv':
                raise ValueError(f'Unexpected peak archive member: {member.filename}')
            cell = path.stem
            if cell in members:
                raise ValueError(f'duplicate peak archive member: {cell}')
            members[cell] = member
        if set(members) != set(expected_cells):
            raise ValueError(f'Peak cells missing or unexpected: {sorted(set(members) ^ set(expected_cells))}')
        for cell in sorted(members):
            target = output / peak_relative(cell)
            target.parent.mkdir(parents=True, exist_ok=True)
            temporary = target.with_name(target.name + '.partial')
            count = 0
            try:
                with archive.open(members[cell]) as binary, io.TextIOWrapper(binary, encoding='utf-8-sig', newline='') as stream, \
                        temporary.open('w', encoding='utf-8', newline='') as result:
                    first = stream.readline()
                    is_bed = first.startswith('chr') and '\t' in first
                    if is_bed != (cell in BED_CELL_TYPES):
                        raise ValueError(f'Peak format differs from the pinned release: {cell}')
                    if is_bed:
                        rows = csv.reader(itertools.chain([first], stream), delimiter='\t', quoting=csv.QUOTE_NONE)
                    else:
                        header = next(csv.reader([first]))
                        has_index = header[:1] in [[''], ['Unnamed: 0']]
                        if (header[1:] if has_index else header) != CSV_COLUMNS:
                            raise ValueError(f'Unexpected peak CSV header: {cell}')
                        rows = csv.reader(stream)
                    for row in rows:
                        if not is_bed:
                            if len(row) != 10 + int(has_index):
                                raise ValueError(f'Invalid CSV column count for {cell}')
                            row = row[1:] if has_index else row
                        validate_peak(row, count + 1)
                        result.write('\t'.join(row) + '\n')
                        count += 1
                if count == 0:
                    raise ValueError(f'Peak file is empty: {cell}')
                temporary.replace(target)
            finally:
                temporary.unlink(missing_ok=True)
            report[cell] = {'relative': peak_relative(cell), 'source_member': members[cell].filename,
                            'peak_count': count, 'centering': 'midpoint' if is_bed else 'summit',
                            'sha256': sha256_file(target)}
            log(f'Checked {count} peaks: {cell}')
    return report


def write_json(path, value):
    temporary = path.with_name(path.name + '.partial')
    temporary.write_text(json.dumps(value, indent=2) + '\n', encoding='utf-8')
    temporary.replace(path)


def write_manifest(path, models, prefix):
    temporary = path.with_name(path.name + '.partial')
    with temporary.open('w', encoding='utf-8', newline='') as stream:
        stream.write('\t'.join(FIELDS) + '\n')
        for item in models:
            cell, fold = item['cell_type'], item['fold']
            stream.write('\t'.join([f'TenK10K_ATAC_{cell}_fold{fold}', cell,
                                     prefix + '/' + item['relative'], prefix + '/' + peak_relative(cell)]) + '\n')
    temporary.replace(path)


def prepare(output, prefix, skip_reference=False, source_index=SOURCE_INDEX, plan_only=False):
    if any(ord(c) < 32 or ord(c) == 127 for c in prefix):
        raise ValueError('GCS destination contains control characters')
    prefix = validate_prefix(prefix)
    output = Path(output).expanduser().resolve()
    t7 = Path('/Volumes/T7')
    if t7 in output.parents and not t7.is_mount():
        raise ValueError('T7 must be mounted before writing its output directory')
    if any(c in str(output) for c in '*?[]{}') or any(ord(c) < 32 or ord(c) == 127 for c in str(output)):
        raise ValueError('Local output path contains wildcards or control characters')
    catalog = json.loads(Path(source_index).read_text(encoding='utf-8'))
    models, sources = select_assets(catalog)
    output.mkdir(parents=True, exist_ok=True)
    model_cells = {item['cell_type'] for item in models}
    unmatched = sorted(set(catalog['peak_cell_types']) - model_cells)
    log(f'{len(models)} models, {len(model_cells)} model cell types, {len(catalog["peak_cell_types"])} peak cell types')
    log(f'Hugging Face revision: {catalog["revision"]}')
    log(f'Models and source peaks: {sum(item["size"] for item in sources.values()) / 1e9:.2f} GB to download')
    if unmatched:
        log('Peak cell types without released models: ' + ', '.join(unmatched))
    if plan_only:
        write_manifest(output / 'models.planned.tsv', models, prefix)
        write_json(output / 'download_plan.json', {**catalog, 'gcs_prefix': prefix,
                   'peak_cell_types_without_models': unmatched, 'include_reference': not skip_reference})
        log('Plan only. Files have not been downloaded or uploaded; models.planned.tsv is not ready for Terra.')
        return
    paths = []
    for item in models:
        target = output / item['relative']
        source = item['source']
        download(asset_url(catalog, source['path']), target, source['sha256'], source['size'])
        with target.open('rb') as stream:
            if stream.read(8) != HDF5_MAGIC:
                raise ValueError(f'Model is not HDF5: {source["path"]}')
        paths.append(target)
    archive = output / 'downloads' / PurePosixPath(PEAK_ARCHIVE).name
    combined = output / 'peaks' / PurePosixPath(COMBINED_PEAKS).name
    for source_name, target in [(PEAK_ARCHIVE, archive), (COMBINED_PEAKS, combined)]:
        source = sources[source_name]
        download(asset_url(catalog, source_name), target, source['sha256'], source['size'])
        paths.append(target)
    peaks = extract_peaks(archive, output, catalog['peak_cell_types'])
    paths.extend(output / peaks[cell]['relative'] for cell in sorted(peaks))
    if not skip_reference:
        paths.extend(prepare_reference(output))
    provenance = output / 'sources.json'
    write_json(provenance, {**catalog, 'assembly': 'GRCh38', 'gcs_prefix': prefix,
               'peak_cell_types_without_models': unmatched, 'normalized_peaks': peaks,
               'reference': None if skip_reference else {
                   'source': 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz',
                   'sha256': sha256_file(output / 'reference/hg38.fa.gz')},
               'files_uploaded': [{'relative': path.relative_to(output).as_posix(),
                                   'sha256': sha256_file(path)} for path in paths]})
    paths.append(provenance)
    for path in paths:
        destination = prefix + '/' + path.relative_to(output).as_posix()
        log(f'Copy {path} to {destination}')
        subprocess.run(['gsutil', 'cp', str(path), destination], check=True)
    # The four-column manifest is published only after all data uploads succeed.
    manifest = output / 'models.tsv'
    write_manifest(manifest, models, prefix)
    subprocess.run(['gsutil', 'cp', str(manifest), prefix + '/models.tsv'], check=True)
    log(f'Complete. WDL model_manifest: {prefix}/models.tsv')
    if not skip_reference:
        log(f'WDL genome: {prefix}/reference/hg38.fa.gz')
        log(f'WDL chrom_sizes: {prefix}/reference/hg38.chrom.sizes')
        log('Leave genome_index unset; the WDL creates it after decompression')


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--output-dir', required=True, type=Path)
    parser.add_argument('--gcs-prefix', required=True, help='Destination gs://bucket/prefix')
    parser.add_argument('--skip-reference', action='store_true', help='Use your existing hg38 reference in Terra')
    parser.add_argument('--plan-only', action='store_true', help='Write a download plan without network access or uploads')
    parser.add_argument('--source-index', type=Path, default=SOURCE_INDEX, help='Versioned source inventory; defaults to the bundled release')
    args = parser.parse_args()
    if not args.plan_only and shutil.which('gsutil') is None:
        parser.error('Install and authenticate gsutil before downloading and uploading')
    prepare(args.output_dir, args.gcs_prefix, args.skip_reference, args.source_index, args.plan_only)


if __name__ == '__main__':
    try:
        main()
    except (ValueError, KeyError, OSError, EOFError, zipfile.BadZipFile, subprocess.CalledProcessError) as error:
        raise SystemExit(f'[tenk10k] ERROR: {error}')
