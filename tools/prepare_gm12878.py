#!/usr/bin/env python3
"""Download GM12878 ATAC inputs, upload with gsutil, and write a WDL TSV manifest.

Sources: https://www.encodeproject.org/annotations/ENCSR389HIH/
Requires Python 3.9+ and an authenticated gsutil with write access to the bucket.
"""
import argparse
import gzip
import hashlib
import json
import re
import shutil
import subprocess
import tarfile
import time
import urllib.error
import urllib.request
from pathlib import Path, PurePosixPath
from urllib.parse import urlsplit

ENCODE = 'https://www.encodeproject.org'
SOURCES = [
    (ENCODE + '/files/ENCFF142IOR/@@download/ENCFF142IOR.tar.gz',
     'models.tar.gz', '474bf4752d917573552e95e5e6a8612f'),
    (ENCODE + '/files/ENCFF971WEQ/@@download/ENCFF971WEQ.tar.gz',
     'regions.tar.gz', '6b7ccb352e0fc2f76853144ed1cd846d'),
]
UCSC = 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips'
FASTA_MD5 = '1c9dcaddfa41027f17cd8f7a82c7293b'
HDF5_MAGIC = b'\x89HDF\r\n\x1a\n'


def log(message):
    print('[prepare] ' + message, flush=True)


def md5_file(path):
    digest = hashlib.md5()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def download(url, destination, expected_md5):
    """Reuse verified downloads; replace a file only after its checksum passes."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.is_file() and md5_file(destination) == expected_md5:
        log(f'Reuse verified download: {destination}')
        return
    temporary = destination.with_name(destination.name + '.partial')
    log(f'Download {url} to {destination}')
    try:
        for attempt in range(3):
            try:
                request = urllib.request.Request(url, headers={'User-Agent': 'ChromBPNet-input-preparation/1.0'})
                with urllib.request.urlopen(request, timeout=60) as source, temporary.open('wb') as target:
                    shutil.copyfileobj(source, target, length=1024 * 1024)
                break
            except (urllib.error.URLError, TimeoutError, OSError):
                if attempt == 2:
                    raise
                log('Download interrupted; retry from the start')
                time.sleep(2)
        if md5_file(temporary) != expected_md5:
            raise ValueError(f'Download checksum mismatch: {url}')
        temporary.replace(destination)
    finally:
        temporary.unlink(missing_ok=True)


def member_for(archive, basename):
    """Select an exact regular-file basename without extracting archive paths."""
    matches = [member for member in archive.getmembers()
               if PurePosixPath(member.name).name == basename]
    if len(matches) != 1 or not matches[0].isfile():
        raise ValueError(f'Archive must contain exactly one regular file: {basename}')
    member = matches[0]
    path = PurePosixPath(member.name)
    if path.is_absolute() or '..' in path.parts:
        raise ValueError(f'Unsafe archive path: {member.name}')
    return member


def extract_inputs(models, regions, output):
    paths = []
    with tarfile.open(models, 'r:gz') as archive:
        for fold in range(5):
            name = f'model.chrombpnet_nobias.fold_{fold}.ENCSR637XSC.h5'
            member = member_for(archive, name)
            target = output / 'models' / f'fold_{fold}' / 'chrombpnet_nobias.h5'
            target.parent.mkdir(parents=True, exist_ok=True)
            temporary = target.with_suffix('.h5.partial')
            log(f'Extract bias-corrected model: fold {fold}')
            try:
                with archive.extractfile(member) as source, temporary.open('wb') as stream:
                    shutil.copyfileobj(source, stream, length=1024 * 1024)
                with temporary.open('rb') as stream:
                    if stream.read(8) != HDF5_MAGIC:
                        raise ValueError(f'Model is not an HDF5 file: {name}')
                temporary.replace(target)
            finally:
                temporary.unlink(missing_ok=True)
            paths.append(target)
    target = output / 'peaks' / 'GM12878_ATAC.narrowPeak'
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_suffix('.narrowPeak.partial')
    count = 0
    log('Extract and decompress the full training input peak set')
    try:
        with tarfile.open(regions, 'r:gz') as archive:
            member = member_for(archive, 'peaks.all_input_regions.ENCSR637XSC.bed.gz')
            with archive.extractfile(member) as source, gzip.GzipFile(fileobj=source) as peaks, temporary.open('wb') as stream:
                for line in peaks:
                    fields = line.decode('utf-8').rstrip('\r\n').split('\t')
                    if len(fields) != 10 or not fields[0].startswith('chr'):
                        raise ValueError('Peak input must have ten narrowPeak columns and chr-prefixed names')
                    start, end, summit = int(fields[1]), int(fields[2]), int(fields[9])
                    if start < 0 or end <= start or not 0 <= summit < end - start:
                        raise ValueError(f'Invalid peak interval or summit on line {count + 1}')
                    stream.write(line)
                    count += 1
        if count == 0:
            raise ValueError('Peak input is empty')
        temporary.replace(target)
    finally:
        temporary.unlink(missing_ok=True)
    log(f'Checked {count} peaks')
    paths.append(target)
    return paths


def prepare_reference(output):
    """Derive chromosome sizes from the downloaded FASTA to keep them matched."""
    genome = output / 'reference' / 'hg38.fa.gz'
    download(UCSC + '/hg38.fa.gz', genome, FASTA_MD5)
    sizes = output / 'reference' / 'hg38.chrom.sizes'
    records = []
    name, length = None, 0
    log('Read compressed FASTA and create matching chromosome sizes')
    with gzip.open(genome, 'rt', encoding='ascii') as stream:
        for line in stream:
            if line.startswith('>'):
                if name is not None:
                    records.append((name, length))
                name, length = line[1:].split()[0], 0
                if not name.startswith('chr'):
                    raise ValueError(f'FASTA name has no chr prefix: {name}')
            else:
                if name is None:
                    raise ValueError('FASTA sequence has no header')
                length += len(line.strip())
    if name is not None:
        records.append((name, length))
    if not records or any(length == 0 for _, length in records) or len({name for name, _ in records}) != len(records):
        raise ValueError('FASTA contains empty or duplicate chromosome records')
    sizes.write_text(''.join(f'{name}\t{length}\n' for name, length in records))
    return [genome, sizes]


def validate_prefix(prefix):
    parsed = urlsplit(prefix)
    if (any(c in prefix for c in '*?[]{}#\r\n\t') or parsed.scheme != 'gs'
            or not re.fullmatch(r'[a-z0-9][a-z0-9._-]{1,220}[a-z0-9]', parsed.netloc)):
        raise ValueError('Use a gs://bucket/prefix destination without wildcards or control characters')
    return prefix.rstrip('/')


def prepare(output, prefix, skip_reference=False):
    prefix = validate_prefix(prefix)
    output = Path(output).expanduser().resolve()
    if any(c in str(output) for c in '*?[]{}\r\n\t'):
        raise ValueError('The local directory must not contain gsutil wildcards or control characters')
    output.mkdir(parents=True, exist_ok=True)
    for url, name, checksum in SOURCES:
        download(url, output / 'downloads' / name, checksum)
    paths = extract_inputs(output / 'downloads' / 'models.tar.gz', output / 'downloads' / 'regions.tar.gz', output)
    if not skip_reference:
        paths += prepare_reference(output)
    for path in paths:
        destination = prefix + '/' + path.relative_to(output).as_posix()
        log(f'Copy {path} to {destination}')
        subprocess.run(['gsutil', 'cp', str(path), destination], check=True)
    # Write the ready manifest only after every data upload succeeds.
    manifest = output / 'models.tsv'
    temporary = output / 'models.tsv.partial'
    text = 'model_id\tcell_type\tmodel\tpeaks\n'
    for fold in range(5):
        text += (f'GM12878_ATAC_ENCSR637XSC_fold{fold}\tGM12878\t'
                 f'{prefix}/models/fold_{fold}/chrombpnet_nobias.h5\t'
                 f'{prefix}/peaks/GM12878_ATAC.narrowPeak\n')
    temporary.write_text(text, encoding='utf-8')
    temporary.replace(manifest)
    subprocess.run(['gsutil', 'cp', str(manifest), prefix + '/models.tsv'], check=True)
    # This JSON is source provenance, not a script-argument wrapper.
    provenance = output / 'sources.json'
    provenance.write_text(json.dumps({'annotation': 'ENCSR389HIH', 'experiment': 'ENCSR637XSC',
        'assembly': 'GRCh38', 'downloads': [{'url': url, 'md5': checksum} for url, _, checksum in SOURCES],
        'reference_url': None if skip_reference else UCSC + '/hg38.fa.gz',
        'reference_md5': None if skip_reference else FASTA_MD5}, indent=2) + '\n')
    log(f'Complete. Local manifest: {manifest}')
    log(f'WDL model_manifest: {prefix}/models.tsv')
    if not skip_reference:
        log(f'WDL genome: {prefix}/reference/hg38.fa.gz')
        log(f'WDL chrom_sizes: {prefix}/reference/hg38.chrom.sizes')
        log('Leave genome_index unset; the WDL creates it after decompression')


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--output-dir', required=True, type=Path, help='Local download and extracted-file directory')
    parser.add_argument('--gcs-prefix', required=True, help='Destination gs://bucket/prefix; files are copied with gsutil')
    parser.add_argument('--skip-reference', action='store_true', help='Download only models and peaks; supply your own hg38 reference in the WDL')
    args = parser.parse_args()
    validate_prefix(args.gcs_prefix)
    if shutil.which('gsutil') is None:
        parser.error('gsutil is not on PATH; install and authenticate the Google Cloud SDK first')
    prepare(args.output_dir, args.gcs_prefix, args.skip_reference)


if __name__ == '__main__':
    try:
        main()
    except (ValueError, OSError, EOFError, tarfile.TarError, subprocess.CalledProcessError) as error:
        raise SystemExit(f'[prepare] ERROR: {error}')
