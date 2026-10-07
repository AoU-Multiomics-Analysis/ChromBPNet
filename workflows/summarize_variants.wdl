version 1.0

# The plan contains metadata and array indices, never serialized File inputs.
struct SummaryGroup {
    String model_group
    String cell_type
    Array[Int] indices
}

workflow ChromBPNetSummarizeVariants {
    input {
        File model_manifest
        Array[File] score_files
        String docker_image
        Int summary_memory_gb = 64
        Int summary_disk_gb = 500
        Int summary_max_retries = 2
        Int merge_memory_gb = 64
        Int merge_disk_gb = 500
        Int merge_max_retries = 2
    }

    call PlanModelSummaries {
        input:
            model_manifest = model_manifest,
            n_score_files = length(score_files),
            docker_image = docker_image
    }
    Array[SummaryGroup] groups = read_json(PlanModelSummaries.groups)
    Array[Array[String]] manifest_rows = read_tsv(PlanModelSummaries.rows)

    scatter (group in groups) {
        scatter (index in group.indices) {
            File group_scores = score_files[index]
            String group_model_ids = manifest_rows[index][0]
        }
        # Coerce URI metadata to File before task localization.
        File group_peaks = manifest_rows[group.indices[0]][3]
        call SummarizeModel {
            input:
                score_files = group_scores,
                model_ids = group_model_ids,
                model_group = group.model_group,
                cell_type = group.cell_type,
                peaks = group_peaks,
                docker_image = docker_image,
                memory_gb = summary_memory_gb,
                disk_gb = summary_disk_gb,
                max_retries = summary_max_retries
        }
    }
    call CombineModelEffects {
        input:
            long_files = SummarizeModel.long_effects,
            wide_files = SummarizeModel.wide_effects,
            summary_files = SummarizeModel.fold_summary,
            docker_image = docker_image,
            memory_gb = merge_memory_gb,
            disk_gb = merge_disk_gb,
            max_retries = merge_max_retries
    }
    output {
        File merged_effects = CombineModelEffects.long_effects
        File wide_effects = CombineModelEffects.wide_effects
        File fold_summary = CombineModelEffects.fold_summary
        File merge_log = CombineModelEffects.log
        Array[File] per_model_fold_summaries = SummarizeModel.fold_summary
        Array[File] summary_logs = SummarizeModel.log
        File grouping_log = PlanModelSummaries.log
    }
}

task PlanModelSummaries {
    input {
        File model_manifest
        Int n_score_files
        String docker_image
    }
    command <<<
        set -euo pipefail
        exec > >(tee groups.log) 2>&1
        echo '[groups] Check manifest and plan model summaries'
        python - --manifest '~{sub(model_manifest, "'", "'\"'\"'")}' \
            --n-score-files ~{n_score_files} <<'PY'
        import argparse
        import json
        import re
        import sys
        sys.path.insert(0, '/opt/chrombpnet/scripts')
        from validate_manifest import read_manifest, write_manifest_rows

        parser = argparse.ArgumentParser()
        parser.add_argument('--manifest', required=True)
        parser.add_argument('--n-score-files', type=int, required=True)
        args = parser.parse_args()
        rows = read_manifest(args.manifest)
        if len(rows) != args.n_score_files:
            raise ValueError('Score file count must equal manifest row count')
        groups = {}
        for index, row in enumerate(rows):
            match = re.fullmatch(r'(.+)_fold([0-9]+)', row['model_id'])
            name, fold = (match[1], int(match[2])) if match else (row['model_id'], None)
            group = groups.setdefault(name, dict(cell_type=row['cell_type'], peaks=row['peaks'], members=[]))
            if group['cell_type'] != row['cell_type']:
                raise ValueError(f'Model group has more than one cell type: {name}')
            if group['peaks'] != row['peaks']:
                raise ValueError(f'Model group must use the same peak URI for all folds: {name}')
            if group['members']:
                folds = [member[0] for member in group['members']]
                if fold is None or None in folds:
                    raise ValueError(f'Model group mixes folds and an unsuffixed model: {name}')
                if fold in folds:
                    raise ValueError(f'Duplicate fold index: {name}')
            group['members'].append((fold, index))
        plan = [dict(model_group=name, cell_type=data['cell_type'],
                     indices=[index for fold, index in sorted(data['members'],
                              key=lambda member: -1 if member[0] is None else member[0])])
                for name, data in sorted(groups.items())]
        with open('summary_groups.json', 'w') as stream:
            json.dump(plan, stream)
        write_manifest_rows(rows, 'models.normalized.tsv')
        print(f'[groups] Planned {len(plan)} model groups from {len(rows)} folds', flush=True)
        PY
        echo '[groups] Model summary plan complete'
    >>>
    output {
        File groups = 'summary_groups.json'
        File rows = 'models.normalized.tsv'
        File log = 'groups.log'
    }
    runtime {
        docker: docker_image
        cpu: 1
        memory: '2 GB'
        disks: 'local-disk 10 SSD'
        bootDiskSizeGb: 50
        preemptible: 0
        maxRetries: 2
    }
}

task SummarizeModel {
    input {
        Array[File] score_files
        Array[String] model_ids
        String model_group
        String cell_type
        File peaks
        String docker_image
        Int memory_gb = 64
        Int disk_gb = 500
        Int max_retries = 2
    }
    command <<<
        set -euo pipefail
        exec > >(tee summary.log) 2>&1
        echo '[summary] Start one model group'
        cat > score_files.list <<'CHROMBPNET_SCORE_FILES'
        ~{sep="\n" score_files}
        CHROMBPNET_SCORE_FILES
        # A period cannot occur in a validated model ID, so this cannot collide.
        cat > model_ids.list <<'CHROMBPNET_MODEL_IDS.END'
        ~{sep="\n" model_ids}
        CHROMBPNET_MODEL_IDS.END
        python - --score-files score_files.list --model-ids model_ids.list \
            --model-group '~{sub(model_group, "'", "'\"'\"'")}' \
            --cell-type '~{sub(cell_type, "'", "'\"'\"'")}' \
            --peaks '~{sub(peaks, "'", "'\"'\"'")}' <<'PY'
        import argparse
        import bisect
        import csv
        import gzip
        import re
        import sys
        sys.path.insert(0, '/opt/chrombpnet/scripts')
        from io_utils import check_label, read_scores, readable
        from merge_scores import merge_scores
        from summarize_folds import summarize_folds

        parser = argparse.ArgumentParser()
        parser.add_argument('--score-files', required=True)
        parser.add_argument('--model-ids', required=True)
        parser.add_argument('--model-group', required=True)
        parser.add_argument('--cell-type', required=True)
        parser.add_argument('--peaks', required=True)
        args = parser.parse_args()
        check_label(args.model_group, args.cell_type)
        paths = readable(args.score_files).read_text().splitlines()
        models = readable(args.model_ids).read_text().splitlines()
        if not paths or len(paths) != len(models) or len(models) != len(set(models)):
            raise ValueError('Score files and model IDs must have equal nonzero counts and unique IDs')

        def validate_members():
            for path, model in zip(paths, models):
                check_label(model, args.cell_type)
                match = re.fullmatch(r'(.+)_fold([0-9]+)', model)
                group = match[1] if match else model
                if group != args.model_group:
                    raise ValueError(f'Unexpected model group: {model}')
                fields, rows = read_scores(path)
                if {(row['model_id'], row['cell_type']) for row in rows} != {(model, args.cell_type)}:
                    raise ValueError(f'Score file model or cell type differs from manifest: {model}')
        validate_members()

        # BED/narrowPeak intervals are zero-based, with an excluded end position.
        peak_path = readable(args.peaks)
        opener = gzip.open if peak_path.suffix == '.gz' else open
        intervals = {}
        with opener(peak_path, 'rt') as stream:
            for number, line in enumerate(stream, 1):
                if not line.strip() or line.startswith(('#', 'track ', 'browser ')):
                    continue
                values = line.rstrip('\r\n').split('\t')
                if len(values) < 3 or not values[0] or not values[1].isdigit() or not values[2].isdigit():
                    raise ValueError(f'Invalid peak row: {number}')
                start, end = int(values[1]), int(values[2])
                if end <= start:
                    raise ValueError(f'Invalid peak interval: {number}')
                intervals.setdefault(values[0], []).append((start, end))
        if not intervals:
            raise ValueError('Peak file is empty')
        peak_index = {}
        for chrom, values in intervals.items():
            starts, ends = [], []
            maximum = -1
            for start, end in sorted(values):
                maximum = max(maximum, end)
                starts.append(start)
                ends.append(maximum)
            peak_index[chrom] = (starts, ends)
        del intervals

        print(f'[summary] Average {len(paths)} folds for {args.model_group}', flush=True)
        merge_scores(paths, 'variant_effects.all_models.tsv', 'variant_effects.wide.tsv')
        summarize_folds('variant_effects.all_models.tsv', 'fold_summary.unannotated.tsv')
        count, overlaps = 0, 0
        with open('fold_summary.unannotated.tsv') as source, open('variant_effects.fold_summary.tsv', 'w') as target:
            reader = csv.DictReader(source, delimiter='\t')
            writer = csv.DictWriter(target, fieldnames=reader.fieldnames + ['in_peak'],
                                    delimiter='\t', lineterminator='\n')
            writer.writeheader()
            for row in reader:
                if row['model_group'] != args.model_group or row['cell_type'] != args.cell_type:
                    raise ValueError('Summary contains an unexpected model group or cell type')
                if not row['pos'].isdigit() or int(row['pos']) < 1 or not re.fullmatch(r'[ACGT]+|-', row['allele1']):
                    raise ValueError(f'Invalid variant coordinates or reference allele: {row["variant_id"]}')
                start = int(row['pos']) - 1
                end = start + max(1, len(row['allele1']) if row['allele1'] != '-' else 0)
                starts, ends = peak_index.get(row['chr'], ([], []))
                index = bisect.bisect_left(starts, end) - 1
                overlap = index >= 0 and ends[index] > start
                row['in_peak'] = 'true' if overlap else 'false'
                writer.writerow(row)
                count += 1
                overlaps += overlap
        print(f'[summary] Wrote {count} variants; {overlaps} overlap peaks', flush=True)
        PY
        echo '[summary] Model group complete'
    >>>
    output {
        File long_effects = 'variant_effects.all_models.tsv'
        File wide_effects = 'variant_effects.wide.tsv'
        File fold_summary = 'variant_effects.fold_summary.tsv'
        File log = 'summary.log'
    }
    runtime {
        docker: docker_image
        cpu: 2
        memory: '~{memory_gb}GB'
        disks: 'local-disk ~{disk_gb} SSD'
        bootDiskSizeGb: 50
        preemptible: 0
        maxRetries: max_retries
    }
}

task CombineModelEffects {
    input {
        Array[File] long_files
        Array[File] wide_files
        Array[File] summary_files
        String docker_image
        Int memory_gb = 64
        Int disk_gb = 500
        Int max_retries = 2
    }
    command <<<
        set -euo pipefail
        exec > >(tee combine.log) 2>&1
        # Put SQLite temporary files on the task's data disk too.
        export SQLITE_TMPDIR="$PWD"
        echo '[combine] Combine completed model groups'
        cat > long_files.list <<'CHROMBPNET_LONG_FILES'
        ~{sep="\n" long_files}
        CHROMBPNET_LONG_FILES
        cat > wide_files.list <<'CHROMBPNET_WIDE_FILES'
        ~{sep="\n" wide_files}
        CHROMBPNET_WIDE_FILES
        cat > summary_files.list <<'CHROMBPNET_SUMMARY_FILES'
        ~{sep="\n" summary_files}
        CHROMBPNET_SUMMARY_FILES
        python - --long-files long_files.list --wide-files wide_files.list \
            --summary-files summary_files.list <<'PY'
        import argparse
        import csv
        import json
        import sqlite3
        import sys
        sys.path.insert(0, '/opt/chrombpnet/scripts')
        from io_utils import VARIANT_FIELDS, readable

        parser = argparse.ArgumentParser()
        for name in ['long-files', 'wide-files', 'summary-files']:
            parser.add_argument('--' + name, required=True)
        args = parser.parse_args()
        lists = [readable(path).read_text().splitlines() for path in
                 [args.long_files, args.wide_files, args.summary_files]]
        if not lists[0] or len({len(paths) for paths in lists}) != 1:
            raise ValueError('Model output lists must have equal nonzero counts')

        def reader_for(stream):
            reader = csv.DictReader(stream, delimiter='\t')
            fields = reader.fieldnames or []
            if len(fields) != len(set(fields)) or not set(VARIANT_FIELDS).issubset(fields):
                raise ValueError('Output table has missing or duplicate columns')
            return reader, fields

        def identity(row):
            if None in row or any(value is None for value in row.values()):
                raise ValueError('Output row has a different number of columns')
            return json.dumps([row[field] for field in VARIANT_FIELDS])

        # Keep variant identities and wide rows on disk instead of in memory.
        db = sqlite3.connect('wide_rows.sqlite')
        db.execute('PRAGMA temp_store=FILE')
        db.execute('PRAGMA cache_size=-8192')
        db.execute('CREATE TABLE variants (id TEXT PRIMARY KEY, identity TEXT NOT NULL, ordinal INTEGER UNIQUE)')
        db.execute('CREATE TABLE scores (id TEXT NOT NULL, group_index INTEGER NOT NULL, payload TEXT NOT NULL, PRIMARY KEY (id, group_index))')
        wide_fields = list(VARIANT_FIELDS)
        n_variants = None
        for group_index, path in enumerate(lists[1]):
            with readable(path).open() as stream:
                reader, fields = reader_for(stream)
                metrics = [field for field in fields if field not in VARIANT_FIELDS]
                if not metrics or set(metrics) & set(wide_fields):
                    raise ValueError('Wide outputs have empty or duplicate model columns')
                wide_fields.extend(metrics)
                count = 0
                for ordinal, row in enumerate(reader):
                    key, ident = row['variant_id'], identity(row)
                    if group_index == 0:
                        db.execute('INSERT INTO variants VALUES (?, ?, ?)', (key, ident, ordinal))
                    else:
                        expected = db.execute('SELECT identity FROM variants WHERE id=?', (key,)).fetchone()
                        if expected is None or expected[0] != ident:
                            raise ValueError(f'Variant identity differs between model groups: {key}')
                    db.execute('INSERT INTO scores VALUES (?, ?, ?)',
                               (key, group_index, json.dumps({field: row[field] for field in metrics})))
                    count += 1
            if count == 0 or (n_variants is not None and count != n_variants):
                raise ValueError('Variant set differs between model groups')
            n_variants = count
            db.commit()

        with open('variant_effects.wide.tsv', 'w') as stream:
            writer = csv.DictWriter(stream, fieldnames=wide_fields, delimiter='\t', lineterminator='\n')
            writer.writeheader()
            query = ('SELECT v.id, v.identity, s.group_index, s.payload FROM variants v '
                     'JOIN scores s ON s.id=v.id ORDER BY v.ordinal, s.group_index')
            current, result = None, None
            for key, ident, group_index, payload in db.execute(query):
                if key != current:
                    if result is not None:
                        writer.writerow(result)
                    result = dict(zip(VARIANT_FIELDS, json.loads(ident)))
                    current = key
                result.update(json.loads(payload))
            if result is not None:
                writer.writerow(result)

        # Concatenation preserves all fold scores and all model-specific means.
        # It does not calculate new scores or average across model groups.
        for paths, output, is_summary in [(lists[0], 'variant_effects.all_models.tsv', False),
                                          (lists[2], 'variant_effects.fold_summary.tsv', True)]:
            expected_fields = None
            seen_groups = set()
            with open(output, 'w') as target:
                for group_index, path in enumerate(paths):
                    with readable(path).open() as source:
                        reader, fields = reader_for(source)
                        required = ['model_group', 'cell_type', 'in_peak'] if is_summary else ['model_id', 'cell_type']
                        if not set(required).issubset(fields):
                            raise ValueError('Output table has missing model columns')
                        if expected_fields is None:
                            expected_fields = fields
                            writer = csv.DictWriter(target, fieldnames=fields, delimiter='\t', lineterminator='\n')
                            writer.writeheader()
                        elif fields != expected_fields:
                            raise ValueError('Output column sets differ between model groups')
                        count, group_label = 0, None
                        # A disk table validates duplicates without storing variant IDs in RAM.
                        db.execute('CREATE TEMP TABLE seen (model TEXT, id TEXT, PRIMARY KEY (model, id))')
                        for row in reader:
                            ident = identity(row)
                            expected = db.execute('SELECT identity FROM variants WHERE id=?', (row['variant_id'],)).fetchone()
                            if expected is None or expected[0] != ident:
                                raise ValueError(f'Variant identity differs in combined output: {row["variant_id"]}')
                            label = row['model_group'] if is_summary else row['model_id']
                            db.execute('INSERT INTO seen VALUES (?, ?)', (label, row['variant_id']))
                            if is_summary:
                                current_label = (label, row['cell_type'])
                                if group_label is not None and current_label != group_label:
                                    raise ValueError('Each summary must contain one model group and cell type')
                                group_label = current_label
                                if row['in_peak'] not in ['true', 'false']:
                                    raise ValueError('Invalid in_peak value')
                            writer.writerow(row)
                            count += 1
                        if count == 0 or (is_summary and count != n_variants):
                            raise ValueError('Variant set differs in model summary')
                        if is_summary:
                            if group_label in seen_groups:
                                raise ValueError('Duplicate model group in summaries')
                            seen_groups.add(group_label)
                        for model, model_count in db.execute('SELECT model, COUNT(*) FROM seen GROUP BY model'):
                            if model_count != n_variants:
                                raise ValueError(f'Variant set differs in model output: {model}')
                        db.execute('DROP TABLE seen')
                        db.commit()
        db.close()
        print(f'[combine] Wrote {n_variants} variants across {len(lists[0])} model groups', flush=True)
        PY
        echo '[combine] Model outputs complete'
    >>>
    output {
        File long_effects = 'variant_effects.all_models.tsv'
        File wide_effects = 'variant_effects.wide.tsv'
        File fold_summary = 'variant_effects.fold_summary.tsv'
        File log = 'combine.log'
    }
    runtime {
        docker: docker_image
        cpu: 2
        memory: '~{memory_gb}GB'
        disks: 'local-disk ~{disk_gb} SSD'
        bootDiskSizeGb: 50
        preemptible: 0
        maxRetries: max_retries
    }
}
