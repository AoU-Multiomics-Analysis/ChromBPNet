version 1.0

import "discover_motifs.wdl" as motifs

# This plan contains labels and row indices. Input paths remain TSV metadata
# until the workflow declares each model and peak path as File.
struct MotifGroup {
    String model_group
    String cell_type
    Array[Int] indices
}

workflow ChromBPNetMotifDiscoveryManifest {
    input {
        File model_manifest
        File genome
        File? genome_index
        File chrom_sizes
        File motif_database
        String docker_image
        String head = "counts"
        Int input_length = 2114
        Int random_seed = 1234
        Int num_backgrounds = 20
        Int batch_size = 64
        Int? max_peaks
        Int discovery_window = 400
        Int max_seqlets = 1000000
        Int n_leiden = 2
        Float match_qvalue = 0.05
        Int n_matches = 5
        Int contribution_memory_gb = 64
        Int contribution_disk_gb = 100
        Int averaging_disk_gb = 150
        Int motif_cpu = 16
        Int motif_memory_gb = 128
        Int motif_disk_gb = 100
        Int num_preempt = 0
    }

    call PlanMotifGroups {
        input:
            model_manifest = model_manifest,
            head = head,
            docker_image = docker_image
    }
    Array[MotifGroup] groups = read_json(PlanMotifGroups.groups)
    Array[Array[String]] manifest_rows = read_tsv(PlanMotifGroups.rows)

    scatter (group in groups) {
        String model_group = group.model_group
        Array[Int] group_indices = group.indices
        Int peak_index = group_indices[0]
        scatter (index in group_indices) {
            File group_models = manifest_rows[index][2]
        }
        File group_peaks = manifest_rows[peak_index][3]

        call motifs.ChromBPNetMotifDiscovery as ModelMotifs {
            input:
                fold_models = group_models,
                peaks = group_peaks,
                cell_type = group.cell_type,
                genome = genome,
                genome_index = genome_index,
                chrom_sizes = chrom_sizes,
                motif_database = motif_database,
                docker_image = docker_image,
                head = head,
                input_length = input_length,
                random_seed = random_seed,
                num_backgrounds = num_backgrounds,
                batch_size = batch_size,
                max_peaks = max_peaks,
                discovery_window = discovery_window,
                max_seqlets = max_seqlets,
                n_leiden = n_leiden,
                match_qvalue = match_qvalue,
                n_matches = n_matches,
                contribution_memory_gb = contribution_memory_gb,
                contribution_disk_gb = contribution_disk_gb,
                averaging_disk_gb = averaging_disk_gb,
                motif_cpu = motif_cpu,
                motif_memory_gb = motif_memory_gb,
                motif_disk_gb = motif_disk_gb,
                num_preempt = num_preempt
        }

        # These paths are output-manifest metadata only. String coercion retains
        # Cromwell's output URIs and avoids localizing large result files again.
        String modisco_uri = ModelMotifs.modisco_motifs
        String contributions_uri = ModelMotifs.averaged_contributions
        String regions_uri = ModelMotifs.interpreted_regions
        String report_uri = ModelMotifs.report_bundle
        String meme_uri = ModelMotifs.discovered_motifs_meme
        String matches_uri = ModelMotifs.candidate_tf_matches
        String inventory_uri = ModelMotifs.motif_inventory
        String tomtom_uri = ModelMotifs.tomtom_results
        String result_row = group.model_group + "\t" + group.cell_type + "\t" + head + "\t" + modisco_uri + "\t" + contributions_uri + "\t" + regions_uri + "\t" + report_uri + "\t" + meme_uri + "\t" + matches_uri + "\t" + inventory_uri + "\t" + tomtom_uri
        String result_row_arg = "--result-row='" + sub(result_row, "'", "'\"'\"'") + "'"
    }

    call WriteMotifManifests {
        input:
            model_manifest = model_manifest,
            group_plan = PlanMotifGroups.groups,
            head = head,
            result_row_args = result_row_arg,
            docker_image = docker_image
    }

    output {
        File model_manifest_with_modisco = WriteMotifManifests.models
        File modisco_manifest = WriteMotifManifests.results
        Array[String] model_groups = model_group
        Array[File] modisco_motifs = ModelMotifs.modisco_motifs
        Array[File] averaged_contributions = ModelMotifs.averaged_contributions
        Array[File] interpreted_regions = ModelMotifs.interpreted_regions
        Array[File] report_bundle = ModelMotifs.report_bundle
        Array[File] discovered_motifs_meme = ModelMotifs.discovered_motifs_meme
        Array[File] candidate_tf_matches = ModelMotifs.candidate_tf_matches
        Array[File] motif_inventory = ModelMotifs.motif_inventory
        Array[File] tomtom_results = ModelMotifs.tomtom_results
        Array[Array[File]] per_fold_contributions = ModelMotifs.per_fold_contributions
        Array[Array[File]] per_fold_metadata = ModelMotifs.per_fold_metadata
        Array[Array[File]] per_fold_logs = ModelMotifs.per_fold_logs
        Array[File] preparation_metadata = ModelMotifs.preparation_metadata
        Array[File] averaging_metadata = ModelMotifs.averaging_metadata
        Array[File] discovery_metadata = ModelMotifs.discovery_metadata
        Array[File] report_metadata = ModelMotifs.report_metadata
        Array[File] preparation_logs = ModelMotifs.preparation_log
        Array[File] averaging_logs = ModelMotifs.averaging_log
        Array[File] discovery_logs = ModelMotifs.discovery_log
        Array[File] report_logs = ModelMotifs.report_log
        File grouping_log = PlanMotifGroups.log
        File manifest_log = WriteMotifManifests.log
    }
}

task PlanMotifGroups {
    input {
        File model_manifest
        String head
        String docker_image
    }
    command <<<
        set -euo pipefail
        exec > >(tee motif_groups.log) 2>&1
        echo '[manifest] Check model manifest and plan five-fold model groups'
        python - --manifest '~{sub(model_manifest, "'", "'\"'\"'")}' \
            --head '~{sub(head, "'", "'\"'\"'")}' <<'PY'
        import argparse
        import csv
        import json
        import re
        from pathlib import Path

        parser = argparse.ArgumentParser()
        parser.add_argument('--manifest', required=True)
        parser.add_argument('--head', choices=['counts', 'profile'], required=True)
        args = parser.parse_args()
        if '://' in args.manifest or any(c in args.manifest for c in '\r\n'):
            raise ValueError('File localization error: unresolved manifest URI or invalid local path')
        fields = ['model_id', 'cell_type', 'model', 'peaks']
        with Path(args.manifest).open(encoding='utf-8-sig', newline='') as stream:
            reader = csv.DictReader(stream, delimiter='\t', quoting=csv.QUOTE_NONE)
            if reader.fieldnames != fields:
                raise ValueError('TSV header must be model_id, cell_type, model, peaks in that order')
            rows = list(reader)
        if not rows:
            raise ValueError('Model manifest must contain at least one five-fold model group')
        seen, grouped = set(), {}
        for index, row in enumerate(rows):
            if set(row) != set(fields) or any(not isinstance(value, str) or not value.strip()
                    or any(c in value for c in '\t\r\n') for value in row.values()):
                raise ValueError(f'Invalid manifest values at data row {index + 1}')
            match = re.fullmatch(r'([A-Za-z0-9_-]+)_fold([0-4])', row['model_id'])
            if not match or row['model_id'] in seen:
                raise ValueError(f'Model IDs must be unique and end in _fold0 through _fold4: {row["model_id"]}')
            seen.add(row['model_id'])
            name, fold = match[1], int(match[2])
            group = grouped.setdefault(name, dict(cell_type=row['cell_type'], peaks=row['peaks'], members={}))
            if group['cell_type'] != row['cell_type'] or group['peaks'] != row['peaks']:
                raise ValueError(f'All folds must share a cell type and peak URI: {name}')
            group['members'][fold] = index
        groups = []
        for name, group in sorted(grouped.items()):
            if set(group['members']) != set(range(5)):
                raise ValueError(f'Model group needs all five folds 0 through 4: {name}')
            indices = [group['members'][fold] for fold in range(5)]
            if len({rows[index]['model'] for index in indices}) != 5:
                raise ValueError(f'Model URIs must be distinct within a group: {name}')
            groups.append(dict(model_group=name, cell_type=group['cell_type'], indices=indices))
        # The TSV retains input URI metadata; WDL declares File values before use.
        with open('models.normalized.tsv', 'w') as stream:
            for row in rows:
                stream.write('\t'.join(row[field] for field in fields) + '\n')
        Path('motif_groups.json').write_text(json.dumps(groups) + '\n')
        print(f'[manifest] Planned {len(groups)} model groups from {len(rows)} folds', flush=True)
        PY
        echo '[manifest] Group plan complete'
    >>>
    output {
        File groups = 'motif_groups.json'
        File rows = 'models.normalized.tsv'
        File log = 'motif_groups.log'
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

task WriteMotifManifests {
    input {
        File model_manifest
        File group_plan
        String head
        # Each argument is safely quoted by the workflow from output metadata.
        Array[String] result_row_args
        String docker_image
    }
    command <<<
        set -euo pipefail
        exec > >(tee motif_manifests.log) 2>&1
        echo '[manifest] Write model and MoDISco result manifests'
        python - --manifest '~{sub(model_manifest, "'", "'\"'\"'")}' \
            --group-plan '~{sub(group_plan, "'", "'\"'\"'")}' \
            --head '~{sub(head, "'", "'\"'\"'")}' \
            ~{sep=" " result_row_args} <<'PY'
        import argparse
        import csv
        import json
        from pathlib import Path

        parser = argparse.ArgumentParser()
        parser.add_argument('--manifest', required=True)
        parser.add_argument('--group-plan', required=True)
        parser.add_argument('--head', choices=['counts', 'profile'], required=True)
        parser.add_argument('--result-row', action='append', required=True)
        args = parser.parse_args()
        input_fields = ['model_id', 'cell_type', 'model', 'peaks']
        file_fields = ['modisco_motifs', 'averaged_contributions', 'interpreted_regions',
                       'report_bundle', 'discovered_motifs_meme', 'candidate_tf_matches',
                       'motif_inventory', 'tomtom_results']
        result_fields = ['model_group', 'cell_type', 'head'] + file_fields

        def local(value):
            if '://' in value or any(c in value for c in '\r\n'):
                raise ValueError('File localization error: unresolved manifest URI or invalid local path')
            return Path(value)

        def valid_values(row):
            return all(isinstance(value, str) and value.strip() and not any(c in value for c in '\t\r\n')
                       for value in row.values())

        with local(args.manifest).open(encoding='utf-8-sig', newline='') as stream:
            reader = csv.DictReader(stream, delimiter='\t', quoting=csv.QUOTE_NONE)
            if reader.fieldnames != input_fields:
                raise ValueError('Input model manifest header differs from the four-column schema')
            rows = list(reader)
        if not rows or any(set(row) != set(input_fields) or not valid_values(row) for row in rows):
            raise ValueError('Input model manifest has invalid or empty rows')
        groups = json.loads(local(args.group_plan).read_text())
        expected, row_groups = {}, {}
        for group in groups:
            name = group['model_group']
            if name in expected or len(group['indices']) != 5:
                raise ValueError('Invalid or duplicate model group in the plan')
            expected[name] = group['cell_type']
            for fold, index in enumerate(group['indices']):
                if type(index) is not int or index < 0 or index >= len(rows) or index in row_groups:
                    raise ValueError('Invalid or duplicate manifest row index in the plan')
                if rows[index]['model_id'] != f'{name}_fold{fold}' or rows[index]['cell_type'] != group['cell_type']:
                    raise ValueError('Group plan differs from the original model manifest')
                row_groups[index] = name
        if len(row_groups) != len(rows):
            raise ValueError('Group plan must include every model manifest row')
        results = {}
        for text in args.result_row:
            values = text.split('\t')
            if len(values) != len(result_fields):
                raise ValueError('MoDISco result row has the wrong number of columns')
            row = dict(zip(result_fields, values))
            if not valid_values(row):
                raise ValueError('MoDISco result values must be nonempty single-line strings')
            name = row['model_group']
            if name not in expected or name in results:
                raise ValueError(f'Unexpected or duplicate MoDISco model group: {name}')
            if row['cell_type'] != expected[name] or row['head'] != args.head:
                raise ValueError(f'MoDISco result cell type or head differs from the plan: {name}')
            results[name] = row
        if set(results) != set(expected):
            raise ValueError('MoDISco results must include every model group exactly once')

        def write_table(path, fields, table):
            with open(path, 'w', newline='') as stream:
                writer = csv.DictWriter(stream, fields, delimiter='\t', quoting=csv.QUOTE_NONE, quotechar=None,
                                        lineterminator='\n', extrasaction='ignore')
                writer.writeheader()
                writer.writerows(table)

        # Result paths are metadata; do not open, resolve, or rewrite these URIs.
        enriched = [dict(row, **{field: results[row_groups[index]][field]
                    for field in ['model_group', 'head'] + file_fields}) for index, row in enumerate(rows)]
        write_table('models.with_modisco.tsv', input_fields + ['model_group', 'head'] + file_fields, enriched)
        write_table('modisco_manifest.tsv', result_fields, [results[group['model_group']] for group in groups])
        print(f'[manifest] Wrote {len(rows)} model rows and {len(results)} MoDISco result rows', flush=True)
        PY
        echo '[manifest] Output manifests complete'
    >>>
    output {
        File models = 'models.with_modisco.tsv'
        File results = 'modisco_manifest.tsv'
        File log = 'motif_manifests.log'
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
