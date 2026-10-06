version 1.0

struct ModelSpec {
    String model_id
    String cell_type
    File model
    File peaks
}

workflow ChromBPNetVariantScoring {
    input {
        File variants
        File model_manifest
        File genome
        File? genome_index
        File chrom_sizes
        String docker_image
        Int batch_size = 128
        Int num_shuf = 0
        Int? max_peaks
        Int random_seed = 1234
        Int scoring_memory_gb = 64
        Int scoring_disk_gb = 100
        Int num_preempt = 0
    }

    call ValidateManifest {
        input:
            manifest = model_manifest,
            docker_image = docker_image
    }

    Array[Array[String]] manifest_rows = read_tsv(ValidateManifest.rows)

    scatter (row in manifest_rows) {
        # Coerce URI metadata to typed File fields before task localization.
        ModelSpec entry = object {
            model_id: row[0],
            cell_type: row[1],
            model: row[2],
            peaks: row[3]
        }
        call ScoreVariants {
            input:
                model_id = entry.model_id,
                cell_type = entry.cell_type,
                model = entry.model,
                peaks = entry.peaks,
                variants = variants,
                genome = genome,
                genome_index = genome_index,
                chrom_sizes = chrom_sizes,
                manifest_validation = ValidateManifest.validation,
                docker_image = docker_image,
                batch_size = batch_size,
                num_shuf = num_shuf,
                max_peaks = max_peaks,
                random_seed = random_seed,
                memory_gb = scoring_memory_gb,
                disk_gb = scoring_disk_gb,
                num_preempt = num_preempt
        }
    }

    call MergeVariantEffects {
        input:
            score_files = ScoreVariants.variant_effects,
            docker_image = docker_image
    }

    output {
        Array[File] per_model_effects = ScoreVariants.variant_effects
        Array[File] per_model_peak_scores = ScoreVariants.peak_scores
        Array[File] per_model_metadata = ScoreVariants.metadata
        Array[File] per_model_logs = ScoreVariants.log
        Array[File] shuffled_scores = flatten(ScoreVariants.shuffled_scores)
        File merged_effects = MergeVariantEffects.long_effects
        File wide_effects = MergeVariantEffects.wide_effects
        File merge_log = MergeVariantEffects.log
        File manifest_validation = ValidateManifest.validation
    }
}

task ValidateManifest {
    input {
        File manifest
        String docker_image
    }
    command <<<
        set -euo pipefail
        echo '[manifest] Start input checks'
        python /opt/chrombpnet/scripts/validate_manifest.py \
            --manifest '~{sub(manifest, "'", "'\"'\"'")}' \
            --output manifest.validated.txt \
            --rows-output models.normalized.tsv
        echo '[manifest] Input checks complete'
    >>>
    output {
        File validation = 'manifest.validated.txt'
        File rows = 'models.normalized.tsv'
    }
    runtime {
        docker: docker_image
        cpu: 1
        memory: '2 GB'
        disks: 'local-disk 10 SSD'
        bootDiskSizeGb: 50
        preemptible: 0
    }
}

task ScoreVariants {
    input {
        String model_id
        String cell_type
        File model
        File peaks
        File variants
        File genome
        File? genome_index
        File chrom_sizes
        File manifest_validation
        String docker_image
        Int batch_size
        Int num_shuf
        Int? max_peaks
        Int random_seed
        Int memory_gb
        Int disk_gb
        Int num_preempt
    }
    command <<<
        set -euo pipefail
        exec > >(tee score.log) 2>&1
        echo '[score] Check manifest validation and GPU driver'
        test -s '~{sub(manifest_validation, "'", "'\"'\"'")}'
        nvidia-smi
        echo '[score] Start variant effects'
        python /opt/chrombpnet/scripts/score_variants.py \
            --model-id '~{sub(model_id, "'", "'\"'\"'")}' \
            --cell-type '~{sub(cell_type, "'", "'\"'\"'")}' \
            --model '~{sub(model, "'", "'\"'\"'")}' \
            --peaks '~{sub(peaks, "'", "'\"'\"'")}' \
            --variants '~{sub(variants, "'", "'\"'\"'")}' \
            --genome '~{sub(genome, "'", "'\"'\"'")}' \
            ~{if defined(genome_index) then "--genome-index '" + sub(select_first([genome_index]), "'", "'\"'\"'") + "'" else ""} \
            --chrom-sizes '~{sub(chrom_sizes, "'", "'\"'\"'")}' \
            --batch-size ~{batch_size} \
            --num-shuf ~{num_shuf} \
            --random-seed ~{random_seed} \
            --threads 16 \
            ~{if defined(max_peaks) then "--max-peaks " + select_first([max_peaks]) else ""} \
            --output-dir effects
        echo '[score] Variant effects complete'
    >>>
    output {
        File variant_effects = 'effects/variant_effects.tsv'
        File peak_scores = 'effects/scored.peak_scores.tsv'
        File metadata = 'effects/run_metadata.json'
        File log = 'score.log'
        Array[File] shuffled_scores = glob('effects/scored.variant_scores.shuffled.tsv')
    }
    runtime {
        docker: docker_image
        memory: '~{memory_gb}GB'
        disks: 'local-disk ~{disk_gb} SSD'
        bootDiskSizeGb: 50
        cpu: 16
        preemptible: num_preempt
        predefinedMachineType: 'g2-standard-16'
        gpuType: 'nvidia-l4'
        gpuCount: 1
        zones: ['us-central1-a', 'us-central1-b', 'us-central1-c', 'us-central1-f']
    }
}

task MergeVariantEffects {
    input {
        Array[File] score_files
        String docker_image
    }
    command <<<
        set -euo pipefail
        exec > >(tee merge.log) 2>&1
        echo '[merge] Start cell-type merge'
        # Create the list in the execution directory from localized File inputs.
        # A write_lines result inside a String expression can remain a cloud URI.
        printf '%s\n' '~{sub(sep("\n", score_files), "'", "'\"'\"'")}' > score_files.list
        echo '[merge] Created task-local score file list'
        python /opt/chrombpnet/scripts/merge_scores.py \
            --score-files score_files.list \
            --long-output variant_effects.all_models.tsv \
            --wide-output variant_effects.wide.tsv
        echo '[merge] Cell-type merge complete'
    >>>
    output {
        File long_effects = 'variant_effects.all_models.tsv'
        File wide_effects = 'variant_effects.wide.tsv'
        File log = 'merge.log'
    }
    runtime {
        docker: docker_image
        cpu: 2
        memory: '8 GB'
        disks: 'local-disk 20 SSD'
        bootDiskSizeGb: 50
        preemptible: 0
    }
}
