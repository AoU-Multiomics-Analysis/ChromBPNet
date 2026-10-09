version 1.0

workflow ChromBPNetMotifDiscovery {
    input {
        Array[File] fold_models
        File peaks
        File genome
        File? genome_index
        File chrom_sizes
        File motif_database
        String cell_type
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

    call PrepareMotifInputs {
        input:
            fold_models = fold_models,
            peaks = peaks,
            genome = genome,
            genome_index = genome_index,
            chrom_sizes = chrom_sizes,
            cell_type = cell_type,
            head = head,
            input_length = input_length,
            discovery_window = discovery_window,
            random_seed = random_seed,
            max_peaks = max_peaks,
            docker_image = docker_image
    }

    scatter (fold_index in range(length(fold_models))) {
        call FoldContributions {
            input:
                model = fold_models[fold_index],
                fold_index = fold_index,
                peaks = PrepareMotifInputs.prepared_peaks,
                genome = PrepareMotifInputs.reference,
                genome_index = PrepareMotifInputs.reference_index,
                preparation = PrepareMotifInputs.metadata,
                cell_type = cell_type,
                head = head,
                random_seed = random_seed,
                num_backgrounds = num_backgrounds,
                batch_size = batch_size,
                threads = 16,
                memory_gb = contribution_memory_gb,
                disk_gb = contribution_disk_gb,
                num_preempt = num_preempt,
                docker_image = docker_image
        }
    }

    call AverageContributions {
        input:
            score_files = FoldContributions.scores,
            region_files = FoldContributions.regions,
            expected_folds = 5,
            chunk_rows = batch_size,
            memory_gb = 8,
            disk_gb = averaging_disk_gb,
            docker_image = docker_image
    }

    call DiscoverMotifs {
        input:
            contributions = AverageContributions.averaged,
            window = discovery_window,
            max_seqlets = max_seqlets,
            n_leiden = n_leiden,
            random_seed = random_seed,
            cpu_count = motif_cpu,
            memory_gb = motif_memory_gb,
            disk_gb = motif_disk_gb,
            docker_image = docker_image
    }

    call ReportMotifs {
        input:
            motifs = DiscoverMotifs.motifs,
            motif_database = motif_database,
            match_qvalue = match_qvalue,
            n_matches = n_matches,
            docker_image = docker_image
    }

    output {
        File modisco_motifs = DiscoverMotifs.motifs
        File averaged_contributions = AverageContributions.averaged
        File interpreted_regions = AverageContributions.regions
        File report_html = ReportMotifs.html
        File report_bundle = ReportMotifs.bundle
        File discovered_motifs_meme = ReportMotifs.meme
        File candidate_tf_matches = ReportMotifs.matches
        File tomtom_results = ReportMotifs.tomtom_results
        File motif_inventory = ReportMotifs.inventory
        File preparation_metadata = PrepareMotifInputs.metadata
        File averaging_metadata = AverageContributions.metadata
        File discovery_metadata = DiscoverMotifs.metadata
        File report_metadata = ReportMotifs.metadata
        Array[File] per_fold_contributions = FoldContributions.scores
        Array[File] per_fold_metadata = FoldContributions.metadata
        Array[File] per_fold_logs = FoldContributions.log
        File preparation_log = PrepareMotifInputs.log
        File averaging_log = AverageContributions.log
        File discovery_log = DiscoverMotifs.log
        File report_log = ReportMotifs.log
    }
}

task PrepareMotifInputs {
    input {
        Array[File] fold_models
        File peaks
        File genome
        File? genome_index
        File chrom_sizes
        String cell_type
        String head
        Int input_length
        Int discovery_window
        Int random_seed
        Int? max_peaks
        String docker_image
    }
    Boolean five_fold_models = length(fold_models) == 5
    command <<<
        set -euo pipefail
        exec > >(tee prepare.log) 2>&1
        echo '[prepare] Start five-fold validation and common peak preparation'
        python /opt/motif_pipeline/scripts/write_file_list.py --output models.list \
            ~{if five_fold_models then "--file '" + sub(fold_models[0], "'", "'\"'\"'") + "'" else ""} \
            ~{if five_fold_models then "--file '" + sub(fold_models[1], "'", "'\"'\"'") + "'" else ""} \
            ~{if five_fold_models then "--file '" + sub(fold_models[2], "'", "'\"'\"'") + "'" else ""} \
            ~{if five_fold_models then "--file '" + sub(fold_models[3], "'", "'\"'\"'") + "'" else ""} \
            ~{if five_fold_models then "--file '" + sub(fold_models[4], "'", "'\"'\"'") + "'" else ""}
        python /opt/motif_pipeline/scripts/prepare_motif_inputs.py \
            --models-list models.list \
            --peaks '~{sub(peaks, "'", "'\"'\"'")}' \
            --genome '~{sub(genome, "'", "'\"'\"'")}' \
            ~{if defined(genome_index) then "--genome-index '" + sub(select_first([genome_index]), "'", "'\"'\"'") + "'" else ""} \
            --chrom-sizes '~{sub(chrom_sizes, "'", "'\"'\"'")}' \
            --cell-type '~{sub(cell_type, "'", "'\"'\"'")}' \
            --head '~{sub(head, "'", "'\"'\"'")}' \
            --input-length ~{input_length} \
            --discovery-window ~{discovery_window} \
            --expected-folds 5 \
            --random-seed ~{random_seed} \
            ~{if defined(max_peaks) then "--max-peaks " + select_first([max_peaks]) else ""} \
            --output-dir prepared
        echo '[prepare] Complete'
    >>>
    output {
        File prepared_peaks = 'prepared/peaks.bed'
        File reference = 'prepared/reference.fa'
        File reference_index = 'prepared/reference.fa.fai'
        File metadata = 'prepared/preparation.json'
        File log = 'prepare.log'
    }
    runtime {
        docker: docker_image
        cpu: 2
        memory: '16 GB'
        disks: 'local-disk 100 SSD'
        bootDiskSizeGb: 50
        preemptible: 0
    }
}

task FoldContributions {
    input {
        File model
        File peaks
        File genome
        File genome_index
        File preparation
        Int fold_index
        String cell_type
        String head
        Int random_seed
        Int num_backgrounds
        Int batch_size
        Int threads
        Int memory_gb
        Int disk_gb
        Int num_preempt
        String docker_image
    }
    command <<<
        set -euo pipefail
        exec > >(tee contributions.log) 2>&1
        echo '[contributions] Check GPU driver and start peak contributions'
        nvidia-smi
        python /opt/motif_pipeline/scripts/fold_contributions.py \
            --model '~{sub(model, "'", "'\"'\"'")}' \
            --peaks '~{sub(peaks, "'", "'\"'\"'")}' \
            --genome '~{sub(genome, "'", "'\"'\"'")}' \
            --genome-index '~{sub(genome_index, "'", "'\"'\"'")}' \
            --preparation '~{sub(preparation, "'", "'\"'\"'")}' \
            --cell-type '~{sub(cell_type, "'", "'\"'\"'")}' \
            --head '~{sub(head, "'", "'\"'\"'")}' \
            --fold-index ~{fold_index} \
            --random-seed ~{random_seed} \
            --num-backgrounds ~{num_backgrounds} \
            --batch-size ~{batch_size} \
            --threads ~{threads} \
            --output-dir contributions
        echo '[contributions] Complete'
    >>>
    output {
        File scores = 'contributions/scores.h5'
        File regions = 'contributions/regions.bed'
        File metadata = 'contributions/metadata.json'
        File log = 'contributions.log'
    }
    runtime {
        docker: docker_image
        cpu: threads
        memory: '~{memory_gb}GB'
        disks: 'local-disk ~{disk_gb} SSD'
        bootDiskSizeGb: 50
        preemptible: num_preempt
        predefinedMachineType: 'g2-standard-16'
        gpuType: 'nvidia-l4'
        gpuCount: 1
        zones: ['us-central1-a', 'us-central1-b', 'us-central1-c', 'us-central1-f']
    }
}

task AverageContributions {
    input {
        Array[File] score_files
        Array[File] region_files
        Int expected_folds
        Int chunk_rows
        Int memory_gb
        Int disk_gb
        String docker_image
    }
    Boolean five_score_files = length(score_files) == 5
    Boolean five_region_files = length(region_files) == 5
    command <<<
        set -euo pipefail
        exec > >(tee average.log) 2>&1
        echo '[average] Create task-local lists from localized File inputs'
        python /opt/motif_pipeline/scripts/write_file_list.py --output scores.list \
            ~{if five_score_files then "--file '" + sub(score_files[0], "'", "'\"'\"'") + "'" else ""} \
            ~{if five_score_files then "--file '" + sub(score_files[1], "'", "'\"'\"'") + "'" else ""} \
            ~{if five_score_files then "--file '" + sub(score_files[2], "'", "'\"'\"'") + "'" else ""} \
            ~{if five_score_files then "--file '" + sub(score_files[3], "'", "'\"'\"'") + "'" else ""} \
            ~{if five_score_files then "--file '" + sub(score_files[4], "'", "'\"'\"'") + "'" else ""}
        python /opt/motif_pipeline/scripts/write_file_list.py --output regions.list \
            ~{if five_region_files then "--file '" + sub(region_files[0], "'", "'\"'\"'") + "'" else ""} \
            ~{if five_region_files then "--file '" + sub(region_files[1], "'", "'\"'\"'") + "'" else ""} \
            ~{if five_region_files then "--file '" + sub(region_files[2], "'", "'\"'\"'") + "'" else ""} \
            ~{if five_region_files then "--file '" + sub(region_files[3], "'", "'\"'\"'") + "'" else ""} \
            ~{if five_region_files then "--file '" + sub(region_files[4], "'", "'\"'\"'") + "'" else ""}
        python /opt/motif_pipeline/scripts/average_contributions.py \
            --scores-list scores.list \
            --regions-list regions.list \
            --expected-folds ~{expected_folds} \
            --chunk-rows ~{chunk_rows} \
            --output averaged.h5 \
            --output-regions regions.bed \
            --metadata averaging.json
        echo '[average] Complete'
    >>>
    output {
        File averaged = 'averaged.h5'
        File regions = 'regions.bed'
        File metadata = 'averaging.json'
        File log = 'average.log'
    }
    runtime {
        docker: docker_image
        cpu: 2
        memory: '~{memory_gb}GB'
        disks: 'local-disk ~{disk_gb} SSD'
        bootDiskSizeGb: 50
        preemptible: 0
    }
}

task DiscoverMotifs {
    input {
        File contributions
        Int window
        Int max_seqlets
        Int n_leiden
        Int random_seed
        Int cpu_count
        Int memory_gb
        Int disk_gb
        String docker_image
    }
    command <<<
        set -euo pipefail
        exec > >(tee modisco.log) 2>&1
        export OMP_NUM_THREADS=~{cpu_count}
        export NUMBA_NUM_THREADS=~{cpu_count}
        echo '[modisco] Start TF-MoDISco on averaged scores'
        python /opt/motif_pipeline/scripts/discover_motifs.py \
            --contributions '~{sub(contributions, "'", "'\"'\"'")}' \
            --window ~{window} \
            --max-seqlets ~{max_seqlets} \
            --n-leiden ~{n_leiden} \
            --random-seed ~{random_seed} \
            --output modisco_results.h5 \
            --metadata discovery.json
        echo '[modisco] Complete'
    >>>
    output {
        File motifs = 'modisco_results.h5'
        File metadata = 'discovery.json'
        File log = 'modisco.log'
    }
    runtime {
        docker: docker_image
        cpu: cpu_count
        memory: '~{memory_gb}GB'
        disks: 'local-disk ~{disk_gb} SSD'
        bootDiskSizeGb: 50
        preemptible: 0
    }
}

task ReportMotifs {
    input {
        File motifs
        File motif_database
        Float match_qvalue
        Int n_matches
        String docker_image
    }
    command <<<
        set -euo pipefail
        exec > >(tee report.log) 2>&1
        echo '[report] Start known-motif matching and HTML report'
        python /opt/motif_pipeline/scripts/report_motifs.py \
            --motifs '~{sub(motifs, "'", "'\"'\"'")}' \
            --motif-database '~{sub(motif_database, "'", "'\"'\"'")}' \
            --match-qvalue ~{match_qvalue} \
            --n-matches ~{n_matches} \
            --output-dir annotation
        echo '[report] Complete'
    >>>
    output {
        File html = 'annotation/report/report.html'
        File bundle = 'annotation/report_bundle.zip'
        File meme = 'annotation/discovered_motifs.meme'
        File matches = 'annotation/candidate_tf_matches.tsv'
        File inventory = 'annotation/motif_inventory.tsv'
        File tomtom_results = 'annotation/tomtom/tomtom.tsv'
        File metadata = 'annotation/report_metadata.json'
        File log = 'report.log'
    }
    runtime {
        docker: docker_image
        cpu: 4
        memory: '16 GB'
        disks: 'local-disk 30 SSD'
        bootDiskSizeGb: 50
        preemptible: 0
    }
}
