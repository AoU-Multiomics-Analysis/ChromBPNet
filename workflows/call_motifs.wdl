version 1.0

workflow ChromBPNetMotifCalls {
    input {
        File contributions
        File regions
        File motifs
        File motif_provenance
        File? tf_matches
        String cell_type
        String head = "counts"
        String docker_image
        Int region_width = 0
        String mode = "pp"
        Float global_lambda = 0.7
        Float trim_threshold = 0.3
        Int batch_size = 256
        Int max_steps = 10000
        Int memory_gb = 64
        Int disk_gb = 100
    }
    call PrepareFinemo {
        input: contributions=contributions, regions=regions, motifs=motifs,
               motif_provenance=motif_provenance, cell_type=cell_type, head=head,
               kind="peaks", region_width=region_width, docker_image=docker_image,
               memory_gb=memory_gb, disk_gb=disk_gb
    }
    call CallFinemo {
        input: regions=PrepareFinemo.npz, motifs=motifs, preparation=PrepareFinemo.metadata,
               mode=mode, global_lambda=global_lambda, trim_threshold=trim_threshold,
               batch_size=batch_size, max_steps=max_steps, docker_image=docker_image,
               memory_gb=memory_gb, disk_gb=disk_gb
    }
    call SummarizeFinemo {
        input: hits=CallFinemo.hits, qc=CallFinemo.qc, regions=regions,
               preparation=PrepareFinemo.metadata, call_metadata=CallFinemo.metadata,
               tf_matches=tf_matches, docker_image=docker_image,
               memory_gb=memory_gb, disk_gb=disk_gb
    }
    output {
        File annotated_hits = SummarizeFinemo.annotated_hits
        File motif_summary = SummarizeFinemo.motif_summary
        File report_html = SummarizeFinemo.html
        File report_bundle = SummarizeFinemo.bundle
        File raw_hits = CallFinemo.hits
        File unique_hits = CallFinemo.unique_hits
        File hits_bed = CallFinemo.bed
        File optimizer_qc = CallFinemo.qc
        File zero_contribution_rows = CallFinemo.zero_rows
        File raw_call_bundle = CallFinemo.bundle
        File prepared_npz = PrepareFinemo.npz
        File preparation_metadata = PrepareFinemo.metadata
        File call_metadata = CallFinemo.metadata
        File summary_metadata = SummarizeFinemo.metadata
        Array[File] logs = [PrepareFinemo.log, CallFinemo.log, SummarizeFinemo.log]
    }
}

task PrepareFinemo {
    input {
        File contributions
        File regions
        File motifs
        File motif_provenance
        File? coordinates
        String cell_type
        String head = "counts"
        String kind = "peaks"
        Int region_width = 0
        String docker_image
        Int memory_gb = 16
        Int disk_gb = 100
    }
    command <<<
        set -euo pipefail
        exec > >(tee preparefinemo.log) 2>&1
        echo '[Fi-NeMo] Start PrepareFinemo'
        python /opt/finemo_pipeline/scripts/prepare_finemo.py \
            --contributions '~{sub(contributions, "'", "'\"'\"'")}' \
            --regions '~{sub(regions, "'", "'\"'\"'")}' \
            --motifs '~{sub(motifs, "'", "'\"'\"'")}' \
            --motif-provenance '~{sub(motif_provenance, "'", "'\"'\"'")}' \
            ~{if defined(coordinates) then "--coordinates '" + sub(select_first([coordinates]), "'", "'\"'\"'") + "'" else ""} \
            --cell-type '~{sub(cell_type, "'", "'\"'\"'")}' \
            --head '~{sub(head, "'", "'\"'\"'")}' \
            --kind '~{sub(kind, "'", "'\"'\"'")}' \
            --region-width ~{region_width} \
            --output-dir prepared
        echo '[Fi-NeMo] Complete PrepareFinemo'
    >>>
    output {
        File npz = 'prepared/regions.npz'
        File metadata = 'prepared/metadata.json'
        File log = 'preparefinemo.log'
    }
    runtime {
        docker: docker_image
        cpu: 4
        memory: '~{memory_gb}GB'
        disks: 'local-disk ~{disk_gb} SSD'
        bootDiskSizeGb: 50
        preemptible: 0
    }
}

task CallFinemo {
    input {
        File regions
        File motifs
        File preparation
        String mode = "pp"
        Float global_lambda = 0.7
        Float trim_threshold = 0.3
        Int batch_size = 256
        Int max_steps = 10000
        String docker_image
        Int memory_gb = 64
        Int disk_gb = 100
    }
    command <<<
        set -euo pipefail
        exec > >(tee callfinemo.log) 2>&1
        echo '[Fi-NeMo] Start CallFinemo'
        nvidia-smi
        python /opt/finemo_pipeline/scripts/call_finemo.py \
            --regions '~{sub(regions, "'", "'\"'\"'")}' \
            --motifs '~{sub(motifs, "'", "'\"'\"'")}' \
            --preparation '~{sub(preparation, "'", "'\"'\"'")}' \
            --mode '~{sub(mode, "'", "'\"'\"'")}' \
            --global-lambda ~{global_lambda} \
            --trim-threshold ~{trim_threshold} \
            --batch-size ~{batch_size} \
            --max-steps ~{max_steps} \
            --threads 16 \
            --output-dir calls
        echo '[Fi-NeMo] Complete CallFinemo'
    >>>
    output {
        File hits = 'calls/hits.tsv'
        File unique_hits = 'calls/hits_unique.tsv'
        File bed = 'calls/hits.bed'
        File qc = 'calls/peaks_qc.tsv'
        File zero_rows = 'calls/zero_contribution_rows.tsv'
        File bundle = 'calls/raw_calls.zip'
        File metadata = 'calls/metadata.json'
        File log = 'callfinemo.log'
    }
    runtime {
        docker: docker_image
        cpu: 16
        memory: '~{memory_gb}GB'
        disks: 'local-disk ~{disk_gb} SSD'
        bootDiskSizeGb: 50
        preemptible: 0
        predefinedMachineType: 'g2-standard-16'
        gpuType: 'nvidia-l4'
        gpuCount: 1
        zones: ['us-central1-a', 'us-central1-b', 'us-central1-c', 'us-central1-f']
    }
}

task SummarizeFinemo {
    input {
        File hits
        File qc
        File regions
        File preparation
        File call_metadata
        File? coordinates
        File? tf_matches
        Int pairing_tolerance = 3
        String docker_image
        Int memory_gb = 16
        Int disk_gb = 100
    }
    command <<<
        set -euo pipefail
        exec > >(tee summarizefinemo.log) 2>&1
        echo '[Fi-NeMo] Start SummarizeFinemo'
        python /opt/finemo_pipeline/scripts/summarize_finemo.py \
            --hits '~{sub(hits, "'", "'\"'\"'")}' \
            --qc '~{sub(qc, "'", "'\"'\"'")}' \
            --regions '~{sub(regions, "'", "'\"'\"'")}' \
            --preparation '~{sub(preparation, "'", "'\"'\"'")}' \
            --call-metadata '~{sub(call_metadata, "'", "'\"'\"'")}' \
            ~{if defined(coordinates) then "--coordinates '" + sub(select_first([coordinates]), "'", "'\"'\"'") + "'" else ""} \
            ~{if defined(tf_matches) then "--tf-matches '" + sub(select_first([tf_matches]), "'", "'\"'\"'") + "'" else ""} \
            --pairing-tolerance ~{pairing_tolerance} \
            --output-dir summary
        echo '[Fi-NeMo] Complete SummarizeFinemo'
    >>>
    output {
        File annotated_hits = 'summary/annotated_hits.tsv'
        File changes = 'summary/motif_changes.tsv'
        File variants = 'summary/variant_summary.tsv'
        File motif_summary = 'summary/motif_summary.tsv'
        File html = 'summary/report.html'
        File bundle = 'summary/report.zip'
        File metadata = 'summary/metadata.json'
        File log = 'summarizefinemo.log'
    }
    runtime {
        docker: docker_image
        cpu: 4
        memory: '~{memory_gb}GB'
        disks: 'local-disk ~{disk_gb} SSD'
        bootDiskSizeGb: 50
        preemptible: 0
    }
}
