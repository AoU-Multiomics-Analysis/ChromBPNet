version 1.0

import "score_variants.wdl" as scoring

# Recover a failed merge using existing per-model score files. No scoring calls.
workflow ChromBPNetMergeVariantEffects {
    input {
        Array[File] score_files
        String docker_image
        Int merge_memory_gb = 64
        Int merge_disk_gb = 500
        Int merge_max_retries = 2
    }

    call scoring.MergeVariantEffects {
        input:
            score_files = score_files,
            docker_image = docker_image,
            memory_gb = merge_memory_gb,
            disk_gb = merge_disk_gb,
            max_retries = merge_max_retries
    }

    output {
        File merged_effects = MergeVariantEffects.long_effects
        File wide_effects = MergeVariantEffects.wide_effects
        File fold_summary = MergeVariantEffects.fold_summary
        File merge_log = MergeVariantEffects.log
    }
}
