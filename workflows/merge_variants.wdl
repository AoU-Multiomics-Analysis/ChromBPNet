version 1.0

import "summarize_variants.wdl" as summaries

# Recover a failed merge using existing per-model score files. No scoring calls.
workflow ChromBPNetMergeVariantEffects {
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

    call summaries.ChromBPNetSummarizeVariants as SummarizeVariants {
        input:
            model_manifest = model_manifest,
            score_files = score_files,
            docker_image = docker_image,
            summary_memory_gb = summary_memory_gb,
            summary_disk_gb = summary_disk_gb,
            summary_max_retries = summary_max_retries,
            merge_memory_gb = merge_memory_gb,
            merge_disk_gb = merge_disk_gb,
            merge_max_retries = merge_max_retries
    }

    output {
        File merged_effects = SummarizeVariants.merged_effects
        File wide_effects = SummarizeVariants.wide_effects
        File fold_summary = SummarizeVariants.fold_summary
        File merge_log = SummarizeVariants.merge_log
        Array[File] per_model_fold_summaries = SummarizeVariants.per_model_fold_summaries
        Array[File] summary_logs = SummarizeVariants.summary_logs
        File grouping_log = SummarizeVariants.grouping_log
    }
}
