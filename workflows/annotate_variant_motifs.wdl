version 1.0
import "discover_motifs.wdl" as motif
import "call_motifs.wdl" as finemo

workflow ChromBPNetVariantMotifs {
    input {
        Array[File] fold_models
        File variants
        File genome
        File? genome_index
        File chrom_sizes
        File motifs
        File motif_provenance
        File? tf_matches
        String cell_type
        String head = "counts"
        String contribution_image
        String finemo_image
        Int input_length = 2114
        Int shard_size = 500
        Int num_backgrounds = 20
        Int random_seed = 1234
        Int contribution_batch_size = 64
        Int contribution_memory_gb = 64
        Int contribution_disk_gb = 100
        Int region_width = 0
        String mode = "pp"
        Float global_lambda = 0.7
        Float trim_threshold = 0.3
        Int finemo_batch_size = 256
        Int max_steps = 10000
        Int finemo_memory_gb = 64
        Int finemo_disk_gb = 100
        Int pairing_tolerance = 3
    }
    call PrepareVariantInputs {
        input: fold_models=fold_models, variants=variants, genome=genome,
               genome_index=genome_index, chrom_sizes=chrom_sizes, motif_provenance=motif_provenance,
               cell_type=cell_type, head=head, input_length=input_length, shard_size=shard_size,
               docker_image=contribution_image, disk_gb=contribution_disk_gb
    }
    scatter (shard in range(length(PrepareVariantInputs.sequences))) {
        scatter (fold in range(length(fold_models))) {
            call VariantContributions {
                input: model=fold_models[fold], sequences=PrepareVariantInputs.sequences[shard],
                       regions=PrepareVariantInputs.regions[shard], preparation=PrepareVariantInputs.preparations[shard],
                       fold_index=fold, cell_type=cell_type, head=head, num_backgrounds=num_backgrounds,
                       random_seed=random_seed, batch_size=contribution_batch_size,
                       docker_image=contribution_image, memory_gb=contribution_memory_gb, disk_gb=contribution_disk_gb
            }
        }
        call motif.AverageContributions {
            input: score_files=VariantContributions.scores, region_files=VariantContributions.copied_regions,
                   expected_folds=5, chunk_rows=contribution_batch_size, memory_gb=16,
                   disk_gb=contribution_disk_gb, docker_image=contribution_image
        }
        call finemo.PrepareFinemo {
            input: contributions=AverageContributions.averaged, regions=AverageContributions.regions,
                   motifs=motifs, motif_provenance=motif_provenance,
                   coordinates=PrepareVariantInputs.sequences[shard], kind="variants",
                   cell_type=cell_type, head=head, region_width=region_width,
                   docker_image=finemo_image, memory_gb=finemo_memory_gb, disk_gb=finemo_disk_gb
        }
        call finemo.CallFinemo {
            input: regions=PrepareFinemo.npz, motifs=motifs, preparation=PrepareFinemo.metadata,
                   mode=mode, global_lambda=global_lambda, trim_threshold=trim_threshold,
                   batch_size=finemo_batch_size, max_steps=max_steps,
                   docker_image=finemo_image, memory_gb=finemo_memory_gb, disk_gb=finemo_disk_gb
        }
        call finemo.SummarizeFinemo {
            input: hits=CallFinemo.hits, qc=CallFinemo.qc, regions=AverageContributions.regions,
                   preparation=PrepareFinemo.metadata, call_metadata=CallFinemo.metadata,
                   coordinates=PrepareVariantInputs.sequences[shard], tf_matches=tf_matches,
                   pairing_tolerance=pairing_tolerance, docker_image=finemo_image,
                   memory_gb=finemo_memory_gb, disk_gb=finemo_disk_gb
        }
    }
    call MergeFinemo {
        input: hit_files=SummarizeFinemo.annotated_hits, change_files=SummarizeFinemo.changes,
               variant_files=SummarizeFinemo.variants, docker_image=finemo_image,
               memory_gb=finemo_memory_gb, disk_gb=finemo_disk_gb
    }
    output {
        File annotated_hits = MergeFinemo.hits
        File motif_changes = MergeFinemo.changes
        File variant_summary = MergeFinemo.variants
        File merge_metadata = MergeFinemo.metadata
        File preparation_metadata = PrepareVariantInputs.metadata
        Array[File] allele_sequences = PrepareVariantInputs.sequences
        Array[File] allele_rows = PrepareVariantInputs.regions
        Array[File] averaged_contributions = AverageContributions.averaged
        Array[File] report_bundles = SummarizeFinemo.bundle
        Array[File] raw_call_bundles = CallFinemo.bundle
        Array[File] optimizer_qc = CallFinemo.qc
        Array[File] zero_contribution_rows = CallFinemo.zero_rows
        Array[File] call_metadata = CallFinemo.metadata
        Array[File] summary_metadata = SummarizeFinemo.metadata
        Array[File] score_logs = flatten(VariantContributions.log)
        Array[File] call_logs = CallFinemo.log
        Array[File] summary_logs = SummarizeFinemo.log
        File preparation_log = PrepareVariantInputs.log
        File merge_log = MergeFinemo.log
    }
}

task PrepareVariantInputs {
    input {
        Array[File] fold_models
        File variants
        File genome
        File? genome_index
        File chrom_sizes
        File motif_provenance
        String cell_type
        String head = "counts"
        Int input_length = 2114
        Int shard_size = 500
        String docker_image
        Int memory_gb = 16
        Int disk_gb = 100
    }
    Boolean five_models = length(fold_models) == 5
    command <<<
        set -euo pipefail
        exec > >(tee preparevariantinputs.log) 2>&1
        echo '[Fi-NeMo] Start PrepareVariantInputs'
        python /opt/motif_pipeline/scripts/write_file_list.py --output models.list \
            ~{if five_models then "--file '" + sub(fold_models[0], "'", "'\"'\"'") + "'" else ""} \
            ~{if five_models then "--file '" + sub(fold_models[1], "'", "'\"'\"'") + "'" else ""} \
            ~{if five_models then "--file '" + sub(fold_models[2], "'", "'\"'\"'") + "'" else ""} \
            ~{if five_models then "--file '" + sub(fold_models[3], "'", "'\"'\"'") + "'" else ""} \
            ~{if five_models then "--file '" + sub(fold_models[4], "'", "'\"'\"'") + "'" else ""}
        python /opt/motif_pipeline/scripts/prepare_variant_motifs.py \
            --models-list models.list \
            --variants '~{sub(variants, "'", "'\"'\"'")}' \
            --genome '~{sub(genome, "'", "'\"'\"'")}' \
            ~{if defined(genome_index) then "--genome-index '" + sub(select_first([genome_index]), "'", "'\"'\"'") + "'" else ""} \
            --chrom-sizes '~{sub(chrom_sizes, "'", "'\"'\"'")}' \
            --motif-provenance '~{sub(motif_provenance, "'", "'\"'\"'")}' \
            --cell-type '~{sub(cell_type, "'", "'\"'\"'")}' \
            --head '~{sub(head, "'", "'\"'\"'")}' \
            --input-length ~{input_length} \
            --shard-size ~{shard_size} \
            --output-dir prepared
        echo '[Fi-NeMo] Complete PrepareVariantInputs'
    >>>
    output {
        Array[File] sequences = glob('prepared/shard_*.h5')
        Array[File] regions = glob('prepared/shard_*.tsv')
        Array[File] preparations = glob('prepared/shard_*.json')
        File metadata = 'prepared/metadata.json'
        File log = 'preparevariantinputs.log'
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

task VariantContributions {
    input {
        File model
        File sequences
        File regions
        File preparation
        String cell_type
        String head = "counts"
        Int fold_index
        Int random_seed = 1234
        Int num_backgrounds = 20
        Int batch_size = 64
        String docker_image
        Int memory_gb = 64
        Int disk_gb = 100
    }
    command <<<
        set -euo pipefail
        exec > >(tee variantcontributions.log) 2>&1
        echo '[Fi-NeMo] Start VariantContributions'
        nvidia-smi
        python /opt/motif_pipeline/scripts/variant_contributions.py \
            --model '~{sub(model, "'", "'\"'\"'")}' \
            --sequences '~{sub(sequences, "'", "'\"'\"'")}' \
            --regions '~{sub(regions, "'", "'\"'\"'")}' \
            --preparation '~{sub(preparation, "'", "'\"'\"'")}' \
            --cell-type '~{sub(cell_type, "'", "'\"'\"'")}' \
            --head '~{sub(head, "'", "'\"'\"'")}' \
            --fold-index ~{fold_index} \
            --random-seed ~{random_seed} \
            --num-backgrounds ~{num_backgrounds} \
            --batch-size ~{batch_size} \
            --threads 16 \
            --output-dir contributions
        echo '[Fi-NeMo] Complete VariantContributions'
    >>>
    output {
        File scores = 'contributions/scores.h5'
        File copied_regions = 'contributions/regions.tsv'
        File metadata = 'contributions/metadata.json'
        File log = 'variantcontributions.log'
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

task MergeFinemo {
    input {
        Array[File] hit_files
        Array[File] change_files
        Array[File] variant_files
        String docker_image
        Int memory_gb = 16
        Int disk_gb = 100
    }
    Boolean hits_0 = length(hit_files) > 0
    Boolean hits_1 = length(hit_files) > 1
    Boolean hits_2 = length(hit_files) > 2
    Boolean hits_3 = length(hit_files) > 3
    Boolean hits_4 = length(hit_files) > 4
    Boolean hits_5 = length(hit_files) > 5
    Boolean hits_6 = length(hit_files) > 6
    Boolean hits_7 = length(hit_files) > 7
    Boolean hits_8 = length(hit_files) > 8
    Boolean hits_9 = length(hit_files) > 9
    Boolean hits_10 = length(hit_files) > 10
    Boolean hits_11 = length(hit_files) > 11
    Boolean hits_12 = length(hit_files) > 12
    Boolean hits_13 = length(hit_files) > 13
    Boolean hits_14 = length(hit_files) > 14
    Boolean hits_15 = length(hit_files) > 15
    Boolean hits_16 = length(hit_files) > 16
    Boolean hits_17 = length(hit_files) > 17
    Boolean hits_18 = length(hit_files) > 18
    Boolean hits_19 = length(hit_files) > 19
    Boolean hits_20 = length(hit_files) > 20
    Boolean hits_21 = length(hit_files) > 21
    Boolean hits_22 = length(hit_files) > 22
    Boolean hits_23 = length(hit_files) > 23
    Boolean hits_24 = length(hit_files) > 24
    Boolean hits_25 = length(hit_files) > 25
    Boolean hits_26 = length(hit_files) > 26
    Boolean hits_27 = length(hit_files) > 27
    Boolean hits_28 = length(hit_files) > 28
    Boolean hits_29 = length(hit_files) > 29
    Boolean hits_30 = length(hit_files) > 30
    Boolean hits_31 = length(hit_files) > 31
    Boolean hits_32 = length(hit_files) > 32
    Boolean hits_33 = length(hit_files) > 33
    Boolean hits_34 = length(hit_files) > 34
    Boolean hits_35 = length(hit_files) > 35
    Boolean hits_36 = length(hit_files) > 36
    Boolean hits_37 = length(hit_files) > 37
    Boolean hits_38 = length(hit_files) > 38
    Boolean hits_39 = length(hit_files) > 39
    Boolean hits_40 = length(hit_files) > 40
    Boolean hits_41 = length(hit_files) > 41
    Boolean hits_42 = length(hit_files) > 42
    Boolean hits_43 = length(hit_files) > 43
    Boolean hits_44 = length(hit_files) > 44
    Boolean hits_45 = length(hit_files) > 45
    Boolean hits_46 = length(hit_files) > 46
    Boolean hits_47 = length(hit_files) > 47
    Boolean hits_48 = length(hit_files) > 48
    Boolean hits_49 = length(hit_files) > 49
    Boolean hits_50 = length(hit_files) > 50
    Boolean hits_51 = length(hit_files) > 51
    Boolean hits_52 = length(hit_files) > 52
    Boolean hits_53 = length(hit_files) > 53
    Boolean hits_54 = length(hit_files) > 54
    Boolean hits_55 = length(hit_files) > 55
    Boolean hits_56 = length(hit_files) > 56
    Boolean hits_57 = length(hit_files) > 57
    Boolean hits_58 = length(hit_files) > 58
    Boolean hits_59 = length(hit_files) > 59
    Boolean hits_60 = length(hit_files) > 60
    Boolean hits_61 = length(hit_files) > 61
    Boolean hits_62 = length(hit_files) > 62
    Boolean hits_63 = length(hit_files) > 63
    Boolean hits_64 = length(hit_files) > 64
    Boolean hits_65 = length(hit_files) > 65
    Boolean hits_66 = length(hit_files) > 66
    Boolean hits_67 = length(hit_files) > 67
    Boolean hits_68 = length(hit_files) > 68
    Boolean hits_69 = length(hit_files) > 69
    Boolean hits_70 = length(hit_files) > 70
    Boolean hits_71 = length(hit_files) > 71
    Boolean hits_72 = length(hit_files) > 72
    Boolean hits_73 = length(hit_files) > 73
    Boolean hits_74 = length(hit_files) > 74
    Boolean hits_75 = length(hit_files) > 75
    Boolean hits_76 = length(hit_files) > 76
    Boolean hits_77 = length(hit_files) > 77
    Boolean hits_78 = length(hit_files) > 78
    Boolean hits_79 = length(hit_files) > 79
    Boolean hits_80 = length(hit_files) > 80
    Boolean hits_81 = length(hit_files) > 81
    Boolean hits_82 = length(hit_files) > 82
    Boolean hits_83 = length(hit_files) > 83
    Boolean hits_84 = length(hit_files) > 84
    Boolean hits_85 = length(hit_files) > 85
    Boolean hits_86 = length(hit_files) > 86
    Boolean hits_87 = length(hit_files) > 87
    Boolean hits_88 = length(hit_files) > 88
    Boolean hits_89 = length(hit_files) > 89
    Boolean hits_90 = length(hit_files) > 90
    Boolean hits_91 = length(hit_files) > 91
    Boolean hits_92 = length(hit_files) > 92
    Boolean hits_93 = length(hit_files) > 93
    Boolean hits_94 = length(hit_files) > 94
    Boolean hits_95 = length(hit_files) > 95
    Boolean hits_96 = length(hit_files) > 96
    Boolean hits_97 = length(hit_files) > 97
    Boolean hits_98 = length(hit_files) > 98
    Boolean hits_99 = length(hit_files) > 99
    Boolean hits_100 = length(hit_files) > 100
    Boolean hits_101 = length(hit_files) > 101
    Boolean hits_102 = length(hit_files) > 102
    Boolean hits_103 = length(hit_files) > 103
    Boolean hits_104 = length(hit_files) > 104
    Boolean hits_105 = length(hit_files) > 105
    Boolean hits_106 = length(hit_files) > 106
    Boolean hits_107 = length(hit_files) > 107
    Boolean hits_108 = length(hit_files) > 108
    Boolean hits_109 = length(hit_files) > 109
    Boolean hits_110 = length(hit_files) > 110
    Boolean hits_111 = length(hit_files) > 111
    Boolean hits_112 = length(hit_files) > 112
    Boolean hits_113 = length(hit_files) > 113
    Boolean hits_114 = length(hit_files) > 114
    Boolean hits_115 = length(hit_files) > 115
    Boolean hits_116 = length(hit_files) > 116
    Boolean hits_117 = length(hit_files) > 117
    Boolean hits_118 = length(hit_files) > 118
    Boolean hits_119 = length(hit_files) > 119
    Boolean hits_120 = length(hit_files) > 120
    Boolean hits_121 = length(hit_files) > 121
    Boolean hits_122 = length(hit_files) > 122
    Boolean hits_123 = length(hit_files) > 123
    Boolean hits_124 = length(hit_files) > 124
    Boolean hits_125 = length(hit_files) > 125
    Boolean hits_126 = length(hit_files) > 126
    Boolean hits_127 = length(hit_files) > 127
    Boolean changes_0 = length(change_files) > 0
    Boolean changes_1 = length(change_files) > 1
    Boolean changes_2 = length(change_files) > 2
    Boolean changes_3 = length(change_files) > 3
    Boolean changes_4 = length(change_files) > 4
    Boolean changes_5 = length(change_files) > 5
    Boolean changes_6 = length(change_files) > 6
    Boolean changes_7 = length(change_files) > 7
    Boolean changes_8 = length(change_files) > 8
    Boolean changes_9 = length(change_files) > 9
    Boolean changes_10 = length(change_files) > 10
    Boolean changes_11 = length(change_files) > 11
    Boolean changes_12 = length(change_files) > 12
    Boolean changes_13 = length(change_files) > 13
    Boolean changes_14 = length(change_files) > 14
    Boolean changes_15 = length(change_files) > 15
    Boolean changes_16 = length(change_files) > 16
    Boolean changes_17 = length(change_files) > 17
    Boolean changes_18 = length(change_files) > 18
    Boolean changes_19 = length(change_files) > 19
    Boolean changes_20 = length(change_files) > 20
    Boolean changes_21 = length(change_files) > 21
    Boolean changes_22 = length(change_files) > 22
    Boolean changes_23 = length(change_files) > 23
    Boolean changes_24 = length(change_files) > 24
    Boolean changes_25 = length(change_files) > 25
    Boolean changes_26 = length(change_files) > 26
    Boolean changes_27 = length(change_files) > 27
    Boolean changes_28 = length(change_files) > 28
    Boolean changes_29 = length(change_files) > 29
    Boolean changes_30 = length(change_files) > 30
    Boolean changes_31 = length(change_files) > 31
    Boolean changes_32 = length(change_files) > 32
    Boolean changes_33 = length(change_files) > 33
    Boolean changes_34 = length(change_files) > 34
    Boolean changes_35 = length(change_files) > 35
    Boolean changes_36 = length(change_files) > 36
    Boolean changes_37 = length(change_files) > 37
    Boolean changes_38 = length(change_files) > 38
    Boolean changes_39 = length(change_files) > 39
    Boolean changes_40 = length(change_files) > 40
    Boolean changes_41 = length(change_files) > 41
    Boolean changes_42 = length(change_files) > 42
    Boolean changes_43 = length(change_files) > 43
    Boolean changes_44 = length(change_files) > 44
    Boolean changes_45 = length(change_files) > 45
    Boolean changes_46 = length(change_files) > 46
    Boolean changes_47 = length(change_files) > 47
    Boolean changes_48 = length(change_files) > 48
    Boolean changes_49 = length(change_files) > 49
    Boolean changes_50 = length(change_files) > 50
    Boolean changes_51 = length(change_files) > 51
    Boolean changes_52 = length(change_files) > 52
    Boolean changes_53 = length(change_files) > 53
    Boolean changes_54 = length(change_files) > 54
    Boolean changes_55 = length(change_files) > 55
    Boolean changes_56 = length(change_files) > 56
    Boolean changes_57 = length(change_files) > 57
    Boolean changes_58 = length(change_files) > 58
    Boolean changes_59 = length(change_files) > 59
    Boolean changes_60 = length(change_files) > 60
    Boolean changes_61 = length(change_files) > 61
    Boolean changes_62 = length(change_files) > 62
    Boolean changes_63 = length(change_files) > 63
    Boolean changes_64 = length(change_files) > 64
    Boolean changes_65 = length(change_files) > 65
    Boolean changes_66 = length(change_files) > 66
    Boolean changes_67 = length(change_files) > 67
    Boolean changes_68 = length(change_files) > 68
    Boolean changes_69 = length(change_files) > 69
    Boolean changes_70 = length(change_files) > 70
    Boolean changes_71 = length(change_files) > 71
    Boolean changes_72 = length(change_files) > 72
    Boolean changes_73 = length(change_files) > 73
    Boolean changes_74 = length(change_files) > 74
    Boolean changes_75 = length(change_files) > 75
    Boolean changes_76 = length(change_files) > 76
    Boolean changes_77 = length(change_files) > 77
    Boolean changes_78 = length(change_files) > 78
    Boolean changes_79 = length(change_files) > 79
    Boolean changes_80 = length(change_files) > 80
    Boolean changes_81 = length(change_files) > 81
    Boolean changes_82 = length(change_files) > 82
    Boolean changes_83 = length(change_files) > 83
    Boolean changes_84 = length(change_files) > 84
    Boolean changes_85 = length(change_files) > 85
    Boolean changes_86 = length(change_files) > 86
    Boolean changes_87 = length(change_files) > 87
    Boolean changes_88 = length(change_files) > 88
    Boolean changes_89 = length(change_files) > 89
    Boolean changes_90 = length(change_files) > 90
    Boolean changes_91 = length(change_files) > 91
    Boolean changes_92 = length(change_files) > 92
    Boolean changes_93 = length(change_files) > 93
    Boolean changes_94 = length(change_files) > 94
    Boolean changes_95 = length(change_files) > 95
    Boolean changes_96 = length(change_files) > 96
    Boolean changes_97 = length(change_files) > 97
    Boolean changes_98 = length(change_files) > 98
    Boolean changes_99 = length(change_files) > 99
    Boolean changes_100 = length(change_files) > 100
    Boolean changes_101 = length(change_files) > 101
    Boolean changes_102 = length(change_files) > 102
    Boolean changes_103 = length(change_files) > 103
    Boolean changes_104 = length(change_files) > 104
    Boolean changes_105 = length(change_files) > 105
    Boolean changes_106 = length(change_files) > 106
    Boolean changes_107 = length(change_files) > 107
    Boolean changes_108 = length(change_files) > 108
    Boolean changes_109 = length(change_files) > 109
    Boolean changes_110 = length(change_files) > 110
    Boolean changes_111 = length(change_files) > 111
    Boolean changes_112 = length(change_files) > 112
    Boolean changes_113 = length(change_files) > 113
    Boolean changes_114 = length(change_files) > 114
    Boolean changes_115 = length(change_files) > 115
    Boolean changes_116 = length(change_files) > 116
    Boolean changes_117 = length(change_files) > 117
    Boolean changes_118 = length(change_files) > 118
    Boolean changes_119 = length(change_files) > 119
    Boolean changes_120 = length(change_files) > 120
    Boolean changes_121 = length(change_files) > 121
    Boolean changes_122 = length(change_files) > 122
    Boolean changes_123 = length(change_files) > 123
    Boolean changes_124 = length(change_files) > 124
    Boolean changes_125 = length(change_files) > 125
    Boolean changes_126 = length(change_files) > 126
    Boolean changes_127 = length(change_files) > 127
    Boolean variants_0 = length(variant_files) > 0
    Boolean variants_1 = length(variant_files) > 1
    Boolean variants_2 = length(variant_files) > 2
    Boolean variants_3 = length(variant_files) > 3
    Boolean variants_4 = length(variant_files) > 4
    Boolean variants_5 = length(variant_files) > 5
    Boolean variants_6 = length(variant_files) > 6
    Boolean variants_7 = length(variant_files) > 7
    Boolean variants_8 = length(variant_files) > 8
    Boolean variants_9 = length(variant_files) > 9
    Boolean variants_10 = length(variant_files) > 10
    Boolean variants_11 = length(variant_files) > 11
    Boolean variants_12 = length(variant_files) > 12
    Boolean variants_13 = length(variant_files) > 13
    Boolean variants_14 = length(variant_files) > 14
    Boolean variants_15 = length(variant_files) > 15
    Boolean variants_16 = length(variant_files) > 16
    Boolean variants_17 = length(variant_files) > 17
    Boolean variants_18 = length(variant_files) > 18
    Boolean variants_19 = length(variant_files) > 19
    Boolean variants_20 = length(variant_files) > 20
    Boolean variants_21 = length(variant_files) > 21
    Boolean variants_22 = length(variant_files) > 22
    Boolean variants_23 = length(variant_files) > 23
    Boolean variants_24 = length(variant_files) > 24
    Boolean variants_25 = length(variant_files) > 25
    Boolean variants_26 = length(variant_files) > 26
    Boolean variants_27 = length(variant_files) > 27
    Boolean variants_28 = length(variant_files) > 28
    Boolean variants_29 = length(variant_files) > 29
    Boolean variants_30 = length(variant_files) > 30
    Boolean variants_31 = length(variant_files) > 31
    Boolean variants_32 = length(variant_files) > 32
    Boolean variants_33 = length(variant_files) > 33
    Boolean variants_34 = length(variant_files) > 34
    Boolean variants_35 = length(variant_files) > 35
    Boolean variants_36 = length(variant_files) > 36
    Boolean variants_37 = length(variant_files) > 37
    Boolean variants_38 = length(variant_files) > 38
    Boolean variants_39 = length(variant_files) > 39
    Boolean variants_40 = length(variant_files) > 40
    Boolean variants_41 = length(variant_files) > 41
    Boolean variants_42 = length(variant_files) > 42
    Boolean variants_43 = length(variant_files) > 43
    Boolean variants_44 = length(variant_files) > 44
    Boolean variants_45 = length(variant_files) > 45
    Boolean variants_46 = length(variant_files) > 46
    Boolean variants_47 = length(variant_files) > 47
    Boolean variants_48 = length(variant_files) > 48
    Boolean variants_49 = length(variant_files) > 49
    Boolean variants_50 = length(variant_files) > 50
    Boolean variants_51 = length(variant_files) > 51
    Boolean variants_52 = length(variant_files) > 52
    Boolean variants_53 = length(variant_files) > 53
    Boolean variants_54 = length(variant_files) > 54
    Boolean variants_55 = length(variant_files) > 55
    Boolean variants_56 = length(variant_files) > 56
    Boolean variants_57 = length(variant_files) > 57
    Boolean variants_58 = length(variant_files) > 58
    Boolean variants_59 = length(variant_files) > 59
    Boolean variants_60 = length(variant_files) > 60
    Boolean variants_61 = length(variant_files) > 61
    Boolean variants_62 = length(variant_files) > 62
    Boolean variants_63 = length(variant_files) > 63
    Boolean variants_64 = length(variant_files) > 64
    Boolean variants_65 = length(variant_files) > 65
    Boolean variants_66 = length(variant_files) > 66
    Boolean variants_67 = length(variant_files) > 67
    Boolean variants_68 = length(variant_files) > 68
    Boolean variants_69 = length(variant_files) > 69
    Boolean variants_70 = length(variant_files) > 70
    Boolean variants_71 = length(variant_files) > 71
    Boolean variants_72 = length(variant_files) > 72
    Boolean variants_73 = length(variant_files) > 73
    Boolean variants_74 = length(variant_files) > 74
    Boolean variants_75 = length(variant_files) > 75
    Boolean variants_76 = length(variant_files) > 76
    Boolean variants_77 = length(variant_files) > 77
    Boolean variants_78 = length(variant_files) > 78
    Boolean variants_79 = length(variant_files) > 79
    Boolean variants_80 = length(variant_files) > 80
    Boolean variants_81 = length(variant_files) > 81
    Boolean variants_82 = length(variant_files) > 82
    Boolean variants_83 = length(variant_files) > 83
    Boolean variants_84 = length(variant_files) > 84
    Boolean variants_85 = length(variant_files) > 85
    Boolean variants_86 = length(variant_files) > 86
    Boolean variants_87 = length(variant_files) > 87
    Boolean variants_88 = length(variant_files) > 88
    Boolean variants_89 = length(variant_files) > 89
    Boolean variants_90 = length(variant_files) > 90
    Boolean variants_91 = length(variant_files) > 91
    Boolean variants_92 = length(variant_files) > 92
    Boolean variants_93 = length(variant_files) > 93
    Boolean variants_94 = length(variant_files) > 94
    Boolean variants_95 = length(variant_files) > 95
    Boolean variants_96 = length(variant_files) > 96
    Boolean variants_97 = length(variant_files) > 97
    Boolean variants_98 = length(variant_files) > 98
    Boolean variants_99 = length(variant_files) > 99
    Boolean variants_100 = length(variant_files) > 100
    Boolean variants_101 = length(variant_files) > 101
    Boolean variants_102 = length(variant_files) > 102
    Boolean variants_103 = length(variant_files) > 103
    Boolean variants_104 = length(variant_files) > 104
    Boolean variants_105 = length(variant_files) > 105
    Boolean variants_106 = length(variant_files) > 106
    Boolean variants_107 = length(variant_files) > 107
    Boolean variants_108 = length(variant_files) > 108
    Boolean variants_109 = length(variant_files) > 109
    Boolean variants_110 = length(variant_files) > 110
    Boolean variants_111 = length(variant_files) > 111
    Boolean variants_112 = length(variant_files) > 112
    Boolean variants_113 = length(variant_files) > 113
    Boolean variants_114 = length(variant_files) > 114
    Boolean variants_115 = length(variant_files) > 115
    Boolean variants_116 = length(variant_files) > 116
    Boolean variants_117 = length(variant_files) > 117
    Boolean variants_118 = length(variant_files) > 118
    Boolean variants_119 = length(variant_files) > 119
    Boolean variants_120 = length(variant_files) > 120
    Boolean variants_121 = length(variant_files) > 121
    Boolean variants_122 = length(variant_files) > 122
    Boolean variants_123 = length(variant_files) > 123
    Boolean variants_124 = length(variant_files) > 124
    Boolean variants_125 = length(variant_files) > 125
    Boolean variants_126 = length(variant_files) > 126
    Boolean variants_127 = length(variant_files) > 127
    command <<<
        set -euo pipefail
        exec > >(tee mergefinemo.log) 2>&1
        echo '[Fi-NeMo] Start MergeFinemo'
        if (( ~{length(hit_files)} < 1 || ~{length(hit_files)} > 128 || ~{length(change_files)} != ~{length(hit_files)} || ~{length(variant_files)} != ~{length(hit_files)} )); then echo 'Merge requires 1-128 aligned shards; increase shard_size'; exit 1; fi
        python /opt/finemo_pipeline/scripts/write_local_files.py --output hits.list \
            ~{if hits_0 then "--file '" + sub(hit_files[0], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_1 then "--file '" + sub(hit_files[1], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_2 then "--file '" + sub(hit_files[2], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_3 then "--file '" + sub(hit_files[3], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_4 then "--file '" + sub(hit_files[4], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_5 then "--file '" + sub(hit_files[5], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_6 then "--file '" + sub(hit_files[6], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_7 then "--file '" + sub(hit_files[7], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_8 then "--file '" + sub(hit_files[8], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_9 then "--file '" + sub(hit_files[9], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_10 then "--file '" + sub(hit_files[10], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_11 then "--file '" + sub(hit_files[11], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_12 then "--file '" + sub(hit_files[12], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_13 then "--file '" + sub(hit_files[13], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_14 then "--file '" + sub(hit_files[14], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_15 then "--file '" + sub(hit_files[15], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_16 then "--file '" + sub(hit_files[16], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_17 then "--file '" + sub(hit_files[17], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_18 then "--file '" + sub(hit_files[18], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_19 then "--file '" + sub(hit_files[19], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_20 then "--file '" + sub(hit_files[20], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_21 then "--file '" + sub(hit_files[21], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_22 then "--file '" + sub(hit_files[22], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_23 then "--file '" + sub(hit_files[23], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_24 then "--file '" + sub(hit_files[24], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_25 then "--file '" + sub(hit_files[25], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_26 then "--file '" + sub(hit_files[26], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_27 then "--file '" + sub(hit_files[27], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_28 then "--file '" + sub(hit_files[28], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_29 then "--file '" + sub(hit_files[29], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_30 then "--file '" + sub(hit_files[30], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_31 then "--file '" + sub(hit_files[31], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_32 then "--file '" + sub(hit_files[32], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_33 then "--file '" + sub(hit_files[33], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_34 then "--file '" + sub(hit_files[34], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_35 then "--file '" + sub(hit_files[35], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_36 then "--file '" + sub(hit_files[36], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_37 then "--file '" + sub(hit_files[37], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_38 then "--file '" + sub(hit_files[38], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_39 then "--file '" + sub(hit_files[39], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_40 then "--file '" + sub(hit_files[40], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_41 then "--file '" + sub(hit_files[41], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_42 then "--file '" + sub(hit_files[42], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_43 then "--file '" + sub(hit_files[43], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_44 then "--file '" + sub(hit_files[44], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_45 then "--file '" + sub(hit_files[45], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_46 then "--file '" + sub(hit_files[46], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_47 then "--file '" + sub(hit_files[47], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_48 then "--file '" + sub(hit_files[48], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_49 then "--file '" + sub(hit_files[49], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_50 then "--file '" + sub(hit_files[50], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_51 then "--file '" + sub(hit_files[51], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_52 then "--file '" + sub(hit_files[52], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_53 then "--file '" + sub(hit_files[53], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_54 then "--file '" + sub(hit_files[54], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_55 then "--file '" + sub(hit_files[55], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_56 then "--file '" + sub(hit_files[56], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_57 then "--file '" + sub(hit_files[57], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_58 then "--file '" + sub(hit_files[58], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_59 then "--file '" + sub(hit_files[59], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_60 then "--file '" + sub(hit_files[60], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_61 then "--file '" + sub(hit_files[61], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_62 then "--file '" + sub(hit_files[62], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_63 then "--file '" + sub(hit_files[63], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_64 then "--file '" + sub(hit_files[64], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_65 then "--file '" + sub(hit_files[65], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_66 then "--file '" + sub(hit_files[66], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_67 then "--file '" + sub(hit_files[67], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_68 then "--file '" + sub(hit_files[68], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_69 then "--file '" + sub(hit_files[69], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_70 then "--file '" + sub(hit_files[70], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_71 then "--file '" + sub(hit_files[71], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_72 then "--file '" + sub(hit_files[72], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_73 then "--file '" + sub(hit_files[73], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_74 then "--file '" + sub(hit_files[74], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_75 then "--file '" + sub(hit_files[75], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_76 then "--file '" + sub(hit_files[76], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_77 then "--file '" + sub(hit_files[77], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_78 then "--file '" + sub(hit_files[78], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_79 then "--file '" + sub(hit_files[79], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_80 then "--file '" + sub(hit_files[80], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_81 then "--file '" + sub(hit_files[81], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_82 then "--file '" + sub(hit_files[82], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_83 then "--file '" + sub(hit_files[83], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_84 then "--file '" + sub(hit_files[84], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_85 then "--file '" + sub(hit_files[85], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_86 then "--file '" + sub(hit_files[86], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_87 then "--file '" + sub(hit_files[87], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_88 then "--file '" + sub(hit_files[88], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_89 then "--file '" + sub(hit_files[89], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_90 then "--file '" + sub(hit_files[90], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_91 then "--file '" + sub(hit_files[91], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_92 then "--file '" + sub(hit_files[92], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_93 then "--file '" + sub(hit_files[93], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_94 then "--file '" + sub(hit_files[94], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_95 then "--file '" + sub(hit_files[95], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_96 then "--file '" + sub(hit_files[96], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_97 then "--file '" + sub(hit_files[97], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_98 then "--file '" + sub(hit_files[98], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_99 then "--file '" + sub(hit_files[99], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_100 then "--file '" + sub(hit_files[100], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_101 then "--file '" + sub(hit_files[101], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_102 then "--file '" + sub(hit_files[102], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_103 then "--file '" + sub(hit_files[103], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_104 then "--file '" + sub(hit_files[104], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_105 then "--file '" + sub(hit_files[105], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_106 then "--file '" + sub(hit_files[106], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_107 then "--file '" + sub(hit_files[107], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_108 then "--file '" + sub(hit_files[108], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_109 then "--file '" + sub(hit_files[109], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_110 then "--file '" + sub(hit_files[110], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_111 then "--file '" + sub(hit_files[111], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_112 then "--file '" + sub(hit_files[112], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_113 then "--file '" + sub(hit_files[113], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_114 then "--file '" + sub(hit_files[114], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_115 then "--file '" + sub(hit_files[115], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_116 then "--file '" + sub(hit_files[116], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_117 then "--file '" + sub(hit_files[117], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_118 then "--file '" + sub(hit_files[118], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_119 then "--file '" + sub(hit_files[119], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_120 then "--file '" + sub(hit_files[120], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_121 then "--file '" + sub(hit_files[121], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_122 then "--file '" + sub(hit_files[122], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_123 then "--file '" + sub(hit_files[123], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_124 then "--file '" + sub(hit_files[124], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_125 then "--file '" + sub(hit_files[125], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_126 then "--file '" + sub(hit_files[126], "'", "'\"'\"'") + "'" else ""} \
            ~{if hits_127 then "--file '" + sub(hit_files[127], "'", "'\"'\"'") + "'" else ""}
        python /opt/finemo_pipeline/scripts/write_local_files.py --output changes.list \
            ~{if changes_0 then "--file '" + sub(change_files[0], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_1 then "--file '" + sub(change_files[1], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_2 then "--file '" + sub(change_files[2], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_3 then "--file '" + sub(change_files[3], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_4 then "--file '" + sub(change_files[4], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_5 then "--file '" + sub(change_files[5], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_6 then "--file '" + sub(change_files[6], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_7 then "--file '" + sub(change_files[7], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_8 then "--file '" + sub(change_files[8], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_9 then "--file '" + sub(change_files[9], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_10 then "--file '" + sub(change_files[10], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_11 then "--file '" + sub(change_files[11], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_12 then "--file '" + sub(change_files[12], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_13 then "--file '" + sub(change_files[13], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_14 then "--file '" + sub(change_files[14], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_15 then "--file '" + sub(change_files[15], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_16 then "--file '" + sub(change_files[16], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_17 then "--file '" + sub(change_files[17], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_18 then "--file '" + sub(change_files[18], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_19 then "--file '" + sub(change_files[19], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_20 then "--file '" + sub(change_files[20], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_21 then "--file '" + sub(change_files[21], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_22 then "--file '" + sub(change_files[22], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_23 then "--file '" + sub(change_files[23], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_24 then "--file '" + sub(change_files[24], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_25 then "--file '" + sub(change_files[25], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_26 then "--file '" + sub(change_files[26], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_27 then "--file '" + sub(change_files[27], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_28 then "--file '" + sub(change_files[28], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_29 then "--file '" + sub(change_files[29], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_30 then "--file '" + sub(change_files[30], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_31 then "--file '" + sub(change_files[31], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_32 then "--file '" + sub(change_files[32], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_33 then "--file '" + sub(change_files[33], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_34 then "--file '" + sub(change_files[34], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_35 then "--file '" + sub(change_files[35], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_36 then "--file '" + sub(change_files[36], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_37 then "--file '" + sub(change_files[37], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_38 then "--file '" + sub(change_files[38], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_39 then "--file '" + sub(change_files[39], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_40 then "--file '" + sub(change_files[40], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_41 then "--file '" + sub(change_files[41], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_42 then "--file '" + sub(change_files[42], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_43 then "--file '" + sub(change_files[43], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_44 then "--file '" + sub(change_files[44], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_45 then "--file '" + sub(change_files[45], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_46 then "--file '" + sub(change_files[46], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_47 then "--file '" + sub(change_files[47], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_48 then "--file '" + sub(change_files[48], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_49 then "--file '" + sub(change_files[49], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_50 then "--file '" + sub(change_files[50], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_51 then "--file '" + sub(change_files[51], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_52 then "--file '" + sub(change_files[52], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_53 then "--file '" + sub(change_files[53], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_54 then "--file '" + sub(change_files[54], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_55 then "--file '" + sub(change_files[55], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_56 then "--file '" + sub(change_files[56], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_57 then "--file '" + sub(change_files[57], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_58 then "--file '" + sub(change_files[58], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_59 then "--file '" + sub(change_files[59], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_60 then "--file '" + sub(change_files[60], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_61 then "--file '" + sub(change_files[61], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_62 then "--file '" + sub(change_files[62], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_63 then "--file '" + sub(change_files[63], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_64 then "--file '" + sub(change_files[64], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_65 then "--file '" + sub(change_files[65], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_66 then "--file '" + sub(change_files[66], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_67 then "--file '" + sub(change_files[67], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_68 then "--file '" + sub(change_files[68], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_69 then "--file '" + sub(change_files[69], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_70 then "--file '" + sub(change_files[70], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_71 then "--file '" + sub(change_files[71], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_72 then "--file '" + sub(change_files[72], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_73 then "--file '" + sub(change_files[73], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_74 then "--file '" + sub(change_files[74], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_75 then "--file '" + sub(change_files[75], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_76 then "--file '" + sub(change_files[76], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_77 then "--file '" + sub(change_files[77], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_78 then "--file '" + sub(change_files[78], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_79 then "--file '" + sub(change_files[79], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_80 then "--file '" + sub(change_files[80], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_81 then "--file '" + sub(change_files[81], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_82 then "--file '" + sub(change_files[82], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_83 then "--file '" + sub(change_files[83], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_84 then "--file '" + sub(change_files[84], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_85 then "--file '" + sub(change_files[85], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_86 then "--file '" + sub(change_files[86], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_87 then "--file '" + sub(change_files[87], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_88 then "--file '" + sub(change_files[88], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_89 then "--file '" + sub(change_files[89], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_90 then "--file '" + sub(change_files[90], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_91 then "--file '" + sub(change_files[91], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_92 then "--file '" + sub(change_files[92], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_93 then "--file '" + sub(change_files[93], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_94 then "--file '" + sub(change_files[94], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_95 then "--file '" + sub(change_files[95], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_96 then "--file '" + sub(change_files[96], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_97 then "--file '" + sub(change_files[97], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_98 then "--file '" + sub(change_files[98], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_99 then "--file '" + sub(change_files[99], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_100 then "--file '" + sub(change_files[100], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_101 then "--file '" + sub(change_files[101], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_102 then "--file '" + sub(change_files[102], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_103 then "--file '" + sub(change_files[103], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_104 then "--file '" + sub(change_files[104], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_105 then "--file '" + sub(change_files[105], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_106 then "--file '" + sub(change_files[106], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_107 then "--file '" + sub(change_files[107], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_108 then "--file '" + sub(change_files[108], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_109 then "--file '" + sub(change_files[109], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_110 then "--file '" + sub(change_files[110], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_111 then "--file '" + sub(change_files[111], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_112 then "--file '" + sub(change_files[112], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_113 then "--file '" + sub(change_files[113], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_114 then "--file '" + sub(change_files[114], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_115 then "--file '" + sub(change_files[115], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_116 then "--file '" + sub(change_files[116], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_117 then "--file '" + sub(change_files[117], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_118 then "--file '" + sub(change_files[118], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_119 then "--file '" + sub(change_files[119], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_120 then "--file '" + sub(change_files[120], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_121 then "--file '" + sub(change_files[121], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_122 then "--file '" + sub(change_files[122], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_123 then "--file '" + sub(change_files[123], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_124 then "--file '" + sub(change_files[124], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_125 then "--file '" + sub(change_files[125], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_126 then "--file '" + sub(change_files[126], "'", "'\"'\"'") + "'" else ""} \
            ~{if changes_127 then "--file '" + sub(change_files[127], "'", "'\"'\"'") + "'" else ""}
        python /opt/finemo_pipeline/scripts/write_local_files.py --output variants.list \
            ~{if variants_0 then "--file '" + sub(variant_files[0], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_1 then "--file '" + sub(variant_files[1], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_2 then "--file '" + sub(variant_files[2], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_3 then "--file '" + sub(variant_files[3], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_4 then "--file '" + sub(variant_files[4], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_5 then "--file '" + sub(variant_files[5], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_6 then "--file '" + sub(variant_files[6], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_7 then "--file '" + sub(variant_files[7], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_8 then "--file '" + sub(variant_files[8], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_9 then "--file '" + sub(variant_files[9], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_10 then "--file '" + sub(variant_files[10], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_11 then "--file '" + sub(variant_files[11], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_12 then "--file '" + sub(variant_files[12], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_13 then "--file '" + sub(variant_files[13], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_14 then "--file '" + sub(variant_files[14], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_15 then "--file '" + sub(variant_files[15], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_16 then "--file '" + sub(variant_files[16], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_17 then "--file '" + sub(variant_files[17], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_18 then "--file '" + sub(variant_files[18], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_19 then "--file '" + sub(variant_files[19], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_20 then "--file '" + sub(variant_files[20], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_21 then "--file '" + sub(variant_files[21], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_22 then "--file '" + sub(variant_files[22], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_23 then "--file '" + sub(variant_files[23], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_24 then "--file '" + sub(variant_files[24], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_25 then "--file '" + sub(variant_files[25], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_26 then "--file '" + sub(variant_files[26], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_27 then "--file '" + sub(variant_files[27], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_28 then "--file '" + sub(variant_files[28], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_29 then "--file '" + sub(variant_files[29], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_30 then "--file '" + sub(variant_files[30], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_31 then "--file '" + sub(variant_files[31], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_32 then "--file '" + sub(variant_files[32], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_33 then "--file '" + sub(variant_files[33], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_34 then "--file '" + sub(variant_files[34], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_35 then "--file '" + sub(variant_files[35], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_36 then "--file '" + sub(variant_files[36], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_37 then "--file '" + sub(variant_files[37], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_38 then "--file '" + sub(variant_files[38], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_39 then "--file '" + sub(variant_files[39], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_40 then "--file '" + sub(variant_files[40], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_41 then "--file '" + sub(variant_files[41], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_42 then "--file '" + sub(variant_files[42], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_43 then "--file '" + sub(variant_files[43], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_44 then "--file '" + sub(variant_files[44], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_45 then "--file '" + sub(variant_files[45], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_46 then "--file '" + sub(variant_files[46], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_47 then "--file '" + sub(variant_files[47], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_48 then "--file '" + sub(variant_files[48], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_49 then "--file '" + sub(variant_files[49], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_50 then "--file '" + sub(variant_files[50], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_51 then "--file '" + sub(variant_files[51], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_52 then "--file '" + sub(variant_files[52], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_53 then "--file '" + sub(variant_files[53], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_54 then "--file '" + sub(variant_files[54], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_55 then "--file '" + sub(variant_files[55], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_56 then "--file '" + sub(variant_files[56], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_57 then "--file '" + sub(variant_files[57], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_58 then "--file '" + sub(variant_files[58], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_59 then "--file '" + sub(variant_files[59], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_60 then "--file '" + sub(variant_files[60], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_61 then "--file '" + sub(variant_files[61], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_62 then "--file '" + sub(variant_files[62], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_63 then "--file '" + sub(variant_files[63], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_64 then "--file '" + sub(variant_files[64], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_65 then "--file '" + sub(variant_files[65], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_66 then "--file '" + sub(variant_files[66], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_67 then "--file '" + sub(variant_files[67], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_68 then "--file '" + sub(variant_files[68], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_69 then "--file '" + sub(variant_files[69], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_70 then "--file '" + sub(variant_files[70], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_71 then "--file '" + sub(variant_files[71], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_72 then "--file '" + sub(variant_files[72], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_73 then "--file '" + sub(variant_files[73], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_74 then "--file '" + sub(variant_files[74], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_75 then "--file '" + sub(variant_files[75], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_76 then "--file '" + sub(variant_files[76], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_77 then "--file '" + sub(variant_files[77], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_78 then "--file '" + sub(variant_files[78], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_79 then "--file '" + sub(variant_files[79], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_80 then "--file '" + sub(variant_files[80], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_81 then "--file '" + sub(variant_files[81], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_82 then "--file '" + sub(variant_files[82], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_83 then "--file '" + sub(variant_files[83], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_84 then "--file '" + sub(variant_files[84], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_85 then "--file '" + sub(variant_files[85], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_86 then "--file '" + sub(variant_files[86], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_87 then "--file '" + sub(variant_files[87], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_88 then "--file '" + sub(variant_files[88], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_89 then "--file '" + sub(variant_files[89], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_90 then "--file '" + sub(variant_files[90], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_91 then "--file '" + sub(variant_files[91], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_92 then "--file '" + sub(variant_files[92], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_93 then "--file '" + sub(variant_files[93], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_94 then "--file '" + sub(variant_files[94], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_95 then "--file '" + sub(variant_files[95], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_96 then "--file '" + sub(variant_files[96], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_97 then "--file '" + sub(variant_files[97], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_98 then "--file '" + sub(variant_files[98], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_99 then "--file '" + sub(variant_files[99], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_100 then "--file '" + sub(variant_files[100], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_101 then "--file '" + sub(variant_files[101], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_102 then "--file '" + sub(variant_files[102], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_103 then "--file '" + sub(variant_files[103], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_104 then "--file '" + sub(variant_files[104], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_105 then "--file '" + sub(variant_files[105], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_106 then "--file '" + sub(variant_files[106], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_107 then "--file '" + sub(variant_files[107], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_108 then "--file '" + sub(variant_files[108], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_109 then "--file '" + sub(variant_files[109], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_110 then "--file '" + sub(variant_files[110], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_111 then "--file '" + sub(variant_files[111], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_112 then "--file '" + sub(variant_files[112], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_113 then "--file '" + sub(variant_files[113], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_114 then "--file '" + sub(variant_files[114], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_115 then "--file '" + sub(variant_files[115], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_116 then "--file '" + sub(variant_files[116], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_117 then "--file '" + sub(variant_files[117], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_118 then "--file '" + sub(variant_files[118], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_119 then "--file '" + sub(variant_files[119], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_120 then "--file '" + sub(variant_files[120], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_121 then "--file '" + sub(variant_files[121], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_122 then "--file '" + sub(variant_files[122], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_123 then "--file '" + sub(variant_files[123], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_124 then "--file '" + sub(variant_files[124], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_125 then "--file '" + sub(variant_files[125], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_126 then "--file '" + sub(variant_files[126], "'", "'\"'\"'") + "'" else ""} \
            ~{if variants_127 then "--file '" + sub(variant_files[127], "'", "'\"'\"'") + "'" else ""}
        python /opt/finemo_pipeline/scripts/merge_finemo.py \
            --hits-list hits.list \
            --changes-list changes.list \
            --variants-list variants.list \
            --output-dir merged
        echo '[Fi-NeMo] Complete MergeFinemo'
    >>>
    output {
        File hits = 'merged/annotated_hits.tsv'
        File changes = 'merged/motif_changes.tsv'
        File variants = 'merged/variant_summary.tsv'
        File metadata = 'merged/metadata.json'
        File log = 'mergefinemo.log'
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
