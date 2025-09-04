//
nextflow.enable.dsl = 2

process COLLATE_STATISTICS {
    tag "concatenate_stats"
    publishDir "${params.output_dir}/Collated_summary", mode: params.publish_mode
    
    input:
    path stats_files

    output:
    path "combined_statistics.tsv", emit: combined_stats

    script:
    """
    python ${projectDir}/bin/Collate_Stats.py \
        "combined_statistics.tsv" \
        $stats_files
    """
}