// Process to count intersections
nextflow.enable.dsl = 2

process COUNT_INTERSECTIONS {
    tag "${intersected_bed.baseName}"
    publishDir "${params.output_dir}/intersection_counts", mode: params.publish_mode, pattern: "*.count.txt"
    publishDir "${params.output_dir}/intersection_summaries", mode: params.publish_mode, pattern: "*.summary.txt"

    // input tuble of intersected files and Asisi reference
    input:
    tuple path(intersected_bed), path(reference_bed)

    output:
    path "${intersected_bed.baseName}.count.txt", emit: count
    path "${intersected_bed.baseName}.summary.txt", emit: summary
    

    script:
    """
    python ${projectDir}/bin/Count_intersections.py \
        "$intersected_bed" \
        "$reference_bed" \
        "${intersected_bed.baseName}" \
        ${params.binsize}

    """
}