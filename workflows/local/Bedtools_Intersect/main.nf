// Process to perform bedtools intersect
nextflow.enable.dsl = 2

process BEDTOOLS_INTERSECT {
    // tagging process for logs
    tag "${merged_beds.baseName}"
    publishDir "${params.output_dir}/intersections", mode: params.publish_mode, pattern: "*.bed"


    // merged beds and reference AsiSi
    input:
        path merged_beds
        path reference_asisi

    // intersected beds and count of intersections
    output:
        path "${merged_beds.baseName}_intersected_${params.merge_dist}.bed", emit: intersected
        path "${merged_beds.baseName}_nonintersected_${params.merge_dist}.bed", emit: non_intersected

    script:
        def options = []
            if (params.write_a) options << "-wa"
            if (params.write_b) options << "-wb"
            if (params.write_o) options << "-wo"
            if (params.merge_dist > 0) options << "-f ${params.merge_dist}"
        // join multiple options
        def options_str = ""
            if (options) {
                options_str = options.join(" ")
            }

        """
        # Perform intersection
        set -euo pipefail

        bedtools intersect -a "$merged_beds" -b "$reference_asisi" $options_str > "${merged_beds.baseName}_intersected_${params.merge_dist}.bed"
        bedtools intersect -a "$merged_beds" -b "$reference_asisi" $options_str -v > "${merged_beds.baseName}_nonintersected_${params.merge_dist}.bed"
    
        """
}
