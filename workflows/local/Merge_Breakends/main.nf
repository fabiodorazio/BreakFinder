nextflow.enable.dsl = 2

process MERGE_BREAKENDS {
    // tagging process for logs
    tag "${bed_file.baseName}"
    publishDir "${params.output_dir}/merged_breaks", mode: params.publish_mode, pattern: "*.bed"

    input:
        path bed_file

    // returns merged bed files
    output:
        path "${bed_file.baseName}_merged.bed", emit: merged_bed

    script:
        """
        set -euo pipefail
        input_file="$bed_file"
        baseNoExt="\${input_file%.*}"

        # filter on quality (col5)
        # merge adjacent break sites (distance = 0) and count occurrences
        # column 4 in beds has Break IDs
        awk -v min=${params.min_score} '(\$5 >= min)' "\$input_file" \
        | bedtools sort -i - \
        | bedtools merge -i - -d "${params.merge_overlaps}" -c 4 -o count \
        > "\${baseNoExt}_merged.bed"    
        """
}