#!/usr/bin/env nextflow

// Declare syntax version
nextflow.enable.dsl=2

// WORKFLOWS
include { MERGE_BREAKENDS }         from './workflows/local/Merge_Breakends/'
include { BEDTOOLS_INTERSECT }      from './workflows/local/Bedtools_Intersect/' 
include { COUNT_INTERSECTIONS }     from './workflows/local/Count_Intersections/'
//include { COLLECT_STATISTICS }      from './workflows/local/Collect_Statistics/'
include { COLLATE_STATISTICS }      from './workflows/local/Collate_Statistics/'


// SUBWORKFLOWS
// Not included here


workflow {
    
    // Create channels 
    reference_ch = Channel.fromPath(params.reference_asisi, checkIfExists: true)
    samples_ch = Channel.fromPath(params.sample_beds, checkIfExists: true)
    
    // Merge adjacent reads
    merged = MERGE_BREAKENDS(samples_ch)
    
    // Intersect with reference
    // map reference bed to each sample
    reference_for_each_sample = merged.merged_bed.map { file -> params.reference_asisi }
    // intersections
    intersections = BEDTOOLS_INTERSECT(merged.merged_bed, reference_for_each_sample)
    // extract outputs
    overlapping_regions = intersections.intersected
    
    // ensure 1:1 pairing
    counts_input = overlapping_regions.combine(reference_ch)
        .filter { it.size() == 2 } // Ensure pairs
        .map { overlapping_bed, reference_bed -> tuple(overlapping_bed, reference_bed) }

    // print
    //counts_input.view() { overlapping_bed, reference_bed -> "test ${overlapping_bed.name}" }

    counts = COUNT_INTERSECTIONS(counts_input)
    // Collect all statistics files
    all_stats = counts.summary.collect()

     // Concatenate into single file
    combined = COLLATE_STATISTICS(all_stats)
    // View result
    combined.combined_stats.view()


}