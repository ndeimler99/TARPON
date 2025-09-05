/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Import Required Workflows and Processes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BAM_TO_MOD_TABLE } from "../bin/process.nf"
include { METHYLATION_ANALYSIS } from "../bin/process.nf"
//include { plot_clustering_results } from "../bin/process.nf"

workflow methylation_detection {

    take:
        input_ch
        putative_reads
        non_telomeric
        sample_file
        reversed_retained
        reversed_removed 
        adaptor_retained_all_samples 
        adaptor_removed_all_samples
        demuxed_reads
        subtelomeric_retained
        subtelomeric_removed
        telomeric_reads_with_stats
        telomeric_reads_retained
        telomeric_reads_no_start
        telomeric_reads_under_threshold

    main:

        methylation_analysis = METHYLATION_ANALYSIS(modkit_conversion)
        
    emit:
        methylation_analysis.modification_stats

}