/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Import Required Workflows and Processes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BAM_TO_MOD_TABLE } from "../bin/process.nf"
include { METHYLATION_ANALYSIS } from "../bin/process.nf"
include { ISOLATION_BY_READ_ID } from "../bin/process.nf"
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

        isolate_telomeric_reads = ISOLATION_BY_READ_ID(telomeric_reads_with_stats, putative_reads)
        BAM_TO_MOD_TABLE(isolate_telomeric_reads.non_processed_reads)
        methylation_results = METHYLATION_ANALYSIS(BAM_TO_MOD_TABLE.out.mod_table_out)
        //methylation_analysis = METHYLATION_ANALYSIS(modkit_conversion)
        //METHYLATION_ANALYSIS.view()
    
    emit:
       methylation_results.modification_stats

}