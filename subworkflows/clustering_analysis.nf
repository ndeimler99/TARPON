/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Import Required Workflows and Processes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { telogator_clustering } from "../bin/process.nf"
//include { plot_clustering_results } from "../bin/process.nf"

workflow clustering_pipeline {

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

        clustering_results = telogator_clustering(telomeric_reads_with_stats)
        
    emit:
        clustering_results.clustering_stats

}