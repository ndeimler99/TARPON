/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Import Required Workflows and Processes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { INDIVIDUAL_READ_PLOTS } from "../bin/process.nf"
include { BASIC_PLOTS } from "../bin/process.nf"
include { GENERATE_DETAILED_PLOTS } from "../bin/process.nf"
include { SUMMARY_STATS_RUN } from "../bin/process.nf"
include { SUMMARY_STATS_SAMPLE } from "../bin/process.nf"
include { RESTRICTION_DIGEST_ANALYSIS } from "../bin/process.nf"
include { SEPERATE_STRANDS as SEPERATE_STRAND_RUN_FILTERING } from "../bin/process.nf"
include { SEPERATE_STRANDS as SEPERATE_STRAND_SAMPLE_FILTERING } from "../bin/process.nf"
include { SEPERATE_STRANDS as SEPERATE_STRAND_RUN_RETAINED } from "../bin/process.nf"
include { SEPERATE_STRANDS as SEPERATE_STRAND_SAMPLE_RETAINED } from "../bin/process.nf"
include { COMBINE_BAM as COMBINED_RETAINED_BAM } from "../bin/process.nf"
include { COMBINE_BAM as COMBINED_FILTERED_BAM } from "../bin/process.nf"
include { FINAL_TELO_STATS } from "../bin/process.nf"
include { MUTANT_ANALYSIS } from "../bin/process.nf"
include { VARIANT_ANALYSIS } from "../bin/process.nf"
include { PLOT_TELO_GRAPHS } from "../bin/process.nf"


workflow telomere_analysis_pipeline {

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

        //subtelo_ch.view()
        if (params.mutant != "false") {
            mutant_analysis = MUTANT_ANALYSIS(telomeric_reads_with_stats)
        }
        variant_analysis = VARIANT_ANALYSIS(telomeric_reads_with_stats)

        PLOT_TELO_GRAPHS(variant_analysis.sequences)
    // //     // merge all stats relevant to retained telomeric sequences
    //     run_retained = COMBINED_RETAINED_BAM(subtelomeric_retained.multiMap { label, stats -> stats: stats }.collect(). map { it -> [ "subtelo_pass", it ]}.mix(
    //                                                 telomeric_reads_retained.multiMap { label, stats -> stats: stats }.collect().map {it -> [ "telo_retained", it]}))

    // //     // merge all stats relevant to filtered telomeric sequences
    //     run_removed = COMBINED_FILTERED_BAM(subtelomeric_removed.multiMap {label, stats -> stats:stats}.collect().map {it -> ["subtelo_fail", it]}.mix(
    //                                                 telomeric_reads_no_start.multiMap { label, stats -> stats: stats}.collect().map{it -> ["no_telo_start", it]},
    //                                                 telomeric_reads_under_threshold.multiMap { label, stats -> stats: stats }.collect().map{it -> ["below_telo_threshold", it]}))

    //     // separate reads into G and C strand and recalculate all statistics
    //     if (params.strand_comparison){
    //         separate_run_removed = SEPERATE_STRAND_RUN_FILTERING(adaptor_removed_all_samples.mix(run_removed.combined))
    //         separate_run_retained = SEPERATE_STRAND_RUN_RETAINED(reversed_retained.mix(adaptor_retained_all_samples, run_retained.combined))


    //         removed_sample = SEPERATE_STRAND_SAMPLE_FILTERING(subtelomeric_removed.mix(telomeric_reads_no_start, telomeric_reads_under_threshold))
    //         retained_sample = SEPERATE_STRAND_SAMPLE_RETAINED(demuxed_reads.mix(subtelomeric_retained, telomeric_reads_retained))

            
    //         //get run retained stats on all reads input, putative reads, putative strand specific, adaptor found, adaptor strand specific.
    //         run_stats = SUMMARY_STATS_RUN(input_ch.mix(putative_reads, reversed_retained, run_retained.combined, separate_run_retained.g_strand, separate_run_retained.c_strand, adaptor_retained_all_samples).groupTuple(), \
    //                                         input_ch.mix(non_telomeric, reversed_removed, run_removed.combined, adaptor_removed_all_samples, separate_run_removed.c_strand, separate_run_removed.g_strand).groupTuple())


    //         // get sample retained stats on number of reads with adaptor, adaptor strand specific, subtelo pass, subtelo strand specific, telomeric, telomeric strand specfic
    //         sample_stats = SUMMARY_STATS_SAMPLE(demuxed_reads.mix(subtelomeric_retained, telomeric_reads_retained, retained_sample.c_strand, retained_sample.g_strand).groupTuple(), \
    //                                         demuxed_reads.mix(subtelomeric_removed, telomeric_reads_no_start, telomeric_reads_under_threshold, removed_sample.c_strand, removed_sample.g_strand).groupTuple())

    //     }
    //     else {
    //         //get run retained stats on all reads input, putative reads, adaptor found
    //         run_stats = SUMMARY_STATS_RUN(input_ch.mix(putative_reads, reversed_retained, adaptor_retained_all_samples, run_retained.combined).groupTuple(), \
    //                                     input_ch.mix(non_telomeric, reversed_removed, adaptor_removed_all_samples, run_removed.combined).groupTuple())

    //         // // get sample retained stats on number of reads with adaptor, subtelo pass, telomeric
    //         sample_stats = SUMMARY_STATS_SAMPLE(demuxed_reads.mix(subtelomeric_retained, telomeric_reads_retained).groupTuple(), \
    //                                         demuxed_reads.mix(subtelomeric_removed, telomeric_reads_no_start, telomeric_reads_under_threshold).groupTuple())
        
    //     }

        //if individual reads is specified one plot per read showing telomeric percentage will be created
        if (params.indiv_read_plots) {
            INDIVIDUAL_READ_PLOTS(telomeric_reads_with_stats)
        }

        // this will generate basic statistical plots for the end user
        BASIC_PLOTS(telomeric_reads_with_stats)

        // if detailed stats is specified analyze telomeric sequences but include things such as percentage telomeric, etc.
        if (params.detailed_stats) {
            telo_stats = GENERATE_DETAILED_PLOTS(variant_analysis.statistics)
        }

        // if a restriction digest sequence is provided run a restrition digest analysis

        if (params.restriction_digest_analysis != ""){
            restriction_digest = RESTRICTION_DIGEST_ANALYSIS(telomeric_reads_with_stats)
        }
        else {
            Channel.from() \
                .map { it -> tuple(it[0], it[1]) } \
                .set { restriction_digest }
        }

        // variant_analysis.statistics.multiMap{ label, stats ->
        //         label: label
        //         stats: stats
        //     }.set { variant_analysis_stats }

        
       // variant_analysis.statistics.map { it[1] }.collect().view()

        FINAL_TELO_STATS(variant_analysis.statistics.map { it[1] }.collect())

        // combine all stats relevant to a sample into one tuple....
            // pass this new item into GENERATE_FINAL_REPORT and python script
        
        
        // run_stats.retained_stats - global statistics for the analysis of entire flow cell data
        // run_stats.removed_stats - global statistcis for the anlaysis of entire flow cell data
        // FINAL_TELO_STATS.out.vrr_stats - one line per sample regarding telomere statistic stats

        // sample_stats.retained_stats - one tuple per sample with sample .retained.stats.txt
        // sample_stats.removed_stats - one tuple per sample with sample .removed.stats.txt
        // variant_analysis.statistics - one tuple per sample with telomere related statistics
        // variant_analysis.repeat_distribution - one tuple per sample with one nucleotide point mutations of telomeric repeat

        // restriction_digest is either one tuple per sample or just an empty tuple - if it is empty it is not included below
        // sample_stats.retained_stats.mix(sample_stats.removed_stats, variant_analysis.statistics, variant_analysis.repeat_distribution, restriction_digest)
        //     .map { it[1] }
        //     .set { file_conglomerate }

        variant_analysis.statistics.mix(variant_analysis.repeat_distribution, restriction_digest)
            .map { it[1] }
            .set { file_conglomerate }
            
    emit:
        // flowcell_retained = run_stats.retained_stats
        // flowcell_removed = run_stats.removed_stats
        telomere_stats = FINAL_TELO_STATS.out.vrr_stats
        sample_specific_stats = file_conglomerate
}