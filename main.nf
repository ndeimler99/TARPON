#!/usr/bin/env nextflow

/* The following pipeline is intended for research purposes only */
nextflow.enable.dsl=2

/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

println """\
    TArPON - Telomere Analysis Pipeline on Nanopore Sequencing Data
    ================================================
    v0.0.1
    """.stripIndent()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Import Required Workflows and Processes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// include { process_name } from "process_file"
include { PUTATIVE_ISOLATION } from "./bin/process.nf"
include { REVERSE_COMPLEMENTATION } from "./bin/process.nf"
include { IDENTIFY_TAGGING_ADAPTOR_AND_DEMUX } from "./bin/process.nf"
include { TELO_START_IDENTIFICATION } from "./bin/process.nf"
include { INDIVIDUAL_READ_PLOTS } from "./bin/process.nf"
include { GENERATE_PLOTS } from "./bin/process.nf"
include { GENERATE_DETAILED_PLOTS } from "./bin/process.nf"
include { SUMMARY_STATS_RUN as SUMMARY_STATS1 } from "./bin/process.nf"
include { SUMMARY_STATS_SAMPLE as SUMMARY_STATS2 } from "./bin/process.nf"
include { RESTRICTION_DIGEST_ANALYSIS } from "./bin/process.nf"
include { GENERATE_FINAL_REPORT } from "./bin/process.nf"
include { SEPERATE_STRANDS as SEPERATE_STRAND_RUN_FILTERING } from "./bin/process.nf"
include { SEPERATE_STRANDS as SEPERATE_STRAND_SAMPLE_FILTERING } from "./bin/process.nf"
include { SEPERATE_STRANDS as SEPERATE_STRAND_RUN_RETAINED } from "./bin/process.nf"
include { SEPERATE_STRANDS as SEPERATE_STRAND_SAMPLE_RETAINED } from "./bin/process.nf"
include { SUBTELO_FILTERING } from "./bin/process.nf"
// include { cleanUp } from "./bin/process.nf"
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    
    // ###### TO DO #######
    // check input_file is a directory or a file
    // if directory concat all files in directory ending in .fastq.gz
    // #############

    // check if input file exists
    // input_file = file(params.input_file)
    // if (!input_file.exists() ) {
    //     exit 1, "Input File Does Not Exist"
    // }

    // // check if output directory exists
    // outdir = file(params.outdir)
    // if (outdir.exists() && !params.overwrite_outdir) {
    //     exit 1, "Out Directory Already Exists, Please Provide New Out Directory Name or Allow Overwriting of Pre-existing directory"
    // }

    // check that all other parameters are valid...
    
    
    // ###### TO DO #######
    // check if filtered_telo is passed in
    // else - filter telomeric reads

    // reference_file for mapping has to be value channel (just leave it at params.reference)
    // #############
    
    // convert all bam files to fastq

    Channel.fromPath( params.input_file, checkIfExists:true)
        .map{ it -> [params.run_name, it] }
        .set{ input_ch }

    //input_ch.view()

    // putative identification of telomeric sequences to limit dataset size
    putative_ch = PUTATIVE_ISOLATION(input_ch)
    //putative_ch.putative_reads.view()
    //putative_ch.non_telomeric.view()

    // reverse complement C strands into G strands and modify header line to include strand information
    reversed_ch = REVERSE_COMPLEMENTATION(putative_ch.putative_reads)
    //reversed_ch.retained_reads.view()
    //reversed_ch.removed_reads.view()
    
    // isolate reads with adaptor sequences and filter by subtelomere size
    adaptor_ch = IDENTIFY_TAGGING_ADAPTOR_AND_DEMUX(reversed_ch.retained_reads)
    //println "Adaptor Reads"
    //adaptor_ch.retained_reads.view()
    //adaptor_ch.demuxed_reads.view()
    //adaptor_ch.filtered_reads.view()

    adaptor_ch.demuxed_reads
        .map { it -> [it.baseName, it]}
        .set { demuxed_reads }
    // println "Demuxed Reads"
    // demuxed_reads.view()
        
    // filter by subtelo length 
    subtelo_filtered_ch = SUBTELO_FILTERING(demuxed_reads)
    //subtelo_filtered_ch.retained_reads.view()

    //telo start and length determination
    // analyze reads and create stats file containing read_id, strand, read_len, VRR_Start, VRR_length, Telo_length
    telo_stats = TELO_START_IDENTIFICATION(subtelo_filtered_ch.retained_reads)

    telo_stats.removed_reads
        .multiMap { it ->
            no_start: tuple(it[0], it[1][1])
            low_percent: tuple(it[0], it[1][0])
        }
        .set { filtered_telo }
    
    //filtered_telo.no_start.view()
    //filtered_telo.low_percent.view()
    //telo_stats.removed_reads.view()
    //telo_stats.final_telomeric.view()


    if (params.strand_comparison){
        separate_run_filtered = SEPERATE_STRAND_RUN_FILTERING(adaptor_ch.filtered_reads)
        separate_run_retained = SEPERATE_STRAND_RUN_RETAINED(reversed_ch.retained_reads.mix(adaptor_ch.retained_reads))
        //filtered_run.stranded.view()
        //retained_run.stranded.view()

        filtered_sample = SEPERATE_STRAND_SAMPLE_FILTERING(subtelo_filtered_ch.filtered_reads.mix(filtered_telo.no_start, filtered_telo.low_percent))
        retained_sample = SEPERATE_STRAND_SAMPLE_RETAINED(demuxed_reads.mix(subtelo_filtered_ch.retained_reads, telo_stats.retained_reads))
        //filtered_sample.stranded.view()
        //retained_sample.stranded.view()
        
        //get run retained stats on all reads input, putative reads, putative strand specific, adaptor found, adaptor strand specific.
        run_retained = SUMMARY_STATS1(input_ch.mix(putative_ch.putative_reads, separate_run_retained.g_strand, separate_run_retained.c_strand, adaptor_ch.retained_reads).groupTuple(), \
                                        input_ch.mix(putative_ch.non_telomeric, reversed_ch.removed_reads, adaptor_ch.filtered_reads, separate_run_filtered.c_strand, separate_run_filtered.g_strand).groupTuple())
        // get sample retained stats on number of reads with adaptor, adaptor strand specific, subtelo pass, subtelo strand specific, telomeric, telomeric strand specfic
        sample_retained = SUMMARY_STATS2(demuxed_reads.mix(subtelo_filtered_ch.retained_reads, telo_stats.retained_reads, retained_sample.c_strand, retained_sample.g_strand).groupTuple(), \
                                        demuxed_reads.mix(subtelo_filtered_ch.filtered_reads, filtered_telo.no_start, filtered_telo.low_percent, filtered_sample.c_strand, filtered_sample.g_strand).groupTuple())
    }
    else {
        //get run retained stats on all reads input, putative reads, adaptor found
        run_retained = SUMMARY_STATS1(input_ch.mix(putative_ch.putative_reads, adaptor_ch.retained_reads).groupTuple(), "Retained_Reads.stats.txt")
        // get run filtered stats on non-telomeric reads, 20-80 filtered, no adaptor found
        run_filtered = SUMMARY_STATS2(input_ch.mix(putative_ch.non_telomeric, reversed_ch.removed_reads, adaptor_ch.filtered_reads).groupTuple(), "Filtered_Reads.stats.txt")
        // get sample retained stats on number of reads with adaptor, subtelo pass, telomeric
        sample_retained = SUMMARY_STATS3(demuxed_reads.mix(subtelo_filtered_ch.retained_reads, telo_stats.retained_reads).groupTuple(), "Retained_Reads.stats.txt")
        // get sample filtered on subtelo fail, no telo start, below threshold
        sample_filtered = SUMMARY_STATS4(demuxed_reads.mix(subtelo_filtered_ch.filtered_reads, filtered_telo.no_start, filtered_telo.low_percent).groupTuple(), "Filtered_Reads.stats.txt")
    }

    //if individual read, create plots
    if (params.indiv_read_plots) {
        INDIVIDUAL_READ_PLOTS(telo_stats.final_telomeric)
    }

    GENERATE_PLOTS(telo_stats.final_telomeric)

    if (params.detailed_stats) {
        telo_stats = GENERATE_DETAILED_PLOTS(telo_stats.final_telomeric)
    }

    if (params.restriction_digest_analysis != ""){
        RESTRICTION_DIGEST_ANALYSIS(telo_stats.final_telomeric)
    }


    //generate_html_report(file(params.outdir), stats_done[1])
    // GENERATE_FINAL_REPORT(run_retained.stats, run_filtered.stats,
    //     [subtelo_filtered_ch.retained_reads.countFastq(), telo_stats.retained_reads.countFastq()],
    //     [subtelo_filtered_ch.filtered_reads.countFastq(), filtered_telo.no_start.countFastq(), filtered_telo.low_percent.countFastq()]
    //     )
}


workflow.onComplete {
    println "Analysis Complete at: $workflow.complete"
    println "Execution Status: ${ workflow.success ? 'OK' : 'failed' }"
    println "Open the Following Report in your Browser ${ params.outdir }/report.html"

    if (workflow.success){
        if (params.remove_wd) {
            "rm -rf ${baseDir}/work".execute()
        }
       
    }
}
