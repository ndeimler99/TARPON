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
include { SUMMARY_STATS } from "./bin/process.nf"
include { RESTRICTION_DIGEST_ANALYSIS } from "./bin/process.nf"
include { GENERATE_HTML_REPORT } from "./bin/process.nf"

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
    //reversed_ch.reversed_reads.view()
    //reversed_ch.removed_reads.view()
    //reversed_ch.strand_specific.view()
    
    // isolate reads with adaptor sequences and filter by subtelomere size
    adaptor_ch = IDENTIFY_TAGGING_ADAPTOR_AND_DEMUX(reversed_ch.reversed_reads)
    // adaptor_ch.demuxed_reads.view()
    

    adaptor_ch.demuxed_reads
        .map { it -> [it.baseName, it]}
        .set { demuxed_reads }
    // demuxed_reads.view()
        
    //telo start and length determination
    // analyze reads and create stats file containing read_id, strand, read_len, VRR_Start, VRR_length, Telo_length
    telo_stats = TELO_START_IDENTIFICATION(demuxed_reads)
    //telo_stats.retained_reads.view()
    //telo_stats.removed_reads.view()
    //telo_stats.telomeric_sequences.view()
    
    //if individual read, create plots
    if (params.indiv_read_plots) {
        INDIVIDUAL_READ_PLOTS(telo_stats.telomeric_sequences)
    }

    GENERATE_PLOTS(telo_stats.telomeric_sequences)

    //     stats_done = summary_stats(telo_stats.telomeric_stats, file(params.outdir))

    if (params.detailed_stats) {
        telo_stats = GENERATE_DETAILED_PLOTS(telo_stats.telomeric_sequences)
    }

    //     if (params.restriction_digest_analysis != ""){
    //         restriction_digest_analysis(telo_stats.telomeric_sequences)
    //     }


    //     generate_html_report(file(params.outdir), stats_done[1])
}


workflow.onComplete {
    println "Analysis Complete at: $workflow.complete"
    println "Execution Status: ${ workflow.success ? 'OK' : 'failed' }"
    println "Open the Following Report in your Browser ${ params.outdir }/report.html"

    if (workflow.success){
        if (params.remove_wd) {
            "rm -rf ${baseDir}/work".execute()
        }
        
        if (params.remove_intermediate_fastq){
            "rm -rf ${params.outdir}/TELOMERIC".execute()
            "rm -rf ${params.outdir}/FILTERED_READS".execute()
        }
    }
}
