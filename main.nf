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
include { SUMMARY_STATS_RUN } from "./bin/process.nf"
include { SUMMARY_STATS_SAMPLE } from "./bin/process.nf"
include { RESTRICTION_DIGEST_ANALYSIS } from "./bin/process.nf"
include { GENERATE_FINAL_REPORT } from "./bin/process.nf"
include { SEPERATE_STRANDS as SEPERATE_STRAND_RUN_FILTERING } from "./bin/process.nf"
include { SEPERATE_STRANDS as SEPERATE_STRAND_SAMPLE_FILTERING } from "./bin/process.nf"
include { SEPERATE_STRANDS as SEPERATE_STRAND_RUN_RETAINED } from "./bin/process.nf"
include { SEPERATE_STRANDS as SEPERATE_STRAND_SAMPLE_RETAINED } from "./bin/process.nf"
include { SUBTELO_FILTERING } from "./bin/process.nf"
include { getParams } from "./bin/process.nf"
include { getVersions } from "./bin/process.nf"
include { getManifest } from "./bin/process.nf"
include { COMBINE_FASTQ as COMBINED_RETAINED_FASTQ } from "./bin/process.nf"
include { COMBINE_FASTQ as COMBINED_FILTERED_FASTQ } from "./bin/process.nf"
include { COMBINED_INPUT as COMBINED_INPUT } from "./bin/process.nf"
include { validateParameters; paramsHelp; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'

// Print help message, supply typical command line usage for the pipeline


// this will check json --sample_file needs to match appropriate format using schema file - how do I do that?
//validateParameters()



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
WorkflowMain.initialise(workflow, params, log)

workflow {

    Pinguscript.ping_start(nextflow, workflow, params)

    if (params.adaptor_sequence == "" && params.sample_file == "") {
        println ("Adaptor Sequence and Sample File cannot both be empty")
        exit 0
    }

    try {
        if (params.sample_file != ""){
            file(params.sample_file, checkIfExists:true)
        }
    }
    catch (Exception e) {
        println("Error - Sample File not Found")
        exit 0
    }


    // ###### TO DO #######
    // check if provided input exists
    // check input_file is a directory or a file
    // if directory concat all files in directory ending in .fastq.gz
    // #############

    try {
        file(params.input_file, checkIfExists:true)
    }
    catch (Exception e) {
        println("Error - Input File or Directory Does not Exist")
        exit 0
    }


    // // check if output directory exists
    outdir = file(params.outdir)
    if (outdir.exists() && !params.overwrite_outdir) {
        exit 1, "Out Directory Already Exists, Please Provide New Out Directory Name or Allow Overwriting of Pre-existing directory"
    }

    // check that all other parameters are valid...
    
    
    // ###### TO DO #######
    // check if filtered_telo is passed in
    // else - filter telomeric reads

    // reference_file for mapping has to be value channel (just leave it at params.reference)
    // #############
    
    // convert all bam files to fastq

    run = params.run_name 

    if (file(params.input_file).isDirectory()) {
        input_ch = COMBINED_INPUT(Channel.fromPath ( "${params.input_file}/*q.gz" ).collect().map{ it -> ["input", it]})
    }

    else{
        Channel.fromPath( params.input_file, checkIfExists:true)
        .map{ it -> [ run , it] }
        .set{ input_ch }
    }    

    //input_ch.view()

    // putative identification of telomeric sequences to limit dataset size
    putative_ch = PUTATIVE_ISOLATION(input_ch)

    // reverse complement C strands into G strands and modify header line to include strand information
    reversed_ch = REVERSE_COMPLEMENTATION(putative_ch.putative_reads)
    //reversed_ch.retained_reads.view()
    //reversed_ch.removed_reads.view()
    
    // isolate reads with adaptor sequences and filter by subtelomere size
    adaptor_ch = IDENTIFY_TAGGING_ADAPTOR_AND_DEMUX(reversed_ch.retained_reads)

    adaptor_ch.demuxed_reads
        .map { it -> [it.baseName, it]}
        .set { demuxed_reads }

    // filter by subtelo length 
    subtelo_filtered_ch = SUBTELO_FILTERING(demuxed_reads)

    //telo start and length determination
    // analyze reads and create stats file containing read_id, strand, read_len, VRR_Start, VRR_length, Telo_length
    telo_stats = TELO_START_IDENTIFICATION(subtelo_filtered_ch.retained_reads)


    run_retained = COMBINED_RETAINED_FASTQ(subtelo_filtered_ch.retained_reads.multiMap { label, stats -> stats: stats }.collect(). map { it -> [ "subtelo_pass", it ]}.mix(
                                                telo_stats.retained_reads.multiMap { label, stats -> stats: stats }.collect().map {it -> [ "telo_retained", it]}))

    run_filtered = COMBINED_FILTERED_FASTQ(subtelo_filtered_ch.filtered_reads.multiMap {label, stats -> stats:stats}.collect().map {it -> ["subtelo_fail", it]}.mix(
                                                telo_stats.no_telo_start.multiMap { label, stats -> stats: stats}.collect().map{it -> ["no_telo_start", it]},
                                                telo_stats.below_telo_threshold.multiMap { label, stats -> stats: stats }.collect().map{it -> ["below_telo_threshold", it]}))

    if (params.strand_comparison){
        separate_run_filtered = SEPERATE_STRAND_RUN_FILTERING(adaptor_ch.filtered_reads.mix(run_filtered.combined))
        separate_run_retained = SEPERATE_STRAND_RUN_RETAINED(reversed_ch.retained_reads.mix(adaptor_ch.retained_reads, run_retained.combined))


        filtered_sample = SEPERATE_STRAND_SAMPLE_FILTERING(subtelo_filtered_ch.filtered_reads.mix(telo_stats.no_telo_start, telo_stats.below_telo_threshold))
        retained_sample = SEPERATE_STRAND_SAMPLE_RETAINED(demuxed_reads.mix(subtelo_filtered_ch.retained_reads, telo_stats.retained_reads))

        
        //get run retained stats on all reads input, putative reads, putative strand specific, adaptor found, adaptor strand specific.
        run_stats = SUMMARY_STATS_RUN(putative_ch.input_ch.mix(putative_ch.putative_reads, reversed_ch.retained_reads, run_retained.combined, separate_run_retained.g_strand, separate_run_retained.c_strand, adaptor_ch.retained_reads).groupTuple(), \
                                        putative_ch.input_ch.mix(putative_ch.non_telomeric, reversed_ch.removed_reads, run_filtered.combined, adaptor_ch.filtered_reads, separate_run_filtered.c_strand, separate_run_filtered.g_strand).groupTuple())
        // get sample retained stats on number of reads with adaptor, adaptor strand specific, subtelo pass, subtelo strand specific, telomeric, telomeric strand specfic
        sample_stats = SUMMARY_STATS_SAMPLE(demuxed_reads.mix(subtelo_filtered_ch.retained_reads, telo_stats.retained_reads, retained_sample.c_strand, retained_sample.g_strand).groupTuple(), \
                                        demuxed_reads.mix(subtelo_filtered_ch.filtered_reads, telo_stats.no_telo_start, telo_stats.below_telo_threshold, filtered_sample.c_strand, filtered_sample.g_strand).groupTuple())

    }
    else {
        //get run retained stats on all reads input, putative reads, adaptor found
        run_stats = SUMMARY_STATS_RUN(putative_ch.input_ch.mix(putative_ch.putative_reads, reversed_ch.retained_reads, adaptor_ch.retained_reads, run_retained.combined).groupTuple(), \
                                    putative_ch.input_ch.mix(putative_ch.non_telomeric, reversed_ch.removed_reads, adaptor_ch.filtered_reads, run_filtered.combined).groupTuple())
            
        // get sample retained stats on number of reads with adaptor, subtelo pass, telomeric
        sample_stats = SUMMARY_STATS_SAMPLE(demuxed_reads.mix(subtelo_filtered_ch.retained_reads, telo_stats.retained_reads).groupTuple(), \
                                        demuxed_reads.mix(subtelo_filtered_ch.filtered_reads, telo_stats.no_telo_start, telo_stats.below_telo_threshold).groupTuple())
    
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

    params = getParams()
    versions = getVersions()
    manifest = getManifest()

    sample_stats.retained_stats.multiMap{ label, stats ->
            label: label
            stats: stats
        }.set { retained_sample }

    sample_stats.filtered_stats.multiMap{ label, stats ->
            label: label
            stats: stats
        }.set { filtered_sample }


    //generate_html_report(file(params.outdir), stats_done[1])
    GENERATE_FINAL_REPORT(params.params, versions.versions, manifest.manifest, \
                        run_stats.retained_stats, run_stats.filtered_stats, \
                        retained_sample.stats.collect(), filtered_sample.stats.collect() \
        )
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
    Pinguscript.ping_complete(nextflow, workflow, params)
}

workflow.onError {
    Pinguscript.ping_error(nextflow, workflow, params)
}
