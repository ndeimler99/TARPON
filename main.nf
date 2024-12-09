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

include { validate_parameters } from "./subworkflows/parameter_validation.nf"
include { preprocess_data_pipeline } from "./subworkflows/preprocess_and_basecall.nf"
include { telomere_analysis_pipeline } from "./subworkflows/telomere_analysis.nf"
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

    // // check if output directory exists
    // outdir_check = Channel.fromPath("${params.outdir}/report.html").ifEmpty(false)
    
    // outdir_check.view()
    valid_params = validate_parameters()

    if (valid_params.passed.value == false){
        exit 1, "Parameter Validation Failed"
    }

    run = params.run_name

    if (params.real_time) {
        real_time_pipeline()
    }
    else {
        // pipeline that takes input, appropriately process it, and returns telomeric sequences
        preprocess_data_pipeline(params.run_name, params.input_file)
        telomere_analysis_pipeline(preprocess_data_pipeline.out, params.sample_file)


    // ###### TO DO #######
    // check if filtered_telo is passed in
    // else - filter telomeric reads

    // reference_file for mapping has to be value channel (just leave it at params.reference)
    // #############
    
    // convert all bam files to fastq
    }
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
