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
    v1.0.3
    """.stripIndent()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Import Required Workflows and Processes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { validate_parameters } from "./subworkflows/parameter_validation.nf"
include { preprocess_data_pipeline } from "./subworkflows/preprocess_and_basecall.nf"
include { telomere_isolation_pipeline } from "./subworkflows/telomere_isolation.nf"
include { telomere_analysis_pipeline } from "./subworkflows/telomere_analysis.nf"
include { clustering_pipeline } from "./subworkflows/clustering_analysis.nf"
include { paramsHelp; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'
include { getParams; getVersions; getManifest; GENERATE_FINAL_REPORT } from "./bin/process.nf"
include { enrichment_stats_pipeline } from "./subworkflows/enrichment_stats.nf"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Run Workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
WorkflowMain.initialise(workflow, params, log)

workflow {

    // validate parameters and throw an error if invalid parameters
    valid_params = validate_parameters()

    parameters = getParams()
    versions = getVersions()
    manifest = getManifest()
    
    if (valid_params.passed.value == false){
        exit 1, "Parameter Validation Failed"
    }

    // check if pipeline should be running while sequencing - this function currently does not work and will result in no output being generated
    if (params.real_time) {
        // real_time_pipeline()
    }
    else {
        // preprocess data pipeline takes the input files or directory and returns SUP basecalled telomeric sequences
        preprocess_data_pipeline(params.run_name, params.input)

        // takes putative telomeric sequences returned by preprocess data pipeline and runs all relevant processes to generate descriptive stats and report.html
        telomere_isolation_pipeline(preprocess_data_pipeline.out, params.sample_file)

        //telomere_isolation.collect().view()
        telomere_stats = telomere_analysis_pipeline(telomere_isolation_pipeline.out)

        enrichment_stats = enrichment_stats_pipeline(telomere_isolation_pipeline.out)
        
        if (params.clustering) {
            clustering_results = clustering_pipeline(telomere_isolation_pipeline.out)
        }

        report = GENERATE_FINAL_REPORT(parameters.params, versions.versions, manifest.manifest, \
                            enrichment_stats.flowcell_retained, enrichment_stats.flowcell_removed, \
                            telomere_stats.telomere_stats, \
                            telomere_stats.sample_specific_stats.mix(enrichment_stats.sample_specific_stats).collect()
        )
    }

}

// When workflow finishes return basic description of finished or not and if it works remove the work directory if specified during run
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
