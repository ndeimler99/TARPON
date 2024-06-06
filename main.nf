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
include { putative_isolation } from "./bin/process.nf"
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    // check if input file exists
    input_file = file(params.input_file)
    if (!input_file.exists() ) {
        exit 1, "Input File Does Not Exist"
    }

    // check if output directory exists
    outdir = file(params.outdir)
    if (outdir.exists() && !params.overwrite_outdir) {
        exit 1, "Out Directory Already Exists, Please Provide New Out Directory Name or Allow Overwriting of Pre-existing directory"
    }
    
    // ###### TO DO #######
    // check if filtered_telo is passed in
    // else - filter telomeric reads
    // #############


    //standard pipeline
    // putative identification of telomeric sequences to limit dataset size
    putative_sequences = putative_isolation(input_file)


   
    //load in files from channel/somehow split demultiplex output into different channels?
    /*
    Channel.fromPath("${params.metadata}") | splitCsv(header: true) | 
        map { row -> meta = [id:row.id, gff:file(row.gff), fasta_ref:file(row.fasta)]
                    [meta, file(row.gff)] } | set { gff_files }
    

    Channel.fromPath("${params.input_gff}/*.gff") . map { id, reads ->
        tokens = id.tokenize(".gff") } | view
    */


    
    
}

