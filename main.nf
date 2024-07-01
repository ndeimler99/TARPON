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
include { reverse_complementation } from "./bin/process.nf"
include { identify_tagging_adaptor } from "./bin/process.nf"
include { telo_start_identification } from "./bin/process.nf"
include { individual_read_plots } from "./bin/process.nf"
include { generate_plots } from "./bin/process.nf"
include { generate_detailed_plots } from "./bin/process.nf"
include { summary_stats } from "./bin/process.nf"
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

    // check that all other parameters are valid...
    
    
    // ###### TO DO #######
    // check if filtered_telo is passed in
    // else - filter telomeric reads
    // #############


    //standard pipeline
    // putative identification of telomeric sequences to limit dataset size
    putative_sequences = putative_isolation(input_file)

    // reverse complement C strands into G strands and modify header line to include strand information
    reverse_complemented = reverse_complementation(putative_sequences.putative_reads)

    // isolate reads with adaptor sequences and filter by subtelomere size
    adaptor_identified = identify_tagging_adaptor(reverse_complemented.reversed_reads)

    // this is where I would add demultiplexing so further processes can be ran in parallel
    //telo start and length determination
    // analyze reads and create stats file containing read_id, strand, read_len, VRR_Start, VRR_length, Telo_length
    telo_stats = telo_start_identification(adaptor_identified.adaptor_reads)

    //depending on which start is choosen create plots
    //if individual read, create plots
    if (params.indiv_read_plots) {
        individual_read_plots(telo_stats.telomeric_sequences, telo_stats.telomeric_stats)
    }

    generate_plots(telo_stats.telomeric_stats)

    summary_stats(telo_stats.telomeric_stats)
    
    // if (params.detailed_stats) {
    //     telo_stats = generate_detailed_plots(telo_stats.telomeric_stats, telo_stats.telomeric_sequences)
    // }

    // remove workdir if value true in nextflow config




    //load in files from channel/somehow split demultiplex output into different channels?
    /*
    Channel.fromPath("${params.metadata}") | splitCsv(header: true) | 
        map { row -> meta = [id:row.id, gff:file(row.gff), fasta_ref:file(row.fasta)]
                    [meta, file(row.gff)] } | set { gff_files }
    

    Channel.fromPath("${params.input_gff}/*.gff") . map { id, reads ->
        tokens = id.tokenize(".gff") } | view
    */ 
    
}

