/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Import Required Workflows and Processes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CREATE_CLUSTER_SPECIFIC_FASTA } from "../bin/process.nf"
include { GENERATE_CONSENSUS } from "../bin/process.nf"
include { MERGE_CONSENSUS } from "../bin/process.nf"
include { CONCORDANCE_ALIGNMENT } from "../bin/process.nf"
include { CREATE_GRAPH_FILE } from "../bin/process.nf"

// include { ALIGN_TO_PARENT as M_ALIGN 
//           ALIGN_TO_PARENT as P_ALIGN 
//         } from "../bin/process.nf"
        
// include { ASSIGN_PATERNITY } from "../bin/process.nf"

workflow cluster_concordance_pipeline {

    take:
        concordance_sample_file

    main:

        // load sample file to channel
        sample_ch = Channel
            .fromPath(file(concordance_sample_file))
            .splitCsv(header:true)
            .map { row -> tuple(row.sample, file(row.stats), file(row.telobam))}
        
        //sample_ch.view()
        // generate cluster fasta
        fasta_ch = CREATE_CLUSTER_SPECIFIC_FASTA(sample_ch)

        fasta_ch.fasta
            .flatMap { sample, files -> files.collect { file -> tuple(sample, file.baseName, file) } }
            .set { fasta_ch }

        //fasta_ch.view()
        // // generate consensus
        cons_seqs = GENERATE_CONSENSUS(fasta_ch)

        cons_seqs
            .map { sample, cluster, fasta ->
                tuple(sample, fasta)
            }
            .groupTuple()
            .set { cons_seqs }

        //cons_seqs.view()

        // // merge consensus
        merged_consensus = MERGE_CONSENSUS( cons_seqs )


        merged_consensus
            .map { sample, file -> file }
            .collect()
            .flatMap { files ->
                files.subsequences()
                .findAll { it.size() == 2 }
                .collect { pair -> tuple(pair[0], pair[1]) }
            }
        .set{ pairs_ch }

        alignments = CONCORDANCE_ALIGNMENT(pairs_ch)
        // align consensus to each other in pairwise manner

        alignments.collect().view()

        graph_file = CREATE_GRAPH_FILE(alignments.collect(), params.concordance_sample_file)
}