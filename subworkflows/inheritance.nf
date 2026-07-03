/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Import Required Workflows and Processes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include {   CREATE_CLUSTER_SPECIFIC_FASTA as M_CLUS_FA
            CREATE_CLUSTER_SPECIFIC_FASTA as P_CLUS_FA
            CREATE_CLUSTER_SPECIFIC_FASTA as O_CLUS_FA 
        } from "../bin/process.nf"
        
include {   GENERATE_CONSENSUS as M_CONS_GEN
            GENERATE_CONSENSUS as P_CONS_GEN
            GENERATE_CONSENSUS as O_CONS_GEN
        } from "../bin/process.nf"

include {  MERGE_CONSENSUS as M_MERGE_CONS
            MERGE_CONSENSUS as P_MERGE_CONS
            MERGE_CONSENSUS as O_MERGE_CONS
        } from "../bin/process.nf"
       
include { ALIGN_TO_PARENT as M_ALIGN 
          ALIGN_TO_PARENT as P_ALIGN 
        } from "../bin/process.nf"
        
include { ASSIGN_PATERNITY } from "../bin/process.nf"

workflow inheritance_assignment_pipeline {

    take:
        maternal_telobam
        paternal_telobam
        offspring_telobam
        maternal_stats
        paternal_stats
        offspring_stats

    main:

        maternal_fastas = M_CLUS_FA(tuple "maternal", file(maternal_stats), file(maternal_telobam))

        //maternal_fastas.view()
        maternal_fastas.fasta
            .flatMap { sample, files -> files.collect { file -> tuple(sample, file.baseName, file) } }
            .set { maternal_fasta_ch }

        maternal_consensus = M_CONS_GEN(maternal_fasta_ch)

        maternal_consensus
            .map { sample, cluster, fasta ->
                tuple(sample, fasta)
            }
            .groupTuple()
            .set { maternal_consensus }

        maternal_consensus = M_MERGE_CONS( maternal_consensus )

        maternal_consensus.view()

        paternal_fastas = P_CLUS_FA(tuple "paternal", file(paternal_stats), file(paternal_telobam))

        paternal_fastas.fasta
            .flatMap { sample, files -> files.collect { file -> tuple(sample, file.baseName, file) } }
            .set { paternal_fasta_ch }
            
        paternal_consensus = P_CONS_GEN(paternal_fasta_ch)

        paternal_consensus
            .map { sample, cluster, fasta ->
                tuple(sample, fasta)
            }
            .groupTuple()
            .set { paternal_consensus }

        paternal_consensus = P_MERGE_CONS( paternal_consensus )

        offspring_fastas = O_CLUS_FA(tuple "offspring", file(offspring_stats), file(offspring_telobam))

        offspring_fastas.fasta
            .flatMap { sample, files -> files.collect { file -> tuple(sample, file.baseName, file) } }
            .set { offspring_fasta_ch }
            
        offspring_consensus = O_CONS_GEN(offspring_fasta_ch)

        offspring_consensus
            .map { sample, cluster, fasta ->
                tuple(sample, fasta)
            }
            .groupTuple()
            .set { offspring_consensus }

        offspring_consensus = O_MERGE_CONS(offspring_consensus)

        maternal_alignment = M_ALIGN(maternal_consensus.consensus_fa, offspring_consensus.consensus_fa)
        paternal_alignment = P_ALIGN(paternal_consensus.consensus_fa, offspring_consensus.consensus_fa)
        
        ASSIGN_PATERNITY(maternal_alignment.aln_results, paternal_alignment.aln_results,
                            maternal_consensus.consensus_fa, paternal_consensus.consensus_fa,
                            offspring_consensus.consensus_fa, offspring_stats)
        // align, filter, and assign
        // plots
}