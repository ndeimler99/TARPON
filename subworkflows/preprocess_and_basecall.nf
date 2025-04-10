/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Import Required Workflows and Processes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PUTATIVE_ISOLATION } from "../bin/process.nf"
include { COMBINE_FASTQ } from "../bin/process.nf"
include { COMBINE_FASTQ_GZ } from "../bin/process.nf"
include { FASTQ_2_FASTQGZ } from "../bin/process.nf"
include { COMBINED_FASTQ_GZ } from "../bin/process.nf"
include { CONVERT_BAM_2_FASTQ } from "../bin/process.nf"
include { COMBINE_BAM } from "../bin/process.nf"
include { ISOLATE_POD5_SQUIGGLES } from "../bin/process.nf"
include { PUTATIVE_ISOLATION as PUTATIVE_ISOLATION2 } from "../bin/process.nf"
include { BASECALLING } from "../bin/process.nf"
include { FASTQ_TO_BAM } from "../bin/process.nf"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Run Workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow preprocess_data_pipeline {
    take:
        run
        input

    main:
        // check inf input parameters is a directory or a file
            // if it is a directory look for fastq.gz, fastq, and then bam files and concatenate fastq.gz, concatenate and gzip fastq, or merge bams and convert to fastq.gz
            // if input file is a single file - check the file extension.  If it is a bam file - convert to fastq.gz, if it is fastq - gzip. Do nothing if it is a fastq.gz
        if (file(input).isDirectory()) {
            if (file("${input}/*.bam", checkIfExists: true) == []){
                if (file("${input}/*.fastq", checkIfExists: true) == []){
                    if (file("${input}/*.fastq.gz", checkIfExists: true) == []){
                        exit 1, "No Valid File Types Identified in Input Directory"
                    }
                    else {
                        fastq_combined = COMBINE_FASTQ_GZ(Channel.fromPath ( "${input}/*fastq.gz" ).collect().map{ it -> [run, it]})
                        input_ch = FASTQ_TO_BAM(fastq_combined)
                    }
                }
                else{
                    fastq_combined = COMBINE_FASTQ(Channel.fromPath ( "${input}/*fastq" ).collect().map{ it -> [run, it]})
                    //convert to bam
                    input_ch = FASTQ_TO_BAM(fastq_combined)
                }
            }
            else {
                input_ch = COMBINE_BAM(Channel.fromPath ( "${input}/*bam" ).collect().map{ it -> [run, it]})
            }
        }

        else{
            if (file(params.input).extension == "bam") {
                Channel.fromPath( params.input, checkIfExists:true)
                    .map{ it -> [ run , it] }
                    .set{ input_ch }
            }
            else if (file(params.input).extension == "fastq") {
                //fastqgz = FASTQ_2_FASTQGZ(Channel.fromPath("${params.input}", checkIfExists:true).collect().map{it -> ["input", it]})
                //convert to bam
                input_ch = FASTQ_TO_BAM(Channel.fromPath("${params.input}", checkIfExists:true).collect().map{it -> [run, it]})
            }
            else {
                //convert to bam
                input_ch = FASTQ_TO_BAM(Channel.fromPath("${params.input}", checkIfExists:true).collect().map{it -> [run, it]})
            }
        }    



        // putative identification of telomeric sequences to limit dataset size
        putative_ch = PUTATIVE_ISOLATION(input_ch)

        if (params.fast_basecalled) {
            // rebasecall pod5 reads
            // this will use pod5 filter -ids
            reads_to_basecall = ISOLATE_POD5_SQUIGGLES(putative_ch.putative_reads, file(params.pod5_dir))
            reduced_input = BASECALLING(reads_to_basecall.pod5_filtered)
            putative_ch = PUTATIVE_ISOLATION2(reduced_input)
        }

    emit:
        input = putative_ch.input_ch
        putative_reads = putative_ch.putative_reads
        non_telomeric = putative_ch.non_telomeric
}
