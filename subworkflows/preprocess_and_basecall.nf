
include { PUTATIVE_ISOLATION } from "../bin/process.nf"
include { COMBINE_FASTQ } from "../bin/process.nf"
include { FASTQ_2_FASTQGZ } from "../bin/process.nf"
include { COMBINED_FASTQ_GZ as COMBINED_FASTQ_GZ } from "../bin/process.nf"
include { CONVERT_BAM_2_FASTQ } from "../bin/process.nf"

workflow preprocess_data_pipeline {
    take:
        run
        input_file

    main:

        if (file(input_file).isDirectory()) {
            try {
                input_ch = COMBINED_FASTQ_GZ(Channel.fromPath ( "${input_file}/*q.gz" ).collect().map{ it -> ["input", it]})
            }
            
        }

        else{
            if (file(params.input_file).extension == "bam") {
                print(file(params.input_file).extension)
                input_ch = CONVERT_BAM_2_FASTQ(Channel.fromPath("${input_file}", checkIfExists: true).collect().map{ it -> ["input", it]})
            }
            else if (file(params.input_file).extension == "fastq") {
                input_ch = FASTQ_2_FASTQGZ(Channel.fromPath("${input_file}", checkIfExists:true).collect().map{it -> ["input", it]})
            }
            else {
                Channel.fromPath( params.input_file, checkIfExists:true)
                .map{ it -> [ run , it] }
                .set{ input_ch }
            }
        }    


        // putative identification of telomeric sequences to limit dataset size
        putative_ch = PUTATIVE_ISOLATION(input_ch)

    emit:
        input = input_ch
        putative_reads = putative_ch.putative_reads
        non_telomeric = putative_ch.non_telomeric
}