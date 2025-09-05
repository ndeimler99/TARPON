/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Import Required Workflows and Processes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    include { REVERSE_COMPLEMENTATION } from "../bin/process.nf"
    include { IDENTIFY_TAGGING_CAPTURE_PROBE_AND_DEMUX } from "../bin/process.nf"
    include { IDENTIFY_TAGGING_CAPTURE_PROBE } from "../bin/process.nf"
    include { TELO_START_IDENTIFICATION } from "../bin/process.nf"
    include { INDIVIDUAL_READ_PLOTS } from "../bin/process.nf"
    include { SUBTELO_FILTERING } from "../bin/process.nf"
    include { FINAL_TELO_STATS } from "../bin/process.nf"
    include { NO_CAPTURE_PROBE } from "../bin/process.nf"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Run Workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow telomere_isolation_pipeline {
    take:
        input_ch
        putative_reads
        non_telomeric
        sample_file

    main:

        if (params.methylation){
            modkit_conversion = BAM_TO_MOD_TABLE(telomeric_reads_with_stats)
        }
        // convert all C strand telomeric sequences to G strand for further code simplicity (this is tracked in header line)
        reversed_ch = REVERSE_COMPLEMENTATION(putative_reads)

        // check to see if adaptor sequence is provided and if reads need to be demultiplexed
        if (params.no_capture_probe){
            adaptor_ch = NO_CAPTURE_PROBE(reversed_ch.retained_reads)
        }
        else if (params.capture_probe_sequence == ""){
            // if no adaptor sequence is provided - reads were multiplexed in the ONT fashion where the barcode sits immediately adjacent to the telomeric sequence
            adaptor_ch = IDENTIFY_TAGGING_CAPTURE_PROBE_AND_DEMUX(reversed_ch.retained_reads, file(sample_file))
        }
        else if (params.sample_file == ""){
            // if no sample file was provided reads are not multiplexed and only the adaptor sequence needs to be identified
            adaptor_ch = IDENTIFY_TAGGING_CAPTURE_PROBE(reversed_ch.retained_reads)
        }
        else {
            // both an adaptor sequence and sample file were found - telomere end will be determined by adaptor sequence (common among all reads) and then demultiplexed based on downstream sequence
            adaptor_ch = IDENTIFY_TAGGING_CAPTURE_PROBE_AND_DEMUX(reversed_ch.retained_reads, file(sample_file))
        }

        // map output of adaptor_ch to demuxed_reads
        adaptor_ch.demuxed_reads.flatten()
            .map { it -> [it.baseName, it]}
            .set { demuxed_reads }

        // filter demuxed reads by subtelo length 
        subtelo_ch = SUBTELO_FILTERING(demuxed_reads)

        //telo start and length determination
        // analyze reads and create stats file containing read_id, strand, read_len, VRR_Start, VRR_length, Telo_length, and read quality
        telo_stats = TELO_START_IDENTIFICATION(subtelo_ch.retained_reads)

        //input_ch.view()
   emit:
        input_ch = input_ch
        putative_reads = putative_reads
        non_telomeric = non_telomeric
        sample_file = sample_file
        reversed_retained = reversed_ch.retained_reads
        reversed_removed = reversed_ch.removed_reads
        adaptor_retained_all_samples = adaptor_ch.retained_reads
        adaptor_removed_all_samples = adaptor_ch.removed_reads
        demuxed_reads = demuxed_reads
        subtelomeric_retained = subtelo_ch.retained_reads
        subtelomeric_removed = subtelo_ch.removed_reads
        telomeric_reads_with_stats = telo_stats.final_telomeric
        telomeric_reads_retained = telo_stats.retained_reads
        telomeric_reads_no_start = telo_stats.no_telo_start
        telomeric_reads_under_threshold = telo_stats.below_telo_threshold

}
