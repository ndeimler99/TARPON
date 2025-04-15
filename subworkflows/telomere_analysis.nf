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
    include { BASIC_PLOTS } from "../bin/process.nf"
    include { GENERATE_DETAILED_PLOTS } from "../bin/process.nf"
    include { SUMMARY_STATS_RUN } from "../bin/process.nf"
    include { SUMMARY_STATS_SAMPLE } from "../bin/process.nf"
    include { RESTRICTION_DIGEST_ANALYSIS } from "../bin/process.nf"
    include { GENERATE_FINAL_REPORT } from "../bin/process.nf"
    include { SEPERATE_STRANDS as SEPERATE_STRAND_RUN_FILTERING } from "../bin/process.nf"
    include { SEPERATE_STRANDS as SEPERATE_STRAND_SAMPLE_FILTERING } from "../bin/process.nf"
    include { SEPERATE_STRANDS as SEPERATE_STRAND_RUN_RETAINED } from "../bin/process.nf"
    include { SEPERATE_STRANDS as SEPERATE_STRAND_SAMPLE_RETAINED } from "../bin/process.nf"
    include { SUBTELO_FILTERING } from "../bin/process.nf"
    include { getParams } from "../bin/process.nf"
    include { getVersions } from "../bin/process.nf"
    include { getManifest } from "../bin/process.nf"
    include { COMBINE_BAM as COMBINED_RETAINED_BAM } from "../bin/process.nf"
    include { COMBINE_BAM as COMBINED_FILTERED_BAM } from "../bin/process.nf"
    include { validateParameters; paramsHelp; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'
    include { FINAL_TELO_STATS } from "../bin/process.nf"
    include { GET_EMPTY_CHANNEL } from "../bin/process.nf"
    include { MUTANT_ANALYSIS } from "../bin/process.nf"
    include { VARIANT_ANALYSIS } from "../bin/process.nf"
    include { PLOT_TELO_GRAPHS } from "../bin/process.nf"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Run Workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow telomere_analysis_pipeline {
    take:
        input_ch
        putative_reads
        non_telomeric
        sample_file

    main:
        // convert all C strand telomeric sequences to G strand for further code simplicity (this is tracked in header line)
        reversed_ch = REVERSE_COMPLEMENTATION(putative_reads)

        // check to see if adaptor sequence is provided and if reads need to be demultiplexed
        if (params.capture_probe_sequence == ""){
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

        if (params.mutant != "false") {
            mutant_analysis = MUTANT_ANALYSIS(telo_stats.final_telomeric)
        }
        else {
            mutant_analysis = VARIANT_ANALYSIS(telo_stats.final_telomeric)
        }

        PLOT_TELO_GRAPHS(mutant_analysis.sequences)
        // merge all stats relevant to retained telomeric sequences
        run_retained = COMBINED_RETAINED_BAM(subtelo_ch.retained_reads.multiMap { label, stats -> stats: stats }.collect(). map { it -> [ "subtelo_pass", it ]}.mix(
                                                    telo_stats.retained_reads.multiMap { label, stats -> stats: stats }.collect().map {it -> [ "telo_retained", it]}))

        // merge all stats relevant to filtered telomeric sequences
        run_removed = COMBINED_FILTERED_BAM(subtelo_ch.removed_reads.multiMap {label, stats -> stats:stats}.collect().map {it -> ["subtelo_fail", it]}.mix(
                                                    telo_stats.no_telo_start.multiMap { label, stats -> stats: stats}.collect().map{it -> ["no_telo_start", it]},
                                                    telo_stats.below_telo_threshold.multiMap { label, stats -> stats: stats }.collect().map{it -> ["below_telo_threshold", it]}))

        // separate reads into G and C strand and recalculate all statistics
        if (params.strand_comparison){
            separate_run_removed = SEPERATE_STRAND_RUN_FILTERING(adaptor_ch.removed_reads.mix(run_removed.combined))
            separate_run_retained = SEPERATE_STRAND_RUN_RETAINED(reversed_ch.retained_reads.mix(adaptor_ch.retained_reads, run_retained.combined))


            removed_sample = SEPERATE_STRAND_SAMPLE_FILTERING(subtelo_ch.removed_reads.mix(telo_stats.no_telo_start, telo_stats.below_telo_threshold))
            retained_sample = SEPERATE_STRAND_SAMPLE_RETAINED(demuxed_reads.mix(subtelo_ch.retained_reads, telo_stats.retained_reads))

            
            //get run retained stats on all reads input, putative reads, putative strand specific, adaptor found, adaptor strand specific.
            run_stats = SUMMARY_STATS_RUN(input_ch.mix(putative_reads, reversed_ch.retained_reads, run_retained.combined, separate_run_retained.g_strand, separate_run_retained.c_strand, adaptor_ch.retained_reads).groupTuple(), \
                                            input_ch.mix(non_telomeric, reversed_ch.removed_reads, run_removed.combined, adaptor_ch.removed_reads, separate_run_removed.c_strand, separate_run_removed.g_strand).groupTuple())


            // get sample retained stats on number of reads with adaptor, adaptor strand specific, subtelo pass, subtelo strand specific, telomeric, telomeric strand specfic
            sample_stats = SUMMARY_STATS_SAMPLE(demuxed_reads.mix(subtelo_ch.retained_reads, telo_stats.retained_reads, retained_sample.c_strand, retained_sample.g_strand).groupTuple(), \
                                            demuxed_reads.mix(subtelo_ch.removed_reads, telo_stats.no_telo_start, telo_stats.below_telo_threshold, removed_sample.c_strand, removed_sample.g_strand).groupTuple())

        }
        else {
            //get run retained stats on all reads input, putative reads, adaptor found
            run_stats = SUMMARY_STATS_RUN(input_ch.mix(putative_reads, reversed_ch.retained_reads, adaptor_ch.retained_reads, run_retained.combined).groupTuple(), \
                                        input_ch.mix(non_telomeric, reversed_ch.removed_reads, adaptor_ch.removed_reads, run_removed.combined).groupTuple())
                
            // get sample retained stats on number of reads with adaptor, subtelo pass, telomeric
            sample_stats = SUMMARY_STATS_SAMPLE(demuxed_reads.mix(subtelo_ch.retained_reads, telo_stats.retained_reads).groupTuple(), \
                                            demuxed_reads.mix(subtelo_ch.removed_reads, telo_stats.no_telo_start, telo_stats.below_telo_threshold).groupTuple())
        
        }

        //if individual reads is specified one plot per read showing telomeric percentage will be created
        if (params.indiv_read_plots) {
            INDIVIDUAL_READ_PLOTS(telo_stats.final_telomeric)
        }

        // this will generate basic statistical plots for the end user
        BASIC_PLOTS(telo_stats.final_telomeric)

        // if detailed stats is specified analyze telomeric sequences but include things such as percentage telomeric, etc.
        if (params.detailed_stats) {
            telo_stats = GENERATE_DETAILED_PLOTS(mutant_analysis.statistics)
        }

        // if a restriction digest sequence is provided run a restrition digest analysis
        if (params.restriction_digest_analysis != ""){
            restriction_digest = RESTRICTION_DIGEST_ANALYSIS(telo_stats.final_telomeric)
        }
        else {
            // simply returns an empty channel for use in report
            restriction_digest = GET_EMPTY_CHANNEL()
        }

        // maintenance details and compliation for report
        params = getParams()
        versions = getVersions()
        manifest = getManifest()

        sample_stats.retained_stats.multiMap{ label, stats ->
                label: label
                stats: stats
            }.set { retained_sample }
 

        sample_stats.removed_stats.multiMap{ label, stats ->
                label: label
                stats: stats
            }.set { removed_sample }


        mutant_analysis.statistics.multiMap{ label, stats ->
                label: label
                stats: stats
            }.set { mutant_analysis_stats }

        mutant_analysis.repeat_distribution.multiMap {label, stats ->
                label: label
                stats: stats
            }.set { mutant_analysis_repeat_distribution }

        mutant_analysis.processivity.multiMap {label, stats ->
                label: label
                stats: stats
            }.set { mutant_analysis_processivity }
        
        FINAL_TELO_STATS(mutant_analysis_stats.stats.collect())

        //generate_html_report(file(params.outdir), stats_done[1])
        report = GENERATE_FINAL_REPORT(params.params, versions.versions, manifest.manifest, \
                            run_stats.retained_stats, run_stats.removed_stats, \
                            retained_sample.stats.collect(), removed_sample.stats.collect(), \
                            mutant_analysis_stats.stats.collect(), \
                            FINAL_TELO_STATS.out.vrr_stats, \
                            restriction_digest.stats.collect(), \
                            mutant_analysis_repeat_distribution.stats.collect(), \
                            mutant_analysis_processivity.stats.collect()
        )
}
