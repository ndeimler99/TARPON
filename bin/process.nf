import groovy.json.JsonOutput


process PUTATIVE_ISOLATION {

    tag "$run_name - Putative Isolation"
    label 'tarpon'

    input:
        tuple val(run_name), path(reads_file, stageAs: "input.fastq.gz")
    
    output:
        tuple val(run_name), path ("putative_reads.fastq"), emit: putative_reads
        tuple val (run_name), path("non_telomeric.fastq"), emit: non_telomeric
        tuple val(run_name), path("input.fastq.gz"), emit: input_ch
    
    //publishDir "${params.outdir}/TELOMERIC/", overwrite: true, mode: 'copy', pattern: "input.fastq.gz"
    //publishDir "${params.outdir}/FILTERED_READS/", overwrite: true, mode: 'copy', pattern: "input.fastq.gz"
    //publishDir "${params.outdir}/TELOMERIC/", mode: 'copy', overwrite: true, pattern: "putative_reads.fastq"
    //publishDir "${params.outdir}/FILTERED_READS/",  overwrite: true, pattern: "non_telomeric.fastq"

    script:
    """
    isolate_putative_telomeric_reads.py --input_file ${reads_file} --repeat ${params.repeat} --repeat_count ${params.repeat_count} --c_strand_only ${params.c_strand_only} --out_file putative_reads.fastq --non_telo non_telomeric.fastq
    """
}

process SEPERATE_STRANDS {

    label 'tarpon'
    tag "$id - Seperating $reads.baseName"
    input:
        tuple val(id), path(reads)
    
    output:
        tuple val(id), path("${reads.baseName}.g_strand.fastq"), emit: g_strand
        tuple val(id), path("${reads.baseName}.c_strand.fastq"), emit: c_strand

    script:
    """
    separate_strands.py --input_file ${reads} --g_file ${reads.baseName}.g_strand.fastq --c_file ${reads.baseName}.c_strand.fastq
    """

}
process REVERSE_COMPLEMENTATION {
    
    label 'tarpon'
    tag "$run_name - Reverse Complementation"
    input:
        tuple val(run_name), path(reads)

    output:
        tuple val(run_name), path("20_80_removed_reads.fastq"), emit: removed_reads
        tuple val(run_name), path("putative_reads.filtered.fastq"), emit: retained_reads
    
    //publishDir "${params.outdir}/FILTERED_READS/", mode: 'copy', overwrite: true, pattern: "20_80_removed*.fastq"
    //publishDir "${params.outdir}/TELOMERIC/", mode: 'copy', overwrite: true, pattern: "putative*.fastq"

    script:
    """
    reverse_complement_reads.py --input_file ${reads} \
        --repeat ${params.repeat} \
        --threshold ${params.reverse_complement_threshold} \
        --c_strand_only ${params.c_strand_only} \
        --out_file putative_reads.filtered.fastq \
        --removed_reads 20_80_removed_reads.fastq
    """
}

process IDENTIFY_TAGGING_ADAPTOR_AND_DEMUX {
    
    label 'tarpon'
    tag "$run_name - Identify Adaptor and Demultiplexing"

    input:
        tuple val(run_name), path(reads)
        path(barcodes_file)

    output:
        path("DEMUX/*.fastq"), emit: demuxed_reads
        tuple val(run_name), path("adaptor.fastq"), emit: retained_reads
        tuple val(run_name), path("adaptor_filtered.fastq"), emit: filtered_reads

    //publishDir "${params.outdir}/", mode: 'copy', overwrite: true, pattern: "FILTERED_READS/*"
    //publishDir "${params.outdir}/", mode: 'copy', overwrite: true, pattern: "TELOMERIC/*"

    script:
    if (params.adaptor_sequence == "")
        // no adaptor sequence provided - using ONT like approach were telomere overhang goes directly into barcode
        // computationally not ideal
        """
        mkdir DEMUX/
        identify_tagging_barcodes.py --input_file ${reads} --sample_file ${barcodes_file} --barcode_errors ${params.barcode_errors} --repeat ${params.repeat} --out_fh DEMUX/ --no_adaptor adaptor_filtered.fastq
        cat DEMUX/* > adaptor.fastq
        """
    else
        // multiple samples
        // run adaptor identification and then demux
        """
        mkdir TELOMERIC
        mkdir TELOMERIC/DEMUX/
        mkdir FILTERED_READS
        """
    // else
        // multiple samples
            // use same scripts as with no sample file with additional step of sample assignment
    
}

process IDENTIFY_TAGGING_ADAPTOR {
    
    label 'tarpon'
    tag "$run_name - Identify Adaptor and Demultiplexing"

    input:
        tuple val(run_name), path(reads)

    output:
        path("${params.sample_name}.fastq"), emit: demuxed_reads
        tuple val(run_name), path("adaptor.fastq"), emit: retained_reads
        tuple val(run_name), path("adaptor_filtered.fastq"), emit: filtered_reads

    //publishDir "${params.outdir}/", mode: 'copy', overwrite: true, pattern: "FILTERED_READS/*"
    //publishDir "${params.outdir}/", mode: 'copy', overwrite: true, pattern: "TELOMERIC/*"

    script:
    """
    mkdir DEMUX/
    identify_tagging_adaptor.py --input_file ${reads} \
        --adaptor_sequence ${params.adaptor_sequence} \
        --adaptor_errors ${params.adaptor_sequence_errors} \
        --repeat ${params.repeat} \
        --adaptor_found ${params.sample_name}.fastq \
        --no_adaptor adaptor_filtered.fastq

    cp ${params.sample_name}.fastq adaptor.fastq
    """
}

process SUBTELO_FILTERING {
    label 'tarpon'
    tag "$sample - Subtelomeric Filtering"
    input:
        tuple val(sample), path(reads)
    
    output:
        tuple val(sample), path("*subtelo_fail.fastq"), emit: filtered_reads
        tuple val(sample), path("*subtelo_pass.fastq"), emit: retained_reads

    script:
    """
    filter_by_subtelo.py --input_file ${reads} \
        --min_subtelo_length ${params.min_subtelo_length} \
        --min_subtelo_threshold ${params.subtelo_threshold} \
        --repeat ${params.repeat} \
        --passes_subtelo ${sample}.subtelo_pass.fastq \
        --fails_subtelo ${sample}.subtelo_fail.fastq

    """
}


process TELO_START_IDENTIFICATION {
    label 'tarpon'
    tag "$sample - Identifying Telomere Start"
    input:
        tuple val(sample), path(reads)

    output:
        tuple val(sample), path("*telomeric.fastq"), path("*telomeric_stats.txt"), emit: final_telomeric
        path("*telomeric_stats.txt"), emit: final_telo_stats
        tuple val(sample), path("*no_telomere_start.fastq"), emit: no_telo_start
        tuple val(sample), path("*.below_telo_%_threshold.fastq"), emit: below_telo_threshold
        tuple val(sample), path("*telomeric.fastq"), emit: retained_reads

    //publishDir "${params.outdir}/FINAL_TELO/", mode:'copy', overwite:true, pattern: "telomeric_stats.txt", saveAs: { filename -> "${sample}.stats.txt" }
    //publishDir "${params.outdir}/FINAL_TELO/", mode:'copy', overwite:true, pattern: "TELOMERIC/*", saveAs: { filename -> "${sample}.telomeric.fastq" }
    //publishDir "${params.outdir}/", mode: 'copy', overwrite: true, pattern: "FILTERED_READS/*"
    //publishDir "${params.outdir}/", mode: 'copy', overwrite: true, pattern: "TELOMERIC/*"

    publishDir "${params.outdir}/${sample}/", mode: 'copy', overwrite:true, pattern:"*.telomeric.fastq"
    publishDir "${params.outdir}/${sample}/", mode: 'copy', overwrite:true, pattern:"telomeric_stats.txt"
    
    script:
    """
    #python script that identified telomere start. Writes out fastq file, stats file, fastq for reads removed because no telo start was found, fastq for reads removed because didnt reach minimum threshold
    mkdir TELOMERIC
    mkdir FILTERED_READS
    identify_telo_start.py --input_file ${reads} --repeat ${params.repeat} --sliding_window ${params.sliding_window_size} \
        --sliding_window_interval ${params.sliding_window_interval} \
        --upper_threshold ${params.upper_threshold} \
        --lower_threshold ${params.lower_threshold} \
        --telomeric_rep_perc ${params.telomeric_repeat_percentage} \
        --consecutive_repeats ${params.consecutive_repeats} \
        --consecutive_threshold ${params.consecutive_threshold} \
        --telomeric_fastq_out ${sample}.telomeric.fastq \
        --no_telomere_out ${sample}.no_telomere_start.fastq \
        --filtered_out ${sample}.below_telo_%_threshold.fastq \
        --stats_fh ${sample}.telomeric_stats.txt
    """
}

process INDIVIDUAL_READ_PLOTS {

    label 'tarpon'
    tag "$sample - Plotting Individual Reads"

    input:
        tuple val(sample), path(reads), path(stats)

    output:
        path("*.pdf")

    publishDir "${params.outdir}/${sample}/FIGURES/INDIVIDUAL_READ_PLOTS/", mode:'move', overwrite: true, pattern: "*.pdf"

    script:
    """
    indiv_read_plots.py --input_file ${reads} --repeat ${params.repeat} --telo_stats ${stats} --sliding_window ${params.sliding_window_size} --sliding_window_interval ${params.sliding_window_interval}
    """
}

process GENERATE_PLOTS {

    label 'tarpon'
    tag "$sample - Generating Output Plots"

    input:
        tuple val(sample), path(telo_reads), path(telo_stats)
    
    output:
        path("*.pdf")
        path("C_G_COMPARISON/*.pdf"), optional:true
    
    publishDir "${params.outdir}/${sample}/FIGURES/", mode:'move', overwrite: true, pattern: "*.pdf"
    publishDir "${params.outdir}/${sample}/FIGURES/", mode:'move', overwrite: true, pattern: "C_G_COMPARISON/*.pdf"

    script:
    """
    basic_plots.R ${telo_stats} ${params.plot_telo_length} ${params.plot_vrr_length} ${params.strand_comparison}
    """
}

process GENERATE_DETAILED_PLOTS {

    label 'tarpon'
    tag "$sample - Generating Detailed Output Plots"

    input:
        tuple val(sample), path(telomeric_reads), path(telo_stats, stageAs: "old_telo_stats.txt")
    
    output:
        tuple val(sample), path(telomeric_reads), path("*telomeric_stats.txt"), emit: final_telomeric
        path("*telomeric_stats.txt"), emit: final_telo_stats
        //path("DETAILED_STATS/*.pdf")
        //path("C_G_COMPARISON/*.pdf"), optional: true
    
    publishDir "${params.outdir}/${sample}/FIGURES/", mode:'move', overwrite: true, pattern: "DETAILED_STATS/*.pdf"
    publishDir "${params.outdir}/${sample}/FIGURES", mode:'move', overwrite: true, pattern: "C_G_COMPARISON/*.pdf"
    publishDir "${params.outdir}/${sample}/", mode:'copy', overwrite: true, pattern: "telomeric_stats.txt", saveAs: { filename -> "${sample}.stats.txt" }

    script:
    """
    detailed_stats.py --fastq_file ${telomeric_reads} --repeat ${params.repeat} --stats_in ${telo_stats} --stats_out ${sample}.telomeric_stats.txt
    #detailed_plots.R telomeric_stats.txt ${params.plot_telo_length} ${params.plot_vrr_length} ${params.strand_comparison}
    """

}

process SUMMARY_STATS_RUN {

    label 'tarpon'
    tag "$id - Collecting Run Summary Statistics"

    input:
        tuple val(id), path(retained, stageAs: "RETAINED/*")
        tuple val(id), path(filtered, stageAs: "FILTERED/*")
        //path(output_dir)
    
    output:
        tuple val(id), path("Retained_Reads.stats.txt"), emit: retained_stats
        tuple val(id), path("Filtered_Reads.stats.txt"), emit: filtered_stats
        path("*.pdf")
       
    publishDir "${params.outdir}/RUN_STATS/", mode:'copy', overwrite:true, pattern:"*stats.txt"
    publishDir "${params.outdir}/RUN_STATS/FIGURES/", mode:'move', overwrite:true, pattern:"*pdf"


    script:
    """
    seqkit stats -a -N 50,90 -T ${retained} > Retained_Reads.stats.txt
    seqkit stats -a -N 50,90 -T ${filtered} > Filtered_Reads.stats.txt
    summary_stats_plots.R Retained_Reads.stats.txt ${params.strand_comparison} telomeric
    summary_stats_plots.R Filtered_Reads.stats.txt ${params.strand_comparison} filtered
    """
}

process SUMMARY_STATS_SAMPLE {

    label 'tarpon'
    tag "$id - Collecting Sample Summary Statistics"

    input:
        tuple val(id), path(retained, stageAs: "RETAINED/*")
        tuple val(id), path(filtered, stageAs: "FILTERED/*")
        //path(output_dir)
    
    output:
        tuple val(id), path("*retained.stats.txt"), emit: retained_stats
        tuple val(id), path("*.filtered.stats.txt"), emit: filtered_stats
    
    publishDir "${params.outdir}/${id}", mode:'copy', overwrite:true, pattern:"*stats.txt"
    
    script:
    """
    seqkit stats -a -N 50,90 -T ${retained} > ${id}.retained.stats.txt
    seqkit stats -a -N 50,90 -T ${filtered} > ${id}.filtered.stats.txt
    """
}

process RESTRICTION_DIGEST_ANALYSIS {

    label 'tarpon'
    tag "$sample - Performing Restriction Digest Analysis"

    input:
        tuple val(sample), path(telo_sequences), path(telo_stats)
    
    output:
        path("*digest_stats.txt")
    
    publishDir "${params.outdir}/${sample}/", mode: 'copy', overwrite: true, pattern: "*digest_stats.txt"

    script:
    """
    for seq in \$(echo "${params.restriction_digest_analysis}" | tr "," "\n"); do seqkit grep -s -p \$seq ${telo_sequences} > \$seq.fastq; done
    seqkit stats -a -N 50,90 -T *.fastq > ${sample}.digest_stats.txt
    """
}

process GENERATE_FINAL_REPORT {

    label 'tarpon'
    tag "Generating Final HTML Report"

    input:
        path("params.json")
        path("versions.txt")
        path("manifest.json")
        tuple val(run), path(stats_run_retained)
        tuple val(run1), path(stats_run_filtered)
        path(sample_stats_retained)
        path(sample_stats_filtered)
        path(telo_stats_per_sample)
        path(telo_descriptive_stats)
        path(vrr_descriptive_stats)
        path(restriction_digest)

    output:
        path("report.html")
    
    publishDir "${params.outdir}/", mode: 'copy', overwrite: true, pattern: "report.html"

    script:
    """
    
    generate_html_report.py --workflow_name TArPON \
                            --report report.html \
                            --template_file ${baseDir}/bin/single_sample_template.html \
                            --params params.json \
                            --versions versions.txt \
                            --manifest manifest.json \
                            --commandLine "${workflow.commandLine}" \
                            --run_stats_retained ${stats_run_retained} \
                            --run_stats_filtered ${stats_run_filtered} \
                            --sample_stats_retained ${sample_stats_retained} \
                            --sample_stats_filtered ${sample_stats_filtered} \
                            --sample_telo_stats ${telo_stats_per_sample} \
                            --run_telo_stats ${telo_descriptive_stats} \
                            --run_vrr_stats ${vrr_descriptive_stats} \
                            --restriction_digest ${restriction_digest} \
                            --plot_telo_length ${params.plot_telo_length} \
                            --plot_vrr_length ${params.plot_vrr_length} \
                            --strand_comparison ${params.strand_comparison} \
                            --detailed_stats ${params.detailed_stats}
    """
}

process getParams {

    label "tarpon"
    tag "Getting Parameters"

    output:
        path "params.json", emit:params

    script:
    json_str = JsonOutput.toJson(params)
    json_indented = JsonOutput.prettyPrint(json_str)
    // NOTE: single quotes are critical here;
    """
    echo '${json_indented}' > "params.json"
    """
}

process GET_EMPTY_CHANNEL {

    label "tarpon"
    tag "Getting Empty Channel"

    output:
        path("false"), emit: stats
    
    script:
    """
    touch false
    """

}

process getVersions {

    label "tarpon"
    tag "Getting Versions"

    output:
        path "versions.txt", emit: versions

    script:
    """
    python --version | sed 's/ /,/' >> versions.txt
    python -c "import regex; print(f'regex,{regex.__version__}')" >> versions.txt
    python -c "import pandas; print(f'pandas,{pandas.__version__}')" >> versions.txt
    seqkit version | sed 's/ /,/' >> versions.txt
    """
}

process getManifest {
    
    label 'tarpon'
    tag "Collecting Manifest Data"

    output:
        path "manifest.json", emit:manifest

    script:
    json_str = JsonOutput.toJson(workflow.manifest)
    json_indented = JsonOutput.prettyPrint(json_str)
    """
    echo '${json_indented}' > "manifest.json"
    """
}

process COMBINE_FASTQ {
    
    label 'tarpon'
    tag "$file_type Concatenating FASTQ Files"

    input:
        tuple val(file_type), path(input_files)

    output:
        tuple val(params.run_name), path("${file_type}.fastq"), emit:combined

    script:
    """
    cat ${input_files} > ${file_type}.fastq
    """
}

process FASTQ_2_FASTQGZ {
    label 'tarpon'
    tag "$file_type Converting FASTQ to FASTQ.GZ"

    input:
        tuple val(file_type), path(input_file)
    
    output:
        tuple val(params.run_name), path("${file_type}.fastq.gz"), emit: combined

    script:
    """
    gzip $input_file > ${file_type}.fastq.gz
    """
}

process CONVERT_BAM_2_FASTQ {

    label 'tarpon'
    tag "$file_type Converting FASTQ to FASTQ.GZ"

    input:
        tuple val(file_type), path(input_file)
    
    output:
        tuple val(params.run_name), path("${file_type}.fastq.gz"), emit: combined

    script:
    """
    samtools fastq -@ 4 $input_file > ${file_type}.fastq.gz
    """

}

process COMBINED_FASTQ_GZ {
    
    label 'tarpon'
    tag "$file_type Concatenating FASTQ Files"

    input:
        tuple val(file_type), path(input_files)

    output:
        tuple val(params.run_name), path("${file_type}.fastq.gz"), emit:combined

    script:
    """
    cat ${input_files} > ${file_type}.fastq.gz
    """
}

process BARCODE_HAMMING_CHECK {

    label 'tarpon'
    tag "$params.run_name Checking Barcode Hamming Distance"

    input:
        path(sample_file)
    
    output:
        path("passed.txt"), optional:true

    script:
    """
    check_hamming_distance.py --sample_file ${sample_file} --barcode_errors ${params.barcode_errors}
    """
}


process FINAL_TELO_STATS {
    
    label 'tarpon'
    tag "$params.run_name Final Plots and Telomere Stats Per Sample"

    input:
        path(input_files)

    output:
        path("sample_stats.txt"), emit: stats
        path("sample_stats.VRR.txt"), emit: vrr_stats
    //    path("combined_df.csv"), emit: combined_df
    //    path("*.pdf")


    publishDir "${params.outdir}/", overwrite: true, mode: 'copy', pattern: "sample_stats.txt"
    publishDir "${params.outdir}/", overwrite: true, mode: 'copy', pattern: "sample_stats.VRR.txt"
    //publishDir "${params.outdir}/", overwrite: true, mode: 'copy', pattern: "*.pdf"


    script:
    """
    processTelomereStats.py --stat_files ${input_files} --vrr_length ${params.plot_vrr_length} --telo_length ${params.plot_telo_length}
    #combinedPlots.R combined_df.csv ${params.plot_vrr_length}
    """

}