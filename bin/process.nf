import groovy.json.JsonOutput

process PUTATIVE_ISOLATION {

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
    isolate_putative_telomeric_reads.py ${reads_file} ${params.repeat} ${params.repeat_count} ${params.c_strand_only} putative_reads.fastq non_telomeric.fastq
    """
}

process SEPERATE_STRANDS {

    label 'tarpon'

    input:
        tuple val(id), path(reads)
    
    output:
        tuple val(id), path("${reads.baseName}.g_strand.fastq"), emit: g_strand
        tuple val(id), path("${reads.baseName}.c_strand.fastq"), emit: c_strand

    script:
    """
    separate_strands.py ${reads} ${reads.baseName}.g_strand.fastq ${reads.baseName}.c_strand.fastq
    """

}
process REVERSE_COMPLEMENTATION {
    
    label 'tarpon'

    input:
        tuple val(run_name), path(reads)

    output:
        tuple val(run_name), path("20_80_removed_reads.fastq"), emit: removed_reads
        tuple val(run_name), path("putative_reads.filtered.fastq"), emit: retained_reads
    
    //publishDir "${params.outdir}/FILTERED_READS/", mode: 'copy', overwrite: true, pattern: "20_80_removed*.fastq"
    //publishDir "${params.outdir}/TELOMERIC/", mode: 'copy', overwrite: true, pattern: "putative*.fastq"

    script:
    """
    reverse_complement_reads.py ${reads} ${params.repeat} ${params.reverse_complement_threshold} ${params.c_strand_only} putative_reads.filtered.fastq 20_80_removed_reads.fastq
    """
}

process IDENTIFY_TAGGING_ADAPTOR_AND_DEMUX {
    
    label 'tarpon'

    input:
        tuple val(run_name), path(reads)
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
        identify_tagging_barcodes.py ${reads} ${barcode_errors} ${params.repeat} DEMUX/ adaptor_filtered.fastq
        """
    else if (params.sample_file == "")
        // single sample provided, not multiplexed
        """
        mkdir DEMUX/
        identify_tagging_adaptor.py ${reads} ${params.adaptor_sequence} ${params.adaptor_sequence_errors} ${params.repeat} DEMUX/${params.sample_name}.fastq adaptor_filtered.fastq
        cat DEMUX/*.fastq > adaptor.fastq
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

process SUBTELO_FILTERING {
    label 'tarpon'

    input:
        tuple val(sample), path(reads)
    
    output:
        tuple val(sample), path("*subtelo_fail.fastq"), emit: filtered_reads
        tuple val(sample), path("*subtelo_pass.fastq"), emit: retained_reads

    script:
    """
    filter_by_subtelo.py ${reads} ${params.min_subtelo_length} ${params.subtelo_threshold} ${params.repeat} ${sample}.subtelo_pass.fastq ${sample}.subtelo_fail.fastq

    """
}


process TELO_START_IDENTIFICATION {
    label 'tarpon'

    input:
        tuple val(sample), path(reads)

    output:
        tuple val(sample), path("*telomeric.fastq"), path("*telomeric_stats.txt"), emit: final_telomeric
        tuple val(sample), path("*no_telomere_start.fastq"), emit: no_telo_start
        tuple val(sample), path("*.below_telo_%_threshold.fastq"), emit: below_telo_threshold
        tuple val(sample), path("*telomeric.fastq"), emit: retained_reads

    //publishDir "${params.outdir}/FINAL_TELO/", mode:'copy', overwite:true, pattern: "telomeric_stats.txt", saveAs: { filename -> "${sample}.stats.txt" }
    //publishDir "${params.outdir}/FINAL_TELO/", mode:'copy', overwite:true, pattern: "TELOMERIC/*", saveAs: { filename -> "${sample}.telomeric.fastq" }
    //publishDir "${params.outdir}/", mode: 'copy', overwrite: true, pattern: "FILTERED_READS/*"
    //publishDir "${params.outdir}/", mode: 'copy', overwrite: true, pattern: "TELOMERIC/*"

    publishDir "${params.outdir}/${sample}/telomeres.fastq", mode: 'copy', overwrite:true, pattern:"*.telomeric.fastq"
    publishDir "${params.outdir}/${sample}/telomere_stats.txt", mode: 'copy', overwrite:true, pattern:"telomeric_stats.txt"
    
    script:
    """
    #python script that identified telomere start. Writes out fastq file, stats file, fastq for reads removed because no telo start was found, fastq for reads removed because didnt reach minimum threshold
    mkdir TELOMERIC
    mkdir FILTERED_READS
    identify_telo_start.py ${reads} ${params.repeat} ${params.sliding_window_size} ${params.sliding_window_interval} ${params.upper_threshold} ${params.lower_threshold} ${params.telomeric_repeat_percentage} ${params.consecutive_repeats} ${sample}.telomeric.fastq ${sample}.no_telomere_start.fastq ${sample}.below_telo_%_threshold.fastq ${sample}.telomeric_stats.txt
    """
}

process INDIVIDUAL_READ_PLOTS {

    label 'tarpon'

    input:
        tuple val(sample), path(reads), path(stats)

    output:
        path("*.pdf")

    publishDir "${params.outdir}/${sample}/FIGURES/INDIVIDUAL_READ_PLOTS/", mode:'copy', overwrite: true, pattern: "*.pdf"

    script:
    """
    indiv_read_plots.py ${reads} ${params.repeat} ${stats} ${params.sliding_window_size} ${params.sliding_window_interval}
    """
}

process GENERATE_PLOTS {

    label 'tarpon'

    input:
        tuple val(sample), path(telo_reads), path(telo_stats)
    
    output:
        path("*.pdf")
        path("C_G_COMPARISON/*.pdf"), optional:true
    
    publishDir "${params.outdir}/${sample}/FIGURES/", mode:'copy', overwrite: true, pattern: "*.pdf"
    publishDir "${params.outdir}/${sample}/FIGURES/", mode:'copy', overwrite: true, pattern: "C_G_COMPARISON/*.pdf"

    script:
    """
    Rscript ${baseDir}/bin/basic_plots.R ${telo_stats} ${params.plot_telo_length} ${params.plot_vrr_length} ${params.strand_comparison}
    """
}

process GENERATE_DETAILED_PLOTS {

    label 'tarpon'

    input:
        tuple val(sample), path(telomeric_reads), path(telo_stats, stageAs: "old_telo_stats.txt")
    
    output:
        tuple val(sample), path(telomeric_reads), path("telomeric_stats.txt"), emit: final_telomeric
        path("DETAILED_STATS/*.pdf")
        path("C_G_COMPARISON/*.pdf"), optional: true
    
    publishDir "${params.outdir}/${sample}/FIGURES/", mode:'copy', overwrite: true, pattern: "DETAILED_STATS/*.pdf"
    publishDir "${params.outdir}/${sample}/FIGURES", mode:'copy', overwrite: true, pattern: "C_G_COMPARISON/*.pdf"
    publishDir "${params.outdir}/${sample}/", mode:'copy', overwrite: true, pattern: "telomeric_stats.txt", saveAs: { filename -> "${sample}.stats.txt" }

    script:
    """
    detailed_stats.py ${telomeric_reads} ${params.repeat} ${telo_stats} telomeric_stats.txt
    Rscript ${baseDir}/bin/detailed_plots.R telomeric_stats.txt ${params.plot_telo_length} ${params.plot_vrr_length} ${params.strand_comparison}
    """

}

process SUMMARY_STATS_RUN {

    label 'tarpon'

    input:
        tuple val(id), path(retained, stageAs: "RETAINED/*")
        tuple val(id), path(filtered, stageAs: "FILTERED/*")
        //path(output_dir)
    
    output:
        tuple val(id), path("Retained_Reads.stats.txt"), emit: retained_stats
        tuple val(id), path("Filtered_Reads.stats.txt"), emit: filtered_stats
        path("*.pdf")
       
    publishDir "${params.outdir}/RUN_STATS/", mode:'copy', overwrite:true, pattern:"*stats.txt"
    publishDir "${params.outdir}/RUN_STATS/FIGURES/", mode:'copy', overwrite:true, pattern:"*pdf"


    script:
    """
    seqkit stats -a -N 50,90 -T ${retained} > Retained_Reads.stats.txt
    seqkit stats -a -N 50,90 -T ${filtered} > Filtered_Reads.stats.txt
    Rscript ${baseDir}/bin/summary_stats_plots.R Retained_Reads.stats.txt ${params.strand_comparison} telomeric
    Rscript ${baseDir}/bin/summary_stats_plots.R Filtered_Reads.stats.txt ${params.strand_comparison} filtered
    """
}

process SUMMARY_STATS_SAMPLE {

    label 'tarpon'

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

    input:
        tuple val(sample), path(telo_sequences), path(telo_stats)
    
    output:
        tuple val(sample), path("digest_stats.txt")
    
    publishDir "${params.outdir}/${sample}/digest_stats.txt", mode: 'copy', overwrite: true, pattern: "digest_stats.txt"

    script:
    """
    for seq in \$(echo "${params.restriction_digest_analysis}" | tr "," "\n"); do seqkit grep -s -p \$seq ${telo_sequences} > \$seq.fastq; done
    seqkit stats -a -N 50,90 -T *.fastq > digest_stats.txt
    """
}

process GENERATE_FINAL_REPORT {

    label 'tarpon'

    input:
        path("params.json")
        path("versions.txt")
        path("manifest.json")
        tuple val(run), path(stats_run_retained)
        tuple val(run1), path(stats_run_filtered)
        path(sample_stats_retained)
        path(sample_stats_filtered)

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
                            --command "${workflow.commandLine}" \
                            --run_stats_retained ${stats_run_retained} \
                            --run_stats_filtered ${stats_run_filtered} \
                            --sample_stats_retained ${sample_stats_retained} \
                            --sample_stats_filtered ${sample_stats_filtered}
    """
}

process getParams {

    label "tarpon"

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


process getVersions {

    label "tarpon"

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

    input:
        tuple val(file_type), path(input_files)

    output:
        tuple val(params.run_name), path("${file_type}.fastq"), emit:combined

    script:
    """
    cat ${input_files} > ${file_type}.fastq
    """
}
