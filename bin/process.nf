process PUTATIVE_ISOLATION {

    label 'tarpon'

    input:
        tuple val(run_name), path(reads_file, stageAs: "input.fastq.gz")
    
    output:
        tuple val(run_name), path ("putative_reads.fastq"), emit: putative_reads
        tuple val (run_name), path("non_telomeric.fastq"), emit: non_telomeric
    
    //publishDir "${params.outdir}/TELOMERIC/", overwrite: true, mode: 'copy', pattern: "input.fastq.gz"
    //publishDir "${params.outdir}/FILTERED_READS/", overwrite: true, mode: 'copy', pattern: "input.fastq.gz"
    //publishDir "${params.outdir}/TELOMERIC/", mode: 'copy', overwrite: true, pattern: "putative_reads.fastq"
    //publishDir "${params.outdir}/FILTERED_READS/",  overwrite: true, pattern: "non_telomeric.fastq"

    script:
    """
    python3 ${baseDir}/bin/isolate_putative_telomeric_reads.py ${reads_file} ${params.repeat} ${params.repeat_count} ${params.c_strand_only} putative_reads.fastq non_telomeric.fastq
    """
}

process REVERSE_COMPLEMENTATION {
    
    label 'tarpon'

    input:
        tuple val(run_name), path(reads)

    output:
        tuple val(run_name), path("20_80_removed_reads.fastq"), emit: removed_reads
        tuple val(run_name), path("putative_reads.filtered.fastq"), emit: reversed_reads
        tuple val(run_name), path("putative.*strand.fastq"), optional: true, emit: strand_specific
    
    //publishDir "${params.outdir}/FILTERED_READS/", mode: 'copy', overwrite: true, pattern: "20_80_removed*.fastq"
    //publishDir "${params.outdir}/TELOMERIC/", mode: 'copy', overwrite: true, pattern: "putative*.fastq"

    script:
    """
    python3 ${baseDir}/bin/reverse_complement_reads.py ${reads} ${params.repeat} ${params.reverse_complement_threshold} ${params.c_strand_only} putative_reads.filtered.fastq 20_80_removed_reads.fastq
    if ${params.strand_comparison}
    then
        python3 ${baseDir}/bin/separate_strands.py putative_reads.filtered.fastq putative.g_strand.fastq putative.c_strand.fastq
        python3 ${baseDir}/bin/separate_strands.py 20_80_removed_reads.fastq 20_80_removed.g_strand.fastq 20_80_removed.c_strand.fastq
    fi
    """
}

process IDENTIFY_TAGGING_ADAPTOR_AND_DEMUX {
    
    label 'tarpon'

    input:
        tuple val(run_name), path(reads)
    output:
        path("TELOMERIC/DEMUX/*.fastq"), emit: demuxed_reads
        tuple val(run_name), path("TELOMERIC/*"), emit: retained_reads
        tuple val(run_name), path("FILTERED_READS/*"), emit: filtered_reads

    //publishDir "${params.outdir}/", mode: 'copy', overwrite: true, pattern: "FILTERED_READS/*"
    //publishDir "${params.outdir}/", mode: 'copy', overwrite: true, pattern: "TELOMERIC/*"

    script:
    if (params.adaptor_sequence == "")
        // no adaptor sequence provided - using ONT like approach were telomere overhang goes directly into barcode
        // computationally not ideal
        """
        mkdir TELOMERIC
        mkdir TELOMERIC/DEMUX/
        mkdir FILTERED_READS
        """
    else if (params.sample_file == "")
        // single sample provided, not multiplexed
        """
        mkdir TELOMERIC
        mkdir TELOMERIC/DEMUX/
        mkdir FILTERED_READS
        python3 ${baseDir}/bin/identify_tagging_adaptor.py ${reads} ${params.adaptor_sequence} ${params.adaptor_sequence_errors} ${params.min_subtelo_length} ${params.subtelo_threshold} ${params.repeat} TELOMERIC/subtelo.fastq TELOMERIC/DEMUX/${params.sample_name}.fastq FILTERED_READS/subtelo_filtered.fastq FILTERED_READS/adaptor_filtered.fastq
        if ${params.strand_comparison}
        then
            python3 ${baseDir}/bin/separate_strands.py TELOMERIC/subtelo.fastq TELOMERIC/subtelo.g_strand.fastq TELOMERIC/subtelo.c_strand.fastq
            python3 ${baseDir}/bin/separate_strands.py TELOMERIC/${params.sample_name}.fastq TELOMERIC/${params.sample_name}.g_strand.fastq TELOMERIC/DEMUX/${params.sample_name}.c_strand.fastq
        
            python3 ${baseDir}/bin/separate_strands.py FILTERED_READS/subtelo_filtered.fastq FILTERED_READS/subtelo_filtered.g_strand.fastq FILTERED_READS/subtelo_filtered.c_strand.fastq
            python3 ${baseDir}/bin/separate_strands.py FILTERED_READS/adaptor_filtered.fastq FILTERED_READS/adaptor_filtered.g_strand.fastq FILTERED_READS/adaptor_filtered.c_strand.fastq
        fi
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

process TELO_START_IDENTIFICATION {
    label 'tarpon'

    input:
        tuple val(sample), path(reads)

    output:
        tuple val(sample), path("TELOMERIC/telomeric.fastq"), path("telomeric_stats.txt"), emit: telomeric_sequences
        tuple val(sample), path("TELOMERIC/*.fastq"), emit: retained_reads
        tuple val(sample), path("FILTERED_READS/*.fastq"), emit: removed_reads

    publishDir "${params.outdir}/FINAL_TELO/", mode:'copy', overwite:true, pattern: "telomeric_stats.txt", saveAs: { filename -> "${sample}.stats.txt" }
    publishDir "${params.outdir}/FINAL_TELO/", mode:'copy', overwite:true, pattern: "TELOMERIC/*", saveAs: { filename -> "${sample}.telomeric.fastq" }
    //publishDir "${params.outdir}/", mode: 'copy', overwrite: true, pattern: "FILTERED_READS/*"
    //publishDir "${params.outdir}/", mode: 'copy', overwrite: true, pattern: "TELOMERIC/*"

    script:
    """
    #python script that identified telomere start. Writes out fastq file, stats file, fastq for reads removed because no telo start was found, fastq for reads removed because didnt reach minimum threshold
    mkdir TELOMERIC
    mkdir FILTERED_READS
    python3 ${baseDir}/bin/identify_telo_start.py ${reads} ${params.repeat} ${params.sliding_window_size} ${params.sliding_window_interval} ${params.upper_threshold} ${params.lower_threshold} ${params.telomeric_repeat_percentage} ${params.consecutive_repeats} TELOMERIC/telomeric.fastq FILTERED_READS/no_telomere_start.fastq FILTERED_READS/below_telo_%_threshold.fastq telomeric_stats.txt

    if ${params.strand_comparison}
    then
        python3 ${baseDir}/bin/separate_strands.py TELOMERIC/telomeric.fastq TELOMERIC/telomeric.g_strand.fastq TELOMERIC/telomeric.c_strand.fastq
        python3 ${baseDir}/bin/separate_strands.py FILTERED_READS/no_telomere_start.fastq FILTERED_READS/no_telomere_start.g_strand.fastq FILTERED_READS/no_telomere_start.c_strand.fastq  
        python3 ${baseDir}/bin/separate_strands.py FILTERED_READS/below_telo_%_threshold.fastq FILTERED_READS/below_telo_%_threshold.g_strand.fastq FILTERED_READS/below_telo_%_threshold.c_strand.fastq
    fi
    """
}

process INDIVIDUAL_READ_PLOTS {

    label 'tarpon'

    input:
        tuple val(sample), path(reads), path(stats)

    output:
        path("*.pdf")

    publishDir "${params.outdir}/FIGURES/INDIVIDUAL_READ_PLOTS/${sample}/", mode:'copy', overwrite: true, pattern: "*.pdf"

    script:
    """
    python3 ${baseDir}/bin/indiv_read_plots.py ${reads} ${params.repeat} ${telo_stats} ${params.sliding_window_size} ${params.sliding_window_interval}
    """
}

process GENERATE_PLOTS {

    label 'tarpon'

    input:
        tuple val(sample), path(telo_reads), path(telo_stats)
    
    output:
        path("*.pdf")
        path("C_G_COMPARISON/*.pdf"), optional:true
    
    publishDir "${params.outdir}/FIGURES/${sample}", mode:'copy', overwrite: true, pattern: "*.pdf"
    publishDir "${params.outdir}/FIGURES/${sample}", mode:'copy', overwrite: true, pattern: "C_G_COMPARISON/*.pdf"

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
        tuple val(sample), path(telomeric_reads), path("telomeric_stats.txt"), emit: telomeric_sequences
        path("DETAILED_STATS/*.pdf")
        path("C_G_COMPARISON/*.pdf")
    
    publishDir "${params.outdir}/FIGURES/${sample}", mode:'copy', overwrite: true, pattern: "DETAILED_STATS/*.pdf"
    publishDir "${params.outdir}/FIGURES/${sample}", mode:'copy', overwrite: true, pattern: "C_G_COMPARISON/*.pdf"
    publishDir "${params.outdir}/FINAL_TELO/", mode:'copy', overwrite: true, pattern: "telomeric_stats.txt", saveAs: { filename -> "${sample}.stats.txt" }

    script:
    """
    python3 ${baseDir}/bin/detailed_stats.py ${telomeric_reads} ${params.repeat} ${telo_stats} telomeric_stats.txt
    Rscript ${baseDir}/bin/detailed_plots.R telomeric_stats.txt ${params.plot_telo_length} ${params.plot_vrr_length} ${params.strand_comparison}
    """

}


process SUMMARY_STATS {

    label 'tarpon'

    input:
        path(telo_stats), stageAs: "telo_stats.txt"
        path(output_dir)
    
    output:
        path("telomeric_stats.txt")
        path("filtered_stats.txt")
        path("*.pdf")
    
    publishDir "${params.outdir}/FIGURES/", mode:'copy', overwrite: true, pattern: "*.pdf"
    publishDir "${params.outdir}/STATS/", mode: 'copy', overwrite: true, pattern: "*.txt"

    script:
    """
    seqkit stats -a -N 50,90 -T ${output_dir}/TELOMERIC/*.fastq* > telomeric_stats.txt
    seqkit stats -a -N 50,90 -T ${output_dir}/FILTERED_READS/*.fastq* > filtered_stats.txt

    Rscript ${baseDir}/bin/seqkit_stats_plots.R telomeric_stats.txt ${params.strand_comparison} telomeric
    Rscript ${baseDir}/bin/seqkit_stats_plots.R filtered_stats.txt ${params.strand_comparison} filtered
    """
}


process RESTRICTION_DIGEST_ANALYSIS {

    label 'tarpon'

    input:
        path(telo_sequences)
    
    output:
        path("digest_stats.txt")
    
    publishDir "${params.outdir}/STATS/", mode: 'copy', overwrite: true, pattern: "*.txt"

    script:
    """
    for seq in \$(echo "${params.restriction_digest_analysis}" | tr "," "\n"); do seqkit grep -s -p \$seq ${telo_sequences} > \$seq.fastq; done
    seqkit stats -a -N 50,90 -T *.fastq > digest_stats.txt
    """
}

process GENERATE_HTML_REPORT {

    label 'tarpon'

    input:
        path(output_dir)
        path(stats_file)
    
    output:
        path("report.html")
    
    publishDir "${params.outdir}/", mode: 'copy', overwrite: true, pattern: "report.html"

    script:
    """
    python3 ${baseDir}/bin/generate_report.py ${output_dir} 0 0 ${params.input_file} ${params.repeat} ${params.outdir} ${baseDir}/bin/single_sample_template.html report.html ${baseDir}/bin/report.css
    """
}
