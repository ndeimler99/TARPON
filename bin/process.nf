process putative_isolation {

    label 'tarpon'

    input:
        path(reads_file), stageAs: "input.fastq.gz"
    
    output:
        path("input.fastq.gz")
        path("putative_reads.fastq"), emit: putative_reads
        path("non_telomeric.fastq")
    
    publishDir "${params.outdir}/TELOMERIC/", overwrite: true, mode: 'copy', pattern: "input.fastq.gz"
    publishDir "${params.outdir}/FILTERED_READS/", overwrite: true, mode: 'copy', pattern: "input.fastq.gz"
    publishDir "${params.outdir}/TELOMERIC/", mode: 'copy', overwrite: true, pattern: "putative_reads.fastq"
    publishDir "${params.outdir}/FILTERED_READS/",  overwrite: true, pattern: "non_telomeric.fastq"

    script:
    """
    python3 ${baseDir}/bin/isolate_putative_telomeric_reads.py ${reads_file} ${params.repeat} ${params.repeat_count} ${params.c_strand_only} putative_reads.fastq non_telomeric.fastq
    """
}

process reverse_complementation {
    
    label 'tarpon'

    input:
        path(reads)

    output:
        path("20_80_removed_reads.fastq")
        path("putative_reads.c_g_filtered.fastq"), emit: reversed_reads
        path("putative.g_strand.fastq"), optional: true
        path("putative.c_strand.fastq"), optional: true
    
    publishDir "${params.outdir}/FILTERED_READS/", mode: 'copy', overwrite: true, pattern: "20_80_removed*.fastq"
    publishDir "${params.outdir}/TELOMERIC/", mode: 'copy', overwrite: true, pattern: "putative*.fastq"

    script:
    """
    python3 ${baseDir}/bin/reverse_complement_reads.py ${reads} ${params.repeat} ${params.reverse_complement_threshold} ${params.c_strand_only} putative_reads.c_g_filtered.fastq 20_80_removed_reads.fastq
    if ${params.strand_comparison}
    then
        python3 ${baseDir}/bin/separate_strands.py putative_reads.c_g_filtered.fastq putative.g_strand.fastq putative.c_strand.fastq
        python3 ${baseDir}/bin/separate_strands.py 20_80_removed_reads.fastq 20_80_removed.g_strand.fastq 20_80_removed.c_strand.fastq
    fi
    """
}

process identify_tagging_adaptor {
    
    label 'tarpon'

    input:
        path(reads)
    output:
        path("FILTERED_READS/*.fastq")
        path("TELOMERIC/*.fastq")
        path("TELOMERIC/adaptor.fastq"), emit: adaptor_reads

    publishDir "${params.outdir}/", mode: 'copy', overwrite: true, pattern: "FILTERED_READS/*"
    publishDir "${params.outdir}/", mode: 'copy', overwrite: true, pattern: "TELOMERIC/*"

    script:
    """
    mkdir TELOMERIC
    mkdir FILTERED_READS
    python3 ${baseDir}/bin/identify_tagging_adaptor.py ${reads} ${params.adaptor_sequence} ${params.adaptor_sequence_errors} ${params.min_subtelo_length} ${params.subtelo_threshold} ${params.repeat} TELOMERIC/subtelo.fastq TELOMERIC/adaptor.fastq FILTERED_READS/subtelo_filtered.fastq FILTERED_READS/adaptor_filtered.fastq
    if ${params.strand_comparison}
    then
        python3 ${baseDir}/bin/separate_strands.py TELOMERIC/subtelo.fastq TELOMERIC/subtelo.g_strand.fastq TELOMERIC/subtelo.c_strand.fastq
        python3 ${baseDir}/bin/separate_strands.py TELOMERIC/adaptor.fastq TELOMERIC/adaptor.g_strand.fastq TELOMERIC/adaptor.c_strand.fastq
        
        python3 ${baseDir}/bin/separate_strands.py FILTERED_READS/subtelo_filtered.fastq FILTERED_READS/subtelo_filtered.g_strand.fastq FILTERED_READS/subtelo_filtered.c_strand.fastq
        python3 ${baseDir}/bin/separate_strands.py FILTERED_READS/adaptor_filtered.fastq FILTERED_READS/adaptor_filtered.g_strand.fastq FILTERED_READS/adaptor_filtered.c_strand.fastq
    fi
    """
}

process telo_start_identification {
    label 'tarpon'
    input:
        path(reads)

    output:
        path("TELOMERIC/telomeric.fastq"), emit: telomeric_sequences
        path("TELOMERIC/*.fastq")
        path("FILTERED_READS/*.fastq")
        path("telomeric_stats.txt"), emit: telomeric_stats

    publishDir "${params.outdir}/FINAL_TELO/", mode:'copy', overwite:true, pattern: "telomeric_stats.txt"
    publishDir "${params.outdir}/FINAL_TELO/", mode:'copy', overwite:true, pattern: "TELOMERIC/*", saveAs: { filename -> "telomeric.fastq" }
    publishDir "${params.outdir}/", mode: 'copy', overwrite: true, pattern: "FILTERED_READS/*"
    publishDir "${params.outdir}/", mode: 'copy', overwrite: true, pattern: "TELOMERIC/*"

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

process individual_read_plots {

    label 'tarpon'

    input:
        path(reads)
        path(telo_stats)
    output:
        path("*.pdf")

    publishDir "${params.outdir}/FIGURES/INDIVIDUAL_READ_PLOTS/", mode:'copy', overwrite: true, pattern: "*.pdf"

    script:
    """
    python3 ${baseDir}/bin/indiv_read_plots.py ${reads} ${params.repeat} ${telo_stats} ${params.sliding_window_size} ${params.sliding_window_interval}
    """
}

process generate_plots {

    label 'tarpon'

    input:
        path(telo_stats)
    
    output:
        path("*.pdf")
        path("C_G_COMPARISON/*.pdf"), optional:true
    
    publishDir "${params.outdir}/FIGURES/", mode:'copy', overwrite: true, pattern: "*.pdf"
    publishDir "${params.outdir}/FIGURES/", mode:'copy', overwrite: true, pattern: "C_G_COMPARISON/*.pdf"

    script:
    """
    Rscript ${baseDir}/bin/basic_plots.R ${telo_stats} ${params.plot_telo_length} ${params.plot_vrr_length} ${params.strand_comparison}
    """
}

process generate_detailed_plots {

    label 'tarpon'

    input:
        path(telo_stats), stageAs: "old_telo_stats.txt"
        path(telomeric_reads)
    
    output:
        path("telomeric.fastq"), emit: telomeric_sequences
        path("telomeric_stats.txt"), emit: telomeric_stats
        path("DETAILED_STATS/*.pdf")
        path("C_G_COMPARISON/*.pdf")
    
    publishDir "${params.outdir}/FIGURES/", mode:'copy', overwrite: true, pattern: "DETAILED_STATS/*.pdf"
    publishDir "${params.outdir}/FIGURES/", mode:'copy', overwrite: true, pattern: "C_G_COMPARISON/*.pdf"
    publishDir "${params.outdir}/FINAL_TELO/", mode:'copy', overwrite: true, pattern: "telomeric_stats.txt"

    script:
    """
    python3 ${baseDir}/bin/detailed_stats.py ${telomeric_reads} ${params.repeat} ${telo_stats} telomeric_stats.txt
    Rscript ${baseDir}/bin/detailed_plots.R telomeric_stats.txt ${params.plot_telo_length} ${params.plot_vrr_length} ${params.strand_comparison}
    """

}


process summary_stats {

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


process restriction_digest_analysis {

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

process cleanUp {

    script:
    """
    rm -rf ${baseDir}/work
    """
}

process generate_html_report {

    label 'tarpon'

    input:
        path(output_dir)
    
    output:
        path("report.html")
    
    publishDir "${params.outdir}/", mode: 'copy', overwrite: true, pattern: "report.html"

    script:
    """
    python3 ${baseDir}/bin/generate_report.py ${output_dir} 0 0 ${params.input_file} ${params.repeat} ${params.outdir} ${baseDir}/bin/single_sample_template.html report.html
    """
}
