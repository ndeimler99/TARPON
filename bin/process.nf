process putative_isolation {

    label 'seqkit'

    input:
        path(reads_file)
    
    output:
        path("putative_reads.fastq"), emit: putative_reads
        path("raw_data.stats.txt")
    
    publishDir "${params.outdir}/STATS/", mode: 'copy', overwrite: true, pattern: "raw_data.stats.txt"
    
    script:
    """
    seqkit stats -a -N 50,90 -T *.gz > raw_data.stats.txt
    python3 ${baseDir}/bin/isolate_putative_telomeric_reads.py ${reads_file} ${params.repeat} ${params.repeat_count} ${params.c_strand_only} putative_reads.fastq
    """
}

process reverse_complementation {
    
    label 'seqkit'

    input:
        path(reads)

    output:
        path("putative_reads.reverse_complemented.fastq"), emit: reversed_reads
        path("putative_reads.stats.txt")
    
    publishDir "${params.outdir}/STATS/", mode: 'copy', overwrite: true, pattern: "putative_reads.stats.txt"

    script:
    """
    python3 ${baseDir}/bin/reverse_complement_reads.py ${reads} ${params.repeat} ${params.reverse_complement_threshold} ${params.c_strand_only} putative_reads.reverse_complemented.fastq putative_reads.g_strand.fastq putative_reads.c_strand.fastq
    seqkit stats -a -N 50,90 -T *.fastq > putative_reads.stats.txt
    """
}

process identify_tagging_adaptor {
    label 'seqkit'

    input:
        path(reads)
    output:
        path("adaptor_present.reads.stats.txt")
        path("subtelo_filtered_reads.fastq"), emit: adaptor_reads

    publishDir "${params.outdir}/STATS/", mode: 'copy', overwrite: true, pattern: "adaptor_present.reads.stats.txt"
    
    script:
    """
    python3 ${baseDir}/bin/identify_tagging_adaptor.py ${reads} ${params.adaptor_sequence} ${params.adaptor_sequence_errors} reads_with_adaptor.fastq ${params.min_subtelo_length} ${params.subtelo_threshold} subtelo_filtered_reads.fastq subtelo_removed.fastq ${params.repeat}
    seqkit stats -a -N 50,90 -T *.fastq > adaptor_present.reads.stats.txt
    """
}

process consecutive_identification {

    input:
        path(reads)

    output:
        path("telomeric_start_identified.fastq"), emit: telomeric_sequences
        path("telo_read_stats.txt")
    publishDir "${params.outdir}/STATS/", mode: 'copy', overwrite: true, pattern: "telo_read_stats.txt"

    script:
    """
    python3 ${baseDir}/bin/identify_consecutive_telo_start.py ${reads} ${params.repeat} ${params.repeat_count} ${params.telomeric_repeat_percentage} "telomeric_reads.fastq" "non_telomeric_reads.fastq"
    seqkit stats -a -N 50,90 -T *.fastq > telo_read_stats.txt
    """
}

process threshold_identification {

    input:
        path(reads)

    output:
        path("telomeric_start_identified.fastq"), emit: telomeric_sequences
        path("telo_read_stats.txt")
    publishDir "${params.outdir}/STATS/", mode: 'copy', overwrite: true, pattern: "telo_read_stats.txt"

    script:
    """
    python3 ${baseDir}/bin/identify_threshold_telo_start.py ${reads} ${params.repeat} ${params.perfect_repeats} ${params.sliding_window_size} ${params.sliding_window_interval} ${params.start_repeat_threshold} ${params.end_repeat_threshold} ${params.telomeric_repeat_percentage} "telomeric_reads.fastq" "non_telomeric_reads.fastq"
    seqkit stats -a -N 50,90 -T *.fastq > telo_read_stats.txt
    """
}