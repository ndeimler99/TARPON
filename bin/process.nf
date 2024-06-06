process putative_isolation {

    label 'seqkit'

    input:
        path(reads_file)
    
    output:
        path("putative_reads.fastq"), emit: putative_reads
    
    publishDir "${params.outdir}/STATS/", mode: 'copy', overwrite: true, pattern: "raw_data.stats.txt"
    
    script:
    """
    seqkit stats -a -N 50,90 -T *.gz > raw_data.stats.txt
    python3 ${baseDir}/bin/isolate_putative_telomeric_reads.py ${reads_file} ${params.repeat} ${params.repeat_count} ${params.c_strand_only} putative_reads.fastq
    """
}

process rev_comp {
    
    label 'seqkit'

    input:
        path(reads)

    output:
        path("reverse_reads.fastq.gz"), emit: reversed_reads
        path("done.txt")
    
    script:
    """
    seqkit seq -r -p ${reads} > reverse_reads.fastq
    gzip reverse_reads.fastq
    touch done.txt
    """
}

