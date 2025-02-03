import groovy.json.JsonOutput


process PUTATIVE_ISOLATION {
    // Process that will take a singular fastq.gz input file and isolate all putative telomeric reads
    // based on the frequency of params.repeat (sequence) being present at least params.repeat_count times

    // It will return 3 output tuples containing the putative telomeric reads, non telomeric reads, and the input channel
    
    // Nothing is published to the output directory

    tag "$run_name - Putative Isolation"
    label 'tarpon'

    input:
        tuple val(run_name), path(reads_file, stageAs: "input.bam")
    
    output:
        tuple val(run_name), path ("putative_reads.bam"), emit: putative_reads
        tuple val (run_name), path("non_telomeric.bam"), emit: non_telomeric
        tuple val(run_name), path("input.bam"), emit: input_ch

    script:
    """
    isolate_putative_telomeric_reads.py --input_file ${reads_file} --repeat ${params.repeat} --repeat_count ${params.repeat_count} --c_strand_only ${params.c_strand_only} --out_file putative_reads.bam --non_telo non_telomeric.bam
    """
}

process SEPERATE_STRANDS {

    // Process that is designed to separate the G and C strand telomeric sequences into their respective categories for statistical analysis purposes

    // Takes an input of a fastq file that is unzipped containg all telomeric reads

    // Returns two tuples for G and C strand telomeric sequences respectively

    // Nothing is published to output directory

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

process ISOLATE_POD5_SQUIGGLES {
    label 'tarpon'
    tag "$run_name - Isolating Pod5 Files"

    input:
        tuple val(run_name), path(reads)
        path(pod5_dir)

    output:
        path("filtered.pod5"), emit: pod5_filtered

    script:
    """
    get_read_ids.py --input_file ${reads} --output_file read_ids.txt
    pod5 filter ${pod5_dir}/*.pod5 --output filtered.pod5 --ids read_ids.txt --missing-ok
    """    
}

process REVERSE_COMPLEMENTATION {
    
    //Process that will take all C strand reads and reverse complement them to G strand reads for pipeline simplicity

    // Input is a fastq file of reads that are then filtered based on the prevalance of G and C strand repeats (chimeras removed)

    // Output is two tuples - one is the reads that are chimeric, the second is all putative telomeric reads that are retained

    // Nothing is published to output directory
    label 'tarpon'
    tag "$run_name - Reverse Complementation"
    input:
        tuple val(run_name), path(reads)

    output:
        tuple val(run_name), path("20_80_removed_reads.fastq"), emit: removed_reads
        tuple val(run_name), path("putative_reads.filtered.fastq"), emit: retained_reads
    
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

    // Process that identifies the end of telomeric repeats from the presence of a barcode or an adaptor sequence
    // and further demultiplexes input fastq based on the sample file

    // takes two inputs : tuple of putative and filtered telomeric sequences and a sample file to demultiplex by
    
    // output a channel containing all reads sucesfully demultiplexed, a tuble containing all reads where an adaptor was succesfully identified
    // and a tuple with all reads not containing an adaptor sequence

    // Nothing is published to output directory

    label 'tarpon'
    tag "$run_name - Identify Adaptor and Demultiplexing"

    input:
        tuple val(run_name), path(reads)
        path(barcodes_file)

    output:
        path("DEMUX/*.fastq"), emit: demuxed_reads
        tuple val(run_name), path("adaptor.fastq"), emit: retained_reads
        tuple val(run_name), path("adaptor_filtered.fastq"), emit: filtered_reads


    script:
    // if no adaptor sequence is provided telomeric reads are demultiplexed and the end is identified simultaneously based on the sequences in the sample file
    if (params.adaptor_sequence == "")
        """
        mkdir DEMUX/
        identify_tagging_barcodes.py --input_file ${reads} --sample_file ${barcodes_file} --barcode_errors ${params.barcode_errors} --repeat ${params.repeat} --out_fh DEMUX/ --no_adaptor adaptor_filtered.fastq
        cat DEMUX/* > adaptor.fastq
        """
    // if both an adaptor sequence and sample file are provided the telomeric end is first identified by the adaptor sequence and then downstream demultiplexed using the sample file
    else
        // run adaptor identification and then demux
        """
        mkdir DEMUX/

        identify_adaptor_and_demux.py --input_file ${reads} \
            --adaptor_sequence ${params.adaptor_sequence} \
            --adaptor_errors ${params.adaptor_sequence_errors} \
            --repeat ${params.repeat} \
            --no_adaptor adaptor_filtered.fastq \
            --sample_file ${barcodes_file} \
            --barcode_errors ${params.barcode_errors} \
            --out_prefix DEMUX 

        cat DEMUX/* > adaptor.fastq

        """
}

process IDENTIFY_TAGGING_ADAPTOR {
    
    // Identical to previous process however with no demultiplexing - this function is run when the adaptor sequence is provided but the sample file is not

    // input and output are identical to previous process

    // Nothing is published to output directory
    label 'tarpon'
    tag "$run_name - Identify Adaptor and Demultiplexing"

    input:
        tuple val(run_name), path(reads)

    output:
        path("${params.sample_name}.fastq"), emit: demuxed_reads
        tuple val(run_name), path("adaptor.fastq"), emit: retained_reads
        tuple val(run_name), path("adaptor_filtered.fastq"), emit: filtered_reads

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

    // Process that looks at start of reasd and removes any read that contains greater than params.subtelo_threshold percentage of telomeric repeats in the first params.min_subtelo_length

    // Input - Takes a fastq file of reads and outputs two tuples of failed and passed filtering reads

    // Nothing is published to output directory

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

process BASECALLING {
    
    label 'basecalling'
    tag "Rebascalling filtered Pod5 File"
    label 'gpu'

    input:
        path(pod5_file)

    output:
        tuple val(params.run_name), path("sup_basecalled.fastq.gz")

    script:
    """
    dorado basecaller sup --no-trim --emit-fastq --recursive ${pod5_file} > sup_basecalled.fastq
    gzip sup_basecalled.fastq
    """
}

process TELO_START_IDENTIFICATION {

    // Process that identifies the start of telomeric reads and performs additional filtering based on sequence prior to start - see manuscript for more details

    // Input fastq file of reads

    // Output - telomeric reads, telomeric read stats, reads where no telomere start was identified, reads that failed filtering

    // Final telomeric reads and telomeric stats are published to output directory

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

    // process that will plot telomeric repeat percentage along read for every finalized telomeric sequence

    // Input : fastq reads and telomeric statistics

    // output: one pdf file per read

    // PDF files are published to output directory

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

process BASIC_PLOTS {

    // process that generates R plots using telomeric reads and telo stats dataframe

    // input = telomeric reads and statistics dataframe

    // output = all relevant plots

    // publishes all relevant plots to output directory

    label 'tarpon'
    tag "$sample - Generating Output Plots"

    input:
        tuple val(sample), path(telo_reads), path(telo_stats)
    
    output:
        path("*.pdf")
        path("STRAND_COMPARISON/*.pdf"), optional:true
    
    publishDir "${params.outdir}/${sample}/FIGURES/", mode:'move', overwrite: true, pattern: "*.pdf"
    publishDir "${params.outdir}/${sample}/FIGURES/", mode:'move', overwrite: true, pattern: "C_G_COMPARISON/*.pdf"

    script:
    """
    basic_plots.R ${telo_stats} ${params.plot_telo_length} ${params.plot_vrr_length} ${params.strand_comparison}
    """
}

process GENERATE_DETAILED_PLOTS {

    // exact same as generate_plots but does so to a much greater extent

    label 'tarpon'
    tag "$sample - Generating Detailed Output Plots"

    input:
        tuple val(sample), path(telomeric_reads), path(telo_stats, stageAs: "old_telo_stats.txt")
    
    output:
        tuple val(sample), path(telomeric_reads), path("*telomeric_stats.txt"), emit: final_telomeric
        path("*telomeric_stats.txt"), emit: final_telo_stats
        path("*.pdf")
    
    publishDir "${params.outdir}/${sample}/FIGURES/DETAILED_STATS/", mode:'move', overwrite: true, pattern: "*.pdf"
    publishDir "${params.outdir}/${sample}/", mode:'copy', overwrite: true, pattern: "telomeric_stats.txt", saveAs: { filename -> "${sample}.stats.txt" }

    // not modified R script
    script:
    """
    detailed_stats.py --fastq_file ${telomeric_reads} --repeat ${params.repeat} --stats_in ${telo_stats} --stats_out ${sample}.telomeric_stats.txt
    detailed_plots.R ${sample}.telomeric_stats.txt ${params.plot_telo_length} ${params.plot_vrr_length} ${params.strand_comparison}
    """

}

process SUMMARY_STATS_RUN {

    // process that takes a list of retained and filtered fastq files and runs seqkit stats on the files for easy plotting

    // input: two tuples composed of retained and filtered read locations

    // output: seqkit stats output files and R plots generated from these statistics

    // publishes stats and figures to output directory

    label 'tarpon'
    tag "$id - Collecting Run Summary Statistics"

    input:
        tuple val(id), path(retained, stageAs: "RETAINED/*")
        tuple val(id), path(filtered, stageAs: "FILTERED/*")
    
    output:
        tuple val(id), path("Retained_Reads.stats.txt"), emit: retained_stats
        tuple val(id), path("Filtered_Reads.stats.txt"), emit: filtered_stats
        tuple val(id), path("retained.quality.txt"), emit: retained_quality_stats
        tuple val(id), path("filtered.quality.txt"), emit: filtered_quality_stats
        path("*.pdf")
       
    publishDir "${params.outdir}/RUN_STATS/", mode:'copy', overwrite:true, pattern:"*stats.txt"
    publishDir "${params.outdir}/RUN_STATS/", mode:'copy', overwrite:true, pattern:"*.quality.txt"
    publishDir "${params.outdir}/RUN_STATS/FIGURES/", mode:'copy', overwrite:true, pattern:"*.pdf"
    publishDir "${params.outdir}/RUN_STATS/FIGURES/STRAND_COMPARISON/", mode:'copy', overwrite:true, pattern: "*strand*.pdf"

    script:
    """
    seqkit stats -a -N 50,90 -T ${retained} > Retained_Reads.stats.txt
    seqkit stats -a -N 50,90 -T ${filtered} > Filtered_Reads.stats.txt
    getQualityDistro.py --fastq_files ${retained} > retained.quality.txt
    getQualityDistro.py --fastq_files ${filtered} > filtered.quality.txt
    summary_stats_plots.R Retained_Reads.stats.txt ${params.strand_comparison} telomeric retained.quality.txt
    summary_stats_plots.R Filtered_Reads.stats.txt ${params.strand_comparison} filtered filtered.quality.txt
    """
}

process SUMMARY_STATS_SAMPLE {

    // does the same thing as summary stats run but for each individual demultiplexed sample

    label 'tarpon'
    tag "$id - Collecting Sample Summary Statistics"

    input:
        tuple val(id), path(retained, stageAs: "RETAINED/*")
        tuple val(id), path(filtered, stageAs: "FILTERED/*")
    
    output:
        tuple val(id), path("*retained.stats.txt"), path("*.retained.quality.txt"), emit: retained_stats
        tuple val(id), path("*.filtered.stats.txt"), path("*.filtered.quality.txt"), emit: filtered_stats
    
    publishDir "${params.outdir}/${id}", mode:'copy', overwrite:true, pattern:"*stats.txt"
    publishDir "${params.outdir}/${id}", mode:'copy', overwrite:true, pattern:"*quality.txt"
    
    script:
    """
    seqkit stats -a -N 50,90 -T ${retained} > ${id}.retained.stats.txt
    seqkit stats -a -N 50,90 -T ${filtered} > ${id}.filtered.stats.txt
    getQualityDistro.py --fastq_files ${retained} > ${id}.retained.quality.txt
    getQualityDistro.py --fastq_files ${filtered} > ${id}.filtered.quality.txt
    """
}

process RESTRICTION_DIGEST_ANALYSIS {

    // if the params.restriction_digest_analysis is set will search for restriction sites in the telomeric sequences based on a comma separated list

    // input telomeric sequences and statistics

    // output: digestion stats - one line per restriction site

    // publishes stats to output directory
     
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
        tuple val(run), path(quality_run_retained)
        tuple val(run), path(quality_run_filtered)
        path(sample_stats_retained)
        path(sample_stats_filtered)
        path(sample_quality_retained)
        path(sample_quality_filtered)
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
                            --minimum_read_count ${params.minimum_telo_reads_per_sample} \
                            --commandLine "${workflow.commandLine}" \
                            --run_stats_retained ${stats_run_retained} \
                            --run_stats_filtered ${stats_run_filtered} \
                            --run_quality_retained ${quality_run_retained} \
                            --run_quality_filtered ${quality_run_filtered} \
                            --sample_stats_retained ${sample_stats_retained} \
                            --sample_stats_filtered ${sample_stats_filtered} \
                            --sample_quality_retained ${sample_quality_retained} \
                            --sample_quality_filtered ${sample_quality_filtered} \
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
    stageInMode "copy"

    input:
        tuple val(file_type), path(input_file)
    
    output:
        tuple val(params.run_name), path("${file_type}.fastq.gz"), emit: combined

    script:
    """
    gzip $input_file 
    mv ${input_file}.gz ${file_type}.fastq.gz
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
    samtools fastq -@ 4 $input_file > ${file_type}.fastq
    gzip ${file_type}.fastq
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

process COMBINE_FASTQ_AND_ZIP {
    label 'tarpon'
    tag "$file_type Concatenating and Zipping FASTQ Files"

    input:
        tuple val(file_type), path(input_files)

    output:
        tuple val(params.run_name), path("${file_type}.fastq.gz"), emit:combined

    script:
    """
    cat ${input_files} > ${file_type}.fastq
    gzip ${file_type}.fastq
    """
}

process COMBINE_BAM {
    label 'tarpon'
    tag "$file_type Combining BAM Files"

    input: 
        tuple val(file_type), path(input_files)

    output:
        tuple val(params.run_name), path("${file_type}.bam"), emit:combined

    script:
    """
    samtools merge -o ${file_type}.bam ${input_files} 
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
        path("combined_stats.txt"), emit: stats
        path("combined_stats.VRR.txt"), emit: vrr_stats
        path("*.pdf")

    publishDir "${params.outdir}/", overwrite: true, mode: 'copy', pattern: "sample_stats.txt"
    publishDir "${params.outdir}/", overwrite: true, mode: 'copy', pattern: "sample_stats.VRR.txt"
    publishDir "${params.outdir}/RUN_STATS/FIGURES/", overwrite: true, mode: 'copy', pattern: "*.pdf"
    
    // add another r script that takes input files

    script:
    """
    processTelomereStats.py --stat_files ${input_files} --vrr_length ${params.plot_vrr_length} --telo_length ${params.plot_telo_length}
    sampleComparison_Plots.R ${params.plot_vrr_length} ${params.plot_telo_length}
    sampleComparison_BarPlot.R ${params.plot_telo_length} ${params.plot_vrr_length} ${input_files}
    """

}
