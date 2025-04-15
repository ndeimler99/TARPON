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
    isolate_putative_telomeric_reads.py --input_file ${reads_file} --repeat ${params.repeat} --repeat_count ${params.repeat_count} --c_strand_only ${params.c_strand_only} --out_file putative_reads.bam --non_telo non_telomeric.bam --mutant ${params.mutant}
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
        tuple val(id), path("${reads.baseName}.g_strand.bam"), emit: g_strand
        tuple val(id), path("${reads.baseName}.c_strand.bam"), emit: c_strand

    script:
    """
    separate_strands.py --input_file ${reads} --g_file ${reads.baseName}.g_strand.bam --c_file ${reads.baseName}.c_strand.bam
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
        tuple val(run_name), path("20_80_removed_reads.bam"), emit: removed_reads
        tuple val(run_name), path("putative_reads.filtered.bam"), emit: retained_reads
    
    script:
    """
    reverse_complement_reads.py --input_file ${reads} \
        --repeat ${params.repeat} \
        --threshold ${params.reverse_complement_threshold} \
        --c_strand_only ${params.c_strand_only} \
        --out_file putative_reads.filtered.bam \
        --removed_reads 20_80_removed_reads.bam \
        --mutant ${params.mutant}
    """
}

process IDENTIFY_TAGGING_CAPTURE_PROBE_AND_DEMUX {

    // Process that identifies the end of telomeric repeats from the presence of a barcode or an adaptor sequence
    // and further demultiplexes input fastq based on the sample file

    // takes two inputs : tuple of putative and filtered telomeric sequences and a sample file to demultiplex by
    
    // output a channel containing all reads sucesfully demultiplexed, a tuble containing all reads where an adaptor was succesfully identified
    // and a tuple with all reads not containing an adaptor sequence

    // Nothing is published to output directory

    label 'tarpon'
    tag "$run_name - Identify Capture Probe and Demultiplexing"

    input:
        tuple val(run_name), path(reads)
        path(barcodes_file)

    output:
        path("DEMUX/*.bam"), emit: demuxed_reads
        tuple val(run_name), path("adaptor.bam"), emit: retained_reads
        tuple val(run_name), path("adaptor_removed.bam"), emit: removed_reads


    script:
    // if no adaptor sequence is provided telomeric reads are demultiplexed and the end is identified simultaneously based on the sequences in the sample file
    if (params.capture_probe_sequence == "")
        """
        mkdir DEMUX/
        identify_tagging_barcodes.py --input_file ${reads} --sample_file ${barcodes_file} \
             --barcode_errors ${params.barcode_errors} --repeat ${params.repeat} \
              --out_fh DEMUX/ --no_adaptor adaptor_removed.bam \
              --mutant ${params.mutant} \
              --overhang_length ${params.capture_probe_overhang_length}
        samtools merge -o adaptor.bam DEMUX/*.bam
        """
    // if both an adaptor sequence and sample file are provided the telomeric end is first identified by the adaptor sequence and then downstream demultiplexed using the sample file
    else
        // run adaptor identification and then demux
        """
        mkdir DEMUX/

        identify_adaptor_and_demux.py --input_file ${reads} \
            --adaptor_sequence ${params.capture_probe_sequence} \
            --adaptor_errors ${params.capture_probe_sequence_errors} \
            --repeat ${params.repeat} \
            --no_adaptor adaptor_removed.bam \
            --sample_file ${barcodes_file} \
            --barcode_errors ${params.barcode_errors} \
            --out_prefix DEMUX \
            --mutant ${params.mutant} \
            --overhang_length ${params.capture_probe_overhang_length}

        samtools merge -o adaptor.bam DEMUX/*.bam

        """
}

process IDENTIFY_TAGGING_CAPTURE_PROBE {
    
    // Identical to previous process however with no demultiplexing - this function is run when the adaptor sequence is provided but the sample file is not

    // input and output are identical to previous process

    // Nothing is published to output directory
    label 'tarpon'
    tag "$run_name - Identify Capture Probe and Demultiplexing"

    input:
        tuple val(run_name), path(reads)

    output:
        path("${params.sample_name}.bam"), emit: demuxed_reads
        tuple val(run_name), path("adaptor.bam"), emit: retained_reads
        tuple val(run_name), path("adaptor_removed.bam"), emit: removed_reads

    script:
    """
    mkdir DEMUX/
    identify_tagging_adaptor.py --input_file ${reads} \
        --adaptor_sequence ${params.capture_probe_sequence} \
        --adaptor_errors ${params.capture_probe_sequence_errors} \
        --repeat ${params.repeat} \
        --adaptor_found ${params.sample_name}.bam \
        --no_adaptor adaptor_removed.bam \
        --mutant ${params.mutant} \
        --overhang_length ${params.capture_probe_overhang_length}

    cp ${params.sample_name}.bam adaptor.bam
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
        tuple val(sample), path("*subtelo_fail.bam"), emit: removed_reads
        tuple val(sample), path("*subtelo_pass.bam"), emit: retained_reads

    script:
    """
    filter_by_subtelo.py --input_file ${reads} \
        --min_subtelo_length ${params.min_subtelo_length} \
        --min_subtelo_threshold ${params.subtelo_threshold} \
        --repeat ${params.repeat} \
        --passes_subtelo ${sample}.subtelo_pass.bam \
        --fails_subtelo ${sample}.subtelo_fail.bam \
        --mutant ${params.mutant}

    """
}

process BASECALLING {
    
    label 'basecalling'
    tag "Rebascalling filtered Pod5 File"
    label 'gpu'

    input:
        path(pod5_file)

    output:
        tuple val(params.run_name), path("sup_basecalled.bam")

    script:
    """
    dorado basecaller sup --no-trim --recursive ${pod5_file} > sup_basecalled.bam
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
        tuple val(sample), path("*telomeric.bam"), path("*telomeric_stats.txt"), emit: final_telomeric
        path("*telomeric_stats.txt"), emit: final_telo_stats
        tuple val(sample), path("*no_telomere_start.bam"), emit: no_telo_start
        tuple val(sample), path("*.below_telo_%_threshold.bam"), emit: below_telo_threshold
        tuple val(sample), path("*telomeric.bam"), emit: retained_reads

    publishDir "${params.outdir}/${sample}/", mode: 'copy', overwrite:true, pattern:"*.telomeric.bam"
    publishDir "${params.outdir}/${sample}/", mode: 'copy', overwrite:true, pattern:"telomeric_stats.txt"
    
    script:
    """
    #python script that identified telomere start. Writes out fastq file, stats file, fastq for reads removed because no telo start was found, fastq for reads removed because didnt reach minimum threshold
    mkdir TELOMERIC
    mkdir REMOVED_READS
    identify_telo_start.py --input_file ${reads} --repeat ${params.repeat} --sliding_window ${params.sliding_window_size} \
        --sliding_window_interval ${params.sliding_window_interval} \
        --upper_threshold ${params.upper_threshold} \
        --lower_threshold ${params.lower_threshold} \
        --telomeric_rep_perc ${params.telomeric_repeat_percentage} \
        --consecutive_repeats ${params.consecutive_repeats} \
        --consecutive_threshold ${params.consecutive_threshold} \
        --telomeric_fastq_out ${sample}.telomeric.bam \
        --no_telomere_out ${sample}.no_telomere_start.bam \
        --filtered_out ${sample}.below_telo_%_threshold.bam \
        --stats_fh ${sample}.telomeric_stats.txt \
        --mutant ${params.mutant} \
        --pre_telomeric_repeat_percentage ${params.pretelomeric_repeat_percentage} \
        --pre_telo_distance ${params.pretelo_start} \
        --minimum_telomere_length ${params.minimum_telomere_length}
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
    basic_plots.R ${telo_stats} ${params.strand_comparison}
    """
}

process GENERATE_DETAILED_PLOTS {

    // exact same as generate_plots but does so to a much greater extent

    label 'tarpon'
    tag "$sample - Generating Detailed Output Plots"

    input:
        tuple val(sample), path(telo_stats)
    
    output:
        path("*.pdf")
    
    publishDir "${params.outdir}/${sample}/FIGURES/DETAILED_STATS/", mode:'move', overwrite: true, pattern: "*.pdf"
    publishDir "${params.outdir}/${sample}/", mode:'copy', overwrite: true, pattern: "telomeric_stats.txt", saveAs: { filename -> "${sample}.stats.txt" }

    // not modified R script
    script:
    """
    detailed_plots.R ${telo_stats} ${params.plot_telo_length} ${params.plot_vrr_length} ${params.strand_comparison}
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
        tuple val(id), path(removed, stageAs: "REMOVED/*")
    
    output:
        tuple val(id), path("Retained_Reads.stats.txt"), emit: retained_stats
        tuple val(id), path("Removed_Reads.stats.txt"), emit: removed_stats
        path("*.pdf")
       
    publishDir "${params.outdir}/RUN_STATS/", mode:'copy', overwrite:true, pattern:"*stats.txt"
    publishDir "${params.outdir}/RUN_STATS/FIGURES/", mode:'copy', overwrite:true, pattern:"*.pdf"
    publishDir "${params.outdir}/RUN_STATS/FIGURES/STRAND_COMPARISON/", mode:'copy', overwrite:true, pattern: "STRAND_COMPARISON/*strand*.pdf"

    script:
    """
    bam_stats.py --bam_files ${retained} --out_file Retained_Reads.stats.txt
    bam_stats.py --bam_files ${removed} --out_file Removed_Reads.stats.txt

    summary_stats_plots.R Retained_Reads.stats.txt ${params.strand_comparison} telomeric
    summary_stats_plots.R Removed_Reads.stats.txt ${params.strand_comparison} removed
    """
}

process SUMMARY_STATS_SAMPLE {

    // does the same thing as summary stats run but for each individual demultiplexed sample

    label 'tarpon'
    tag "$id - Collecting Sample Summary Statistics"

    input:
        tuple val(id), path(retained, stageAs: "RETAINED/*")
        tuple val(id), path(removed, stageAs: "REMOVED/*")
    
    output:
        tuple val(id), path("*retained.stats.txt"), emit: retained_stats
        tuple val(id), path("*.removed.stats.txt"), emit: removed_stats
    
    publishDir "${params.outdir}/${id}", mode:'copy', overwrite:true, pattern:"*stats.txt"
    
    script:
    """
    bam_stats.py --bam_files ${retained} --out_file ${id}.retained.stats.txt
    bam_stats.py --bam_files ${removed} --out_file ${id}.removed.stats.txt
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
        path("*digest_stats.txt"), emit: stats
    
    publishDir "${params.outdir}/${sample}/", mode: 'copy', overwrite: true, pattern: "*digest_stats.txt"

    script:
    """
    for seq in \$(echo "${params.restriction_digest_analysis}" | tr "," "\n"); do samtools view ${telo_sequences} | grep \$seq | samtools view -b > \$seq.bam; done
    for seq in \$(echo "${params.restriction_digest_analysis}" | tr "," "\n"); do samtools reheader ${telo_sequences} \$seq.bam > \$seq.reheaded.bam; done
    bam_stats.py --bam_files *.reheaded.bam --out_file ${sample}.digest_stats.txt
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
        tuple val(run1), path(stats_run_removed)
        path(sample_stats_retained)
        path(sample_stats_removed)
        path(telo_stats_per_sample)
        path(vrr_descriptive_stats)
        path(restriction_digest)
        path(mutant_analysis_repeat_distribution)
        path(mutant_analysis_processivity)

    output:
        path("report.html")
    
    publishDir "${params.outdir}/", mode: 'copy', overwrite: true, pattern: "report.html"
    publishDir "${params.out_dir}/", mode: 'copy', overwrite: true, pattern: "report.html"

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
                            --run_stats_removed ${stats_run_removed} \
                            --sample_stats_retained ${sample_stats_retained} \
                            --sample_stats_removed ${sample_stats_removed} \
                            --sample_telo_stats ${telo_stats_per_sample} \
                            --run_vrr_stats ${vrr_descriptive_stats} \
                            --restriction_digest ${restriction_digest} \
                            --strand_comparison ${params.strand_comparison} \
                            --detailed_stats ${params.detailed_stats} \
                            --mutant ${params.mutant} \
                            --mutant_analysis_repeat_distribution ${mutant_analysis_repeat_distribution} \
                            --mutant_analysis_processivity ${mutant_analysis_processivity} \
                            --repeat ${params.repeat}
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

process COMBINE_FASTQ_GZ {
    
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
        path("combined_stats.VRR.txt"), emit: vrr_stats
        path("*.pdf")

    publishDir "${params.outdir}/", overwrite: true, mode: 'copy', pattern: "sample_stats.txt"
    publishDir "${params.outdir}/", overwrite: true, mode: 'copy', pattern: "sample_stats.VRR.txt"
    publishDir "${params.outdir}/RUN_STATS/FIGURES/", overwrite: true, mode: 'copy', pattern: "*.pdf"
    
    // add another r script that takes input files

    script:
    """
    processTelomereStats.py --stat_files ${input_files}
    sampleComparison_Plots.R
    sampleComparison_BarPlot.R ${input_files}
    """

}

process FASTQ_TO_BAM {
    label 'tarpon'
    tag "$params.run_name Converting FASTQ to BAM"

    input:
         tuple val(run_name), path(fastq_file)

    output:
        tuple val(run_name), path("input.bam")

    script:
    """
    picard -Xmx100G FastqToSam FASTQ=${fastq_file} OUTPUT=input.bam SAMPLE_NAME=input
    """
}

process MUTANT_ANALYSIS {

    label 'tarpon'
    tag "$sample - Mutant Repeat Analysis"

    input:
        tuple val(sample), path(reads), path(stats)

    output:
        tuple val(sample), path("*telo_stats.txt"), emit: statistics
        tuple val(sample), path("telo_sequences.txt"), emit: sequences
        tuple val(sample), path("*processivity_stats.txt"), emit: processivity
        tuple val(sample), path("*repeat_distribution.txt"), emit: repeat_distribution
        path("*.pdf")

    // publish png images
    // publish processivty sand repeat distribution stats
    publishDir "${params.outdir}/${sample}/FIGURES/", overwrite: true, mode: "copy", pattern: "*.pdf"
    publishDir "${params.outdir}/${sample}/", overwrite: true, mode: "copy", pattern: "*.telo_stats.txt"

    script:
    """
    mutantRepeatAnalysis.py --input_file ${reads} --stats_file ${stats} \
                            --wt_processivity ${sample}.wt.processivity_stats.txt \
                            --mt_processivity ${sample}.mt.processivity_stats.txt \
                            --stats_out ${sample}.telo_stats.txt \
                            --repeat_distribution ${sample}.repeat_distribution.txt --repeat ${params.repeat} \
                            --mutant ${params.mutant}
    variantRepeatPlots.R ${sample}.telo_stats.txt ${sample}.repeat_distribution.txt ${params.repeat} ${params.mutant}
    processivityPlots.R ${sample}.wt.processivity_stats.txt ${sample}.mt.processivity_stats.txt ${params.repeat} ${params.mutant}
    """
}

process VARIANT_ANALYSIS {

    label 'tarpon'
    tag "$sample - Mutant Repeat Analysis"

    input:
        tuple val(sample), path(reads), path(stats)

    output:
        tuple val(sample), path("*telo_stats.txt"), emit: statistics
        tuple val(sample), path("telo_sequences.txt"), emit: sequences
        tuple val(sample), path("*processivity.txt"), emit: processivity
        tuple val(sample), path("*repeat_distribution.txt"), emit: repeat_distribution
        //path("*.png")
        path("*.pdf")

    publishDir "${params.outdir}/${sample}/FIGURES/", overwrite: true, mode: "copy", pattern: "*.pdf"
    publishDir "${params.outdir}/${sample}/", overwrite: true, mode: "copy", pattern: "*.telo_stats.txt"

    script:
    """
    variantRepeatAnalysis.py --input_file ${reads} --stats_file ${stats} \
                                --repeat ${params.repeat} \
                                --stats_out ${sample}.telo_stats.txt \
                                --repeat_distribution ${sample}.repeat_distribution.txt

    variantRepeatPlots.R ${sample}.telo_stats.txt ${sample}.repeat_distribution.txt ${params.repeat} ${params.mutant}
    touch ${sample}.processivity.txt
    """
}


process PLOT_TELO_GRAPHS {
    label 'pycairo'

    input:
        tuple val(sample), path(telo_reads)
    
    output:
       path("*.png")

    publishDir "${params.outdir}/${sample}/FIGURES/", overwrite: true, mode: "copy", pattern: "*.png"
    
    script:
    """
    plotTeloGraphs.py --telo_sequences ${telo_reads} --repeat ${params.repeat} --mutant ${params.mutant} --telo_plot ${sample}.telomere_visualization.png
    """
}