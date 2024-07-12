# TArPON
Telomere Analysis Pipeline On Nanopore Sequencing Data


## What is TArPON and Why Should I Use It?

TArPON is a pipeline developed in the lab of Dr. Peter Baumann by PhD candidate Nathaniel Deimler. At the time of the pipeline's creation there were multiple publications (cite publications) that have developed protocols for the enrichment of telomeric sequences to succesfully sequence human telomeres by Nanopore sequencing while remaining cost productive. However, all of these papers focused on the biological aspects of telomeres and the mechanical or chemical enrichment of telomeres for sequencing, not on the computational analysis of individual telomere reads - something never before capable. Unfortunately, all of these papers use different methods to calculate and identify telomeric reads and provide little to no justificiation on why those parameters and techniques were chosen. TArPON aims to provide a one stop shop pipeline for clinicians as well as researchers interested in telomere dynamics using a set of best practices that have been justified in Deimler et al. TArPON is capable of being run in a command line environment or through the ONT GUI Epi2Me and is compatible with all pre-existing methods of telomere sequencing and enrichment including the ONT telo-seq methodology.

## The Pipeline

INSERT SCHEMATIC WORKFLOW HERE

TArPON begins by identifying a set of putative telomeric sequences to reduce computational time and resources at later steps in the pipeline. By default, this is done by isolating reads that contain at least ten perfect telomeric repeats (see parameters: repeat, repeat_count ). These repeats do not have to be consecutive.  The ratio of G telomeric repeats to C telomeric repeats is calculated for all putative telomeric reads. If the ratio falls between 20 and 80% ( see parameter reverse_complement_threshold ) the reads are discarded as they are most likely chimeric sequences.  All reads that are less than 20% G telomeric repeats (these would be C strand reads) are reverse complemented so all telomeric reads are in the same orientation for further analysis. Strandedness information is retaiend throughout the pipeline and strands can be compared using the parameter --strand_comparison.  Telomeric reads are then demultiplexed and the end of the read is identified using either a provided adaptor sequence (--adaptor_sequence) or by the barcodes found in --sample_file. If no sample_file is provided all reads are considered to end in the adaptor_sequence and belong to the same sample: in this case --sample_name should also be provided.  In the case that all telomeric capture probes/sequences are unique and do not contain a shared sequence, --adaptor_sequence should be left blank and the end of the telomeric read will be determined by the barcodes present in sample_file. If both adaptor sequence and sample_file are provided, the end of the telomere will be determined by the presence of the adaptor_sequence and demultiplexed by the barcodes found in sample_file. In all cases of sequence searches, fuzzy searching is conducted based on --adaptor_sequence_errors and --barcode_errrors. It is recomended a minimum sequence length of 12 nucleotides be used with two errors, increasing by one error for every additional 8 base pairs, assuming the sequence difference between all barcodes is greater than three. The pipeline will return an error if the hamming distance of barcodes is less than the number of allowable errors.

INSERT DEMULTIPLEXING SCHEMATIC HERE

After demultiplexing telomeric reads are filtered to ensure the telomere was completely sequenced by removing reads containing more than 30% (--subtelo_threshold) one nucleotide telomeric repeat variants in the first  300 base pairs (--min_subtelo_length). Telomere length determination is divided into two steps. The first is the identification of the Variable Repeat Rich (VRR) Region Start, followed by telomere length calculation.  During development of the pipeline it was clearly seen in both HG002 data and data collected from clinical samples, that while there was rarely a sharp and sudden increase from 0 to 100% of perfect telomeric repeats, there was a region immediately adjacent to the telomere composed of nearly 100% one nucleotide telomere repeat variants (see Deimler et al. for more information). The beginnig of this VRR region was often a very clear and suddent increase from 0% to 100%. Computationally, this region is identified by a sliding window 200bp in size (--sliding_window_size) jumping 10bp (--sliding_window_interval) at a time. Once this window reaches a threshold of being composed of 60% (--upper_threshold) one nucleotide telomere repeat variants the VRR region starts at the first one nucleotide variant repeat sequence within the sliding window, assuming it does not drop below 10% (--lower_threshold). These thresholds were determined by computationally comparing to a manually annotated truth set. The lower threshold is required to ensure that one nucleotide variation islands that exist within the subtelomere (as seen on chromosome arm WHICH CHROMOSOME ARM) do not drastically skew any length estimations. However, this entire VRR region would not be functionally relevant. Therefore telomere length is calculated as the summation of all perfect telomere repeats that exist between the VRR start and the end of the read when they occur at least 3 repeats in a row (--consecutive_repeats). All statistics and plots use this telomere length calculation metric unless --plot_vrr_length is specified. In addition, a easily readable HTML file is generated containing all relevant statistics and plots. See the output section of this documentation for a more detailed description of the output produced by this pipeline.

## Installation

TArPON is a nextflow pipeline and is readily integrable into Epi2Me but can also be installed on the command line. Nextflow must also be installed on the computer which requires a java installation.  Please see https://www.nextflow.io/docs/latest/install.html for more information on installing nextflow and java. For more information on NextFlow please see https://www.nextflow.io/docs/latest/index.html. 

To install Tarpon on the command line simply clone this github repository and ensure Docker is installed on your system. Nextflow will automatically pull the appropriate docker images from dockerhub the first time the program is executed to ensure no dependency issues arise.

    git clone git@github.com:ndeimler99/TArPON.git
    chmod +x TArPON/bin/*
 

## Running TArPON through Epi2Me

Epi2Me can be access and downloaded at https://epi2me.nanoporetech.com/. TArPON can be directly imported into the Epi2Me application by clicking View Workflows > Import Workflows (top right corner) and adding the following url into the dialogue window: https://github.com/ndeimler99/TArPON

## Running TArPON from the Command Line
    
#### Basic Example
#### To activate a boolean such as c_strand_only
#### To include restriction digest

## Input Parameters
| Parameter      | Epi2Me Appearance |Description | Type | Default     |
| :---        |    :----:   | :----:   |         ---: |
| run_name      | Run Name |     | | String| Run | 
| input_file   | Input File |  FastQ File from Nanopore Sequencing for Analysis  |File| ./data/test_fastq.fastq.gz |
| outdir| Output Directory | Location of where you would like the Pipeline to output results | Path | ./output|
| overwrite_outdir | Overwrite Out Directory | If the output directory you have specified already exists, delete the directory |Boolean | False |
| sample_name | Sample Name | Either sample file or sample name must be specified. Use sample name when not demultiplexing and all reads in the input file belong to the same sample. See sample_file | sample |
| sample_file | Sample File | A comma separated values (csv) file containing two columns without headers, "sample_name,barcode_sequence". See test_data/example_sample_file.csv. See sample_name| None |
| trace_dir | | Location for pipeline execution information including CPU usage and time | outdir/pipeline_info |
| c_strand_only | C Strand Only | Boolean value to dictate if the data was collected in a manner at which only C strand telomeric sequences would be expected | False |
| Repeat | Repeat | Telomeric repeat of interest composing Telomeric Sequences | GGTTAG|
| minimum_telo_reads_per_sample | Minimum Telo Reads Per Sample | Minimum number of telomeric reads per sample for the sample to pass filtering criteria and considered valid for drawing conclusions | 3000 |
| repeat_count | Repeat Count | Number of telomeric repeats required in a read for it to be identified as putatively telomeric | 10 |
| reverse_complement_threshold | Reverse Complement Threshold | | 0.80|
| adaptor_sequence | Adaptor Sequence | Sequence that marks the end of telomeric repeats in G strand orientation, even if c_strand_only is activated. See sample file.  If no adaptor sequence is provided the barcodes listed in sample file will be considered the end of the telomeric sequence. If both a sample_file and an adaptor sequence are provided, the adaptor sequence will mark the end of the telomere, while the barcodes within sample_file will be used for demultiplexing | None |
| adaptor_sequence_errors | Adaptor Sequence Errors | Number of Errors to allow within Adaptor Sequence. It is recomended this value is 1/8th of the sequence length. For example, an adaptor of 16 nucleotides should contain 2 errors, while an adaptor of 24 nucleotides should contain 3 | 2 |
| barcode_errors | Barcode Errors | See adaptor sequence errors, but for barcodes provided in sample_file | 3 |
| min_subtelo_length | Min Subtelomere Length | Minimum length of sequence that is less than subtelo_threshold at the beginning of the read to ensure complete sequence of the telomere. Additionally, this removes reads that are less than this length | 300 |
| subtelo_threshold | Subtelomere Threshold | Proportion of read from the start to min_subtelo_length that must not be exceeded in order to consider the read to be completely sequenced through the telomere. For example a subtelo threshold of 0.3 with a min_subtelo_length of 300 indicates the first 300 nucleotides of the read are not telomeric (in G strand conformation). | 0.3 |
|sliding_window_size|Sliding Window Size| To determine telomere start position the read is analyzed using a sliding window approach, this is the size of said sliding window |200|
|sliding_window_interval|Sliding Window Interval| The Interval will not move one bp at a time, but will rather jump. This is the jumping interval of the sliding window|10|
|upper_threshold| Upper Threshold | When the sliding window is composed of this ratio of one nucleotide variants of the telomeric repeat, the VRR region starts assuming it does not drop below lower threshold |0.6| 
|lower_threshold| Lower Threshold|The lower threshold (see upper threshold) at which a the proportion of one nucleotide variants of the telomeric repeat must not drop below in a sliding window |0.1|
|telomeric_repeat_percentage| Telomeric Repeat Percentage | From VRR start to the end of the telomere (identified by the presence of adaptors or barcodes) the telomere must be composed of this percentage of one nucleotide variations of the telomeric repeat |0.7|
|consecutive_repeats| Consecutive Repeats | The entire VRR region is not functional due to the inclusion of one nucleotide variant repeats; however this is highly composed of perfect telomeric repeats. See publication. Telomere length is then calculated by the summation in length of all perfect telomeric repeats that exist in consecutive groups of at least this number or more |3|
|strand_comparison| Strand Comparison | Results in the creation of additional statistics and plots comparing G strand telomeric repeats and C strand telomeric repeats|False|
|indiv_read_plots|Individual Read Plots| Will create one plot for every telomeric read with percentage of telomeric repeats on the Y-axis and location within the read on the X-axis for both perfect repeats and one nucleotide variations. In addition, the start of the VRR region is denoted by a solid red, vertical line |False|
|plot_vrr_length| Plot VRR Length| Output plots and statistics will use VRR length, not telomere length |False|
|plot_telo_length| Plot Telo Length | Output plots and statistics will use telomere length (see consecutive_repeats) instead of VRR lenth |True|
|detailed_stats| Detailed Stats | Results in the output of additional plots and figures for more indepth comparison and characterization of telomeric sequences |False|
|remove_wd| Remove WD| Removes Nextflow Working Directory. Setting this option to True will save space, but not allow for debugging or checking results mid pipeline |False|
|restriction_digest_analysis| Restriction Digest Analysis|A comma separated list of restriction enzyme cut-sites to determine the completeness of digestion of telomeric sequences|None|

## Ouput