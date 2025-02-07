# TARPON
Telomere Analysis and Research Pipeline Optimized for Nanopore Sequencing Data

## Table of Contents
1. [What is TARPON and Why Should I Use It?](#what_and_why)
2. [The Pipeline](#pipeline)
3. [Installation](#installation)
4. [Running Tarpon Through Epi2Me](#epi2me)
5. [Running Tarpon Through the Command Line](#commandline)
6. [Input Paramaters](#input)
7. [Advanced Input Parameteres](#advanced_input)
8. [Output Files](#output)


## [What is TARPON and Why Should I Use It?](#what_and_why)

TARPON is a pipeline developed in the lab of Dr. Peter Baumann by PhD candidate Nathaniel Deimler. At the time of the pipeline's creation there were multiple publications (cite publications) that have developed protocols for the enrichment of telomeric sequences to succesfully sequence human telomeres by Nanopore sequencing while remaining cost productive. However, all of these papers focused on the biological aspects of telomeres and the mechanical or chemical enrichment of telomeres for sequencing, not on the computational analysis of individual telomere reads - something never before capable. Unfortunately, all of these papers use different methods to calculate and identify telomeric reads and provide little to no justificiation on why those parameters and techniques were chosen. TArPON aims to provide a one stop shop pipeline for clinicians as well as researchers interested in telomere dynamics using a set of best practices that have been justified in Deimler et al. TArPON is capable of being run in a command line environment or through the ONT GUI Epi2Me and is compatible with all pre-existing methods of telomere sequencing and enrichment including the ONT telo-seq methodology.

## [The Pipeline](#pipeline)

<img src="./src/TArPON_Pipeline.png" alt="pipeline_schematic" width=200>


TARPON begins by identifying a set of putative telomeric sequences to reduce computational time and resources at later steps in the pipeline. By default, this is done by isolating reads that contain at least ten perfect telomeric repeats (see parameters: repeat, repeat_count ). These repeats do not have to be consecutive.  The ratio of G telomeric repeats to C telomeric repeats is calculated for all putative telomeric reads. If the ratio falls between 20 and 80% ( see parameter reverse_complement_threshold ) the reads are discarded as they are most likely chimeric sequences.  All reads that are less than 20% G telomeric repeats (these would be C strand reads) are reverse complemented so all telomeric reads are in the same orientation for further analysis. Strandedness information is retaiend throughout the pipeline and strands can be compared using the parameter --strand_comparison.  Telomeric reads are then demultiplexed and the end of the read is identified using either a provided adaptor sequence (--adaptor_sequence) or by the barcodes found in --sample_file. If no sample_file is provided all reads are considered to end in the adaptor_sequence and belong to the same sample: in this case --sample_name should also be provided.  In the case that all telomeric capture probes/sequences are unique and do not contain a shared sequence, --adaptor_sequence should be left blank and the end of the telomeric read will be determined by the barcodes present in sample_file. If both adaptor sequence and sample_file are provided, the end of the telomere will be determined by the presence of the adaptor_sequence and demultiplexed by the barcodes found in sample_file. In all cases of sequence searches, fuzzy searching is conducted based on --adaptor_sequence_errors and --barcode_errrors. It is recomended a minimum sequence length of 12 nucleotides be used with two errors, increasing by one error for every additional 8 base pairs, assuming the sequence difference between all barcodes is greater than three. The pipeline will return an error if the hamming distance of barcodes is less than the number of allowable errors.

<img src="./src/demux_figure.png" alt="demux_schematic" width=800>

After demultiplexing telomeric reads are filtered to ensure the telomere was completely sequenced by removing reads containing more than 30% (--subtelo_threshold) one nucleotide telomeric repeat variants in the first  300 base pairs (--min_subtelo_length). Telomere length determination is divided into two steps. The first is the identification of the Variable Repeat Rich (VRR) Region Start, followed by telomere length calculation.  During development of the pipeline it was clearly seen in both HG002 data and data collected from clinical samples, that while there was rarely a sharp and sudden increase from 0 to 100% of perfect telomeric repeats, there was a region immediately adjacent to the telomere composed of nearly 100% one nucleotide telomere repeat variants (see Deimler et al. for more information). The beginnig of this VRR region was often a very clear and suddent increase from 0% to 100%. Computationally, this region is identified by a sliding window 100bp in size (--sliding_window_size) jumping 15bp (--sliding_window_interval) at a time. Once this window reaches a threshold of being composed of 60% (--upper_threshold) one nucleotide telomere repeat variants the VRR region starts at the first one nucleotide variant repeat sequence within the sliding window, assuming it does not drop below 5% (--lower_threshold) for 15 consecutive windows (--consecutive_threshold). These thresholds were determined by computationally comparing to a manually annotated truth set. The lower threshold is required to ensure that one nucleotide variation islands that exist within the subtelomere (as seen on chromosome arm WHICH CHROMOSOME ARM) do not drastically skew any length estimations. However, the functionality of this VRR region is unknown and an additional paramter (--plot_telo_length) calculates the summation of all perfect telomere repeats that exist between the VRR start and the end of the read when they occur at least 3 repeats in a row (--consecutive_repeats). By default the pipeline using VRR Telomere Length for all statistics and plots unles --plot_telo_length is specified. In addition, a easily readable HTML file is generated containing all relevant statistics and plots. See the output section of this documentation for a more detailed description of the output produced by this pipeline.

## [Installation](#installation)

TARPON is a nextflow pipeline and is readily integrable into Epi2Me but can also be installed on the command line. Nextflow must also be installed on the computer which requires a java installation.  Please see https://www.nextflow.io/docs/latest/install.html for more information on installing nextflow and java. For more information on NextFlow please see https://www.nextflow.io/docs/latest/index.html. 

To install TARPON on the command line simply clone this github repository and ensure Docker is installed on your system. Nextflow will automatically pull the appropriate docker images from dockerhub the first time the program is executed to ensure no dependency issues arise.

    git clone git@github.com:ndeimler99/TARPON.git
    chmod +x TARPON/bin/*
 
 Installation can be tested by running the following commands

    nextflow run main.nf --version
<br>
    
    nextflow run main.nf --help
<br>

    nextflow run main.nf --input_file test_data/test_1.fastq.gz --sample_file test_data/sampleFile.csv --adaptor NNNNNNNNNNNN

## [Running TARPON through Epi2Me](#epi2me)

Epi2Me can be access and downloaded at https://epi2me.nanoporetech.com/. TARPON can be directly imported into the Epi2Me application by clicking View Workflows > Import Workflows (top right corner) and adding the following url into the dialogue window: https://github.com/ndeimler99/TARPON. 

## [Running TARPON from the Command Line](#commandline)

In the simplest form TARPON can be run using the following command replacing "file" with the path of a bam file or directory and NNNNNNNNNN with the appropriate adaptor sequence.

    nextflow run main.nf --input_file "file" --adaptor NNNNNNNNNN

To enable demultiplexing an additional sample file will have to be provided

    nextflow run main.nf --input_file "file.bam" --adaptor NNNNNNNNNN --sample_file "sampleFile.csv"

To override any default paramters the parameter must be specified with --parameter_name parameter_value. For example to override the run_name.

    nextflow run main.nf --input_file "file" --adaptor NNNNNNNNNN --run_name TEST_RUN

To activate a boolean parameter such as strand comparison or c_strand_only no value needs to be provided after the parameter name

    nextflow run main.nf --input_file "file.fastq.gz" --sample_file "sampleFile.csv" --c_strand_only --run_name TEST_C_STRAND_ONLY
  
To include a restriction digest analysis use the paramter --restriction_digest_analysis with a comma separated list of cut sites. For example searching for EcoRV and EcoRI cut sites.

    nextflow run main.nf --input_file "file.fastq.gz" --sample_file "sampleFile.csv" --restriction_digest_analysis GATATC,GAATTC

## [Basic Input Parameters](#input)
| Parameter      | Epi2Me Appearance |Description | Type | Default     |
| :---        |    :----:   | :----:   | :---: |       ---: |
| run_name      | Run Name |  Name of Sequencing Run for Overall Statistics and printed on html report  | String | String| Run | 
| input_file   | Input File |  Bam File, FastQ File, compressed fastq file, or Directory from Nanopore Sequencing for Analysis. If a directory, all bam files will be merged together. If fastq files are provided or a directory containing fastq files is provided, the fastq files will first be merged together and converted to a bam file prior to analysis  | File or Directory| ./test_data/test.bam |
| fast_basecalled | Fast Basecalled | TARPON allows for data basecalled using the fast models to be analysed.  If this parameter is specified telomere reads will be isolated from the sequencing dataset and rebasecalled using the specified pod5_directory using a SUP model. This is designed to reduce turn around time. As opposed to SUP basecalling the entire dataset, only relevant reads are SUP basecalled. | Boolean | False |
| pod5_dir | Pod5 Directory | Location of raw pod5 files from sequencing. Must be specified when using --fast_basecalled | Directory | None |
| sample_name | Sample Name | Either sample file or sample name must be specified. Use sample name when not demultiplexing and all reads in the input file belong to the same sample. See sample_file | String | sample | 
| sample_file | Sample File | A comma separated values (csv) file containing two columns without headers, "sample_name,barcode_sequence". See test_data/example_sample_file.csv. See sample_name| CSV File | None |
| Repeat | Repeat | Telomeric repeat of interest composing Telomeric Sequences | String | GGTTAG|
| outdir| Output Directory | Location of where you would like the Pipeline to output results | Path | ./output|
| trace_dir | | Location for pipeline execution information including CPU usage and time | Directory | outdir/pipeline_info |

## [Advanced Input Parameteres](#advanced_input)

| Parameter      | Epi2Me Appearance |Description | Type | Default     |
| :---        |    :----:   | :----:   | :---: |       ---: |
| overwrite_outdir | Overwrite Out Directory | If the output directory you have specified already exists, delete the directory |Boolean | False |
| c_strand_only | C Strand Only | Boolean value to dictate if the data was collected in a manner at which only C strand telomeric sequences would be expected | Boolean | False |
| minimum_telo_reads_per_sample | Minimum Telo Reads Per Sample | Minimum number of telomeric reads per sample for the sample to pass filtering criteria and considered valid for drawing conclusions | Integer | 500 |
| repeat_count | Repeat Count | Number of telomeric repeats required in a read for it to be identified as putatively telomeric |  Integer | 10 |
| reverse_complement_threshold | Reverse Complement Threshold | Percentage of telomeric repeats must be greater than this value to be classified as C strand. If the percentage is less than 1 minus this value, the read is classified as C strand. Anthing greater than 1-value but less than value is removed from analysis |  Float | 0.80|
| adaptor_sequence | Adaptor Sequence | Sequence that marks the end of telomeric repeats in G strand orientation, even if c_strand_only is activated. See sample file.  If no adaptor sequence is provided the barcodes listed in sample file will be considered the end of the telomeric sequence. If both a sample_file and an adaptor sequence are provided, the adaptor sequence will mark the end of the telomere, while the barcodes within sample_file will be used for demultiplexing | String| None |
| adaptor_sequence_errors | Adaptor Sequence Errors | Number of Errors to allow within Adaptor Sequence. It is recomended this value is 1/8th of the sequence length. For example, an adaptor of 16 nucleotides should contain 2 errors, while an adaptor of 24 nucleotides should contain 3 | Integer | 2 |
| barcode_errors | Barcode Errors | See adaptor sequence errors, but for barcodes provided in sample_file | Integer | 3 |
| min_subtelo_length | Min Subtelomere Length | Minimum length of sequence that is less than subtelo_threshold at the beginning of the read to ensure complete sequence of the telomere. Additionally, this removes reads that are less than this length | Integer | 300 |
| subtelo_threshold | Subtelomere Threshold | Proportion of read from the start to min_subtelo_length that must not be exceeded in order to consider the read to be completely sequenced through the telomere. For example a subtelo threshold of 0.3 with a min_subtelo_length of 300 indicates the first 300 nucleotides of the read are not telomeric (in G strand conformation). | Float | 0.3 |
|sliding_window_size|Sliding Window Size| To determine telomere start position the read is analyzed using a sliding window approach, this is the size of said sliding window | Integer | 200|
|sliding_window_interval|Sliding Window Interval| The Interval will not move one bp at a time, but will rather jump. This is the jumping interval of the sliding window| Integer | 10|
|upper_threshold| Upper Threshold | When the sliding window is composed of this ratio of one nucleotide variants of the telomeric repeat, the VRR region starts assuming it does not drop below lower threshold | Float | 0.6| 
|lower_threshold| Lower Threshold|The lower threshold (see upper threshold) at which a the proportion of one nucleotide variants of the telomeric repeat must not drop below in a sliding window | Float | 0.05|
|consecutive_threshold| Consecutive Threshold | Number of jumping windows in which the ratio of one nucleotide substitutions of the telomeric repeat must be below the lower threshold for the telomere start to be considered not identified | Integer | 15 |
|telomeric_repeat_percentage| Telomeric Repeat Percentage | From VRR start to the end of the telomere (identified by the presence of adaptors or barcodes) the telomere must be composed of this percentage of one nucleotide variations of the telomeric repeat | Float | 0.7|
|consecutive_repeats| Consecutive Repeats | The entire VRR region is not functional due to the inclusion of one nucleotide variant repeats; however this is highly composed of perfect telomeric repeats. See publication. Telomere length is then calculated by the summation in length of all perfect telomeric repeats that exist in consecutive groups of at least this number or more | Integer | 3|
|strand_comparison| Strand Comparison | Results in the creation of additional statistics and plots comparing G strand telomeric repeats and C strand telomeric repeats| Boolean | False|
|indiv_read_plots|Individual Read Plots| Will create one plot for every telomeric read with percentage of telomeric repeats on the Y-axis and location within the read on the X-axis for both perfect repeats and one nucleotide variations. In addition, the start of the VRR region is denoted by a solid red, vertical line | Boolean | False|
|plot_vrr_length| Plot VRR Length| Output plots and statistics will use VRR length, not telomere length | Boolean | False|
|plot_telo_length| Plot Telo Length | Output plots and statistics will use telomere length (see consecutive_repeats) instead of VRR lenth | Boolean | True|
|detailed_stats| Detailed Stats | Results in the output of additional plots and figures for more indepth comparison and characterization of telomeric sequences | Boolean | False|
|remove_wd| Remove WD| Removes Nextflow Working Directory. Setting this option to True will save space, but not allow for debugging or checking results mid pipeline | Boolean | False|
|restriction_digest_analysis| Restriction Digest Analysis|A comma separated list of restriction enzyme cut-sites to determine the completeness of digestion of telomeric sequences| String |None|
|mutant| USE WITH CAUTION : If a telomerase RNA mutant incorporates mutant telomeric repeats into the telomere sequence, the mutant repeat of interest can be specified here to return mutant specific calculations of processivity and occurence. Note that this function is in it's preliminary development and is used for exploratory purposes only. It is highly recomended your own analysis is done on the isolated telomeric sequences for higher levels of clarity - I am more than happy to assist and program the analysis specific to your mutant of interest.  Specifying this argument also impacts all other functions of TARPON as the mutant repeat will be included as a telomeric repeat.| String | False |

## [Ouput](#output)

Various output files are produced when the pipeline is run with different parameters. In all cases, the output directory is created containing the following subdirectories.
1. pipeline_info
2. A directory for global run statistics
3. A directory per sample labeled with either --sample_name (simplex) or the sample names listed in --sample_file (multiplex)

|Directory|Contents|
|:--|:--:|
| specified output directory| Contains an HTML report with all relevant statistics and plots for easy viewing to open in your web browser, as well as a directory containing statistics and figures from the entire run, a directory per sample if multiplexed, and a directory containing nextflow relevant workflow files|
| pipeline_info | execution_trace.txt - List of processes performed during the execution of the nextflow pipeline and their system requirement<br/>execution_timeline.html - Timeline of all processes showing memory usage and duration throughout pipeline execution<br/>execution_report.html - Detailed analysis on all pipeline executation statistics and commands<br/>pipeline_dag.html - Directed Acyclic Graph of Pipeline Execution showing how all processes interact |
| RUN_STATS | Global statistics on all reads filtered or retained throughout the analysis independent of sample.  Will contain basic .txt files listing relevant statistics and a FIGURES directory showing higher resolution plots of those found in the html report |
| One Directory per Sample | Statistics on each individual sample including filtered/retained reads, as well as indiivudal read telomere length and strand information. Contains a subdirectory FIGURES containing higher resolution plots of those found in the html report |


## Additional Output Files

Specifying one of the booleans below will result in more detailed or additional results/figures being returned with little increase in computational time or cost.  The parameters are stackable and will often produce additional input when stacked. For example combining --detailed_stats with --strand_comparison will produce more information than either parameter alone.  

| Parameter | Location | Description |
| :--------:|:--:|:--:|
|--strand_comparison| <nobr>sample directory / C_G_Comparison / |An additional directory within each sample containing three pdf files with direct comparisons of telo length and read length between the C and G strands|
|<nobr>--restriction_digest_analysis | <nobr> sample directory / restriction_digest.stats.txt |Will create one additional file per sample containing the summary statistics of all reads containing each specified cut site produced by seqkit stats|
|--plot_telo_length|sample directory / FIGURES |Will produce additional files for all predescribed telomere statistics, but using telomere length instead of the vrr telo length. This compounds with --strand_comparison to produce strand specific telo length plots |
|--detailed_stats| sample directory / DETAILED_STATS / |Will produce multiple additional figures comparing telomere length to VRR length, looking at the percentage of perfect repeats within the sequencing data, looking at the quality of telomeric sequences, etc. All such files can be found in the sample directory under FIGURES/DETAILED_STATS. This compounds with both strand_comparison and plot_vrr_length to produce additional files in the C_G_Comparison directory |
