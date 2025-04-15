#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)
library(stringr)

COLORS=c('#0173b2', '#de8f05', '#029e73', '#d55e00', '#cc78bc', '#ca9161', '#fbafe4', '#949494', '#ece133', '#56b4e9')

# script funtions to generate summary stats of entire sequencing run/experiment
args = commandArgs(trailingOnly=TRUE)

read_stats <- read.table(args[1], header=TRUE, sep="\t")

read_stats$file <- str_split_i(read_stats$file, "/", -1)
read_stats$file <- str_split_i(read_stats$file, ".bam", 1)

if (args[2]){
  g_strand <- read_stats[grepl(".g_strand", read_stats$file),]
  c_strand <- read_stats[grepl(".c_strand", read_stats$file),]
  read_stats <- read_stats[!grepl("strand", read_stats$file),]
  
  g_strand$file <- str_split_i(g_strand$file, "\\.", 1)
  c_strand$file <- str_split_i(c_strand$file, "\\.", 1)
}

if(args[3] == "telomeric"){
  if (args[2]){
    strand_order <- c('putative_reads', 'adaptor', 'subtelo_pass', "telo_retained")
  }
  
  order <- c('input', 'putative_reads', 'putative_reads.filtered', 'adaptor', 'subtelo_pass', 'telo_retained')
} 

if(args[3] == "removed"){
  if (args[2]){
    strand_order <- c('adaptor_filtered', "subtelo_fail", "below_telo_threshold", "no_telo_start")
  }
  
  order <- c("input", "non_telomeric", "20_80_removed_reads", "adaptor_filtered", "subtelo_fail", "no_telo_start", "below_telo_threshold")
}


read_stats$file <- factor(read_stats$file, levels=order)

if (args[2]){
  g_strand$file <- factor(g_strand$file, levels=strand_order)
  c_strand$file <- factor(c_strand$file, levels=strand_order)
}

# number of reads
read_counts <- ggplot(data=read_stats) +
  geom_bar(mapping=aes(x=file, y=num_seqs), fill=COLORS[1:length(read_stats$file)], stat='identity') +
  theme_minimal() +
  theme(axis.title = element_text(size=20),
        axis.title.x = element_blank(),
        axis.text = element_text(size=15),
        axis.text.x = element_text(size=15, angle=45, vjust=1, hjust=1)) +
  ylab("Number of Sequences")

read_counts_no_input <- ggplot(data=read_stats[!read_stats$file %in% c("input", 'non_telomeric'),]) +
  geom_bar(mapping=aes(x=file, y=num_seqs), stat='identity', fill=COLORS[2:(1+length(read_stats[!read_stats$file %in% c("input", 'non_telomeric'),]$file))]) +
  theme_minimal() +
  theme(axis.title = element_text(size=20),
        axis.title.x = element_blank(),
        axis.text = element_text(size=15),
        axis.text.x = element_text(size=15, angle=45, hjust=1, vjust=1)) +
  ylab("Number of Sequences")

ggsave(paste(args[3], ".read_counts.pdf", sep=''), device="pdf", plot=read_counts, width=12, height=10)
ggsave(paste(args[3], ".read_counts.no_input.pdf", sep=''), device="pdf", plot=read_counts_no_input, width=12, height=10)

if (args[2]){
  g_counts <- ggplot(data=g_strand) +
    geom_bar(mapping=aes(x=file, y=num_seqs), stat='identity', fill=COLORS[1:length(g_strand$file)]) +
    theme_minimal() +
    theme(axis.title=element_text(size=20),
          axis.title.x = element_blank(),
          axis.text=element_text(size=15),
          axis.text.x=element_text(size=15, angle=45, hjust=1, vjust=1)) +
    ylab("Number of Sequences")
  
  ggsave(paste("STRAND_COMPARISON/", args[3], ".g_strand.read_counts.pdf", sep=""), device="pdf", width=12, height=10, create.dir = TRUE)

  c_counts <- ggplot(data=c_strand) +
    geom_bar(mapping=aes(x=file, y=num_seqs), stat='identity', fill=COLORS[1:length(c_strand$file)]) +
    theme_minimal() +
    theme(axis.title=element_text(size=20),
          axis.title.x = element_blank(),
          axis.text=element_text(size=15),
          axis.text.x=element_text(size=15, angle=45, hjust=1, vjust=1)) +
    ylab("Number of Sequences")
  
  ggsave(paste("STRAND_COMPARISON/", args[3], ".c_strand.read_counts.pdf", sep=""), device="pdf", width=12, height=10, create.dir = TRUE)
}


# read lengths

length_box <- ggplot(data=read_stats) +
  geom_boxplot(mapping=aes(ymin=min_length,lower=Q1, middle=Q2, upper=Q3, ymax=max_length, x=file), stat='identity', fill=COLORS[1:length(read_stats$file)]) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size=15),
        axis.title = element_text(size=20),
        axis.text.x = element_text(size=15, angle=45, hjust=1, vjust=1)) +
  ylab("Read Length")

length_box_zoomed <- length_box + ylim(0,50000)

ggsave(paste(args[3], ".length_boxplot.pdf", sep=''), device="pdf", plot=length_box, width=12, height=10)
ggsave(paste(args[3], ".length_boxplot_zoomed.pdf", sep=''), device="pdf", plot=length_box_zoomed, width=12, height=10)

if (args[2]){
  g_length_box <- ggplot(data=g_strand) +
    geom_boxplot(mapping=aes(ymin=min_length,lower=Q1, middle=Q2, upper=Q3, ymax=max_length, x=file), stat='identity', fill=COLORS[1:length(g_strand$file)]) +
    theme_minimal() +
    theme(axis.title.x = element_blank(),
          axis.text = element_text(size=15),
          axis.title = element_text(size=20),
          axis.text.x = element_text(size=15, angle=45, vjust=1, hjust=1)) +
    ylab("Read Length") + ylim(0,50000)
  
  c_length_box <- ggplot(data=c_strand) +
    geom_boxplot(mapping=aes(ymin=min_length,lower=Q1, middle=Q2, upper=Q3, ymax=max_length, x=file), stat='identity', fill=COLORS[1:length(c_strand$file)]) +
    theme_minimal() +
    theme(axis.title.x = element_blank(),
          axis.text = element_text(size=15),
          axis.title = element_text(size=20),
          axis.text.x = element_text(size=15, angle=45, hjust=1, vjust=1)) +
    ylab("Read Length") + ylim(0,50000)

  ggsave(paste("STRAND_COMPARISON/", args[3], ".g_strand_length_boxplot.pdf", sep=''), device="pdf", plot=g_length_box, width=12, height=10, create.dir = TRUE)
  ggsave(paste("STRAND_COMPARISON/", args[3], ".c_strand_length_boxplot.pdf", sep=''), device="pdf", plot=g_length_box, width=12, height=10, create.dir = TRUE)
}


# quality score plot
qscore <- ggplot(data=read_stats) +
  geom_boxplot(mapping=aes(x=file, ymin=min_quality,lower=Q1_qual, middle=Q2_qual, upper=Q3_qual, ymax=max_quality),stat='identity', fill=COLORS[1:length(read_stats$file)]) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size=15),
        axis.title=element_text(size=20),
        axis.text.x = element_text(size=15, angle=45, hjust=1, vjust=1)) +
   ylab("Quality Score")

ggsave(paste(args[3], ".quality_score.pdf", sep=''), device="pdf", plot=qscore, width=12, height=10)

if (args[2]){
  g_qscore <- ggplot(data=g_strand) +
    geom_boxplot(mapping=aes(x=file, ymin=min_quality,lower=Q1_qual, middle=Q2_qual, upper=Q3_qual, ymax=max_quality),stat='identity', fill=COLORS[1:length(g_strand$file)]) +
    theme_minimal() +
    theme(axis.title.x = element_blank(),
          axis.text = element_text(size=15),
          axis.title=element_text(size=20),
          axis.text.x = element_text(size=15, angle=45, hjust=1, vjust=1)) +
    ylab("Quality Score")
  
  c_qscore <- ggplot(data=c_strand) +
    geom_boxplot(mapping=aes(x=file, ymin=min_quality,lower=Q1_qual, middle=Q2_qual, upper=Q3_qual, ymax=max_quality),stat='identity', fill=COLORS[1:length(c_strand$file)]) +
    theme_minimal() +
    theme(axis.title.x = element_blank(),
          axis.text = element_text(size=15),
          axis.title=element_text(size=20),
          axis.text.x = element_text(size=15, angle=45, hjust=1, vjust=1)) +
    ylab("Quality Score")
  
  ggsave(paste("STRAND_COMPARISON/", args[3], ".g_strand_quality_score.pdf", sep=''), device="pdf", plot=g_qscore, width=12, height=10, create.dir = TRUE)
  ggsave(paste("STRAND_COMPARISON/", args[3], ".c_strand_quality_score.pdf", sep=''), device="pdf", plot=c_qscore, width=12, height=10, create.dir = TRUE)
}


