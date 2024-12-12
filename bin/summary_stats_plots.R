#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)
library(stringr)

# script funtions to generate summary stats of entire sequencing run/experiment
args = commandArgs(trailingOnly=TRUE)
read_stats <- read.table(args[1], header=TRUE, sep="\t")
quality_stats <- read.table(args[4], header=TRUE, sep="\t")

read_stats$file <- str_split_i(read_stats$file, "/", -1)
read_stats$file <- str_split_i(read_stats$file, ".fastq*", 1)
quality_stats$Sample <- str_split_i(quality_stats$Sample, "/", -1)
quality_stats$Sample <- str_split_i(quality_stats$Sample, ".fastq*", 1)

if(args[3] == "telomeric"){
  if (args[2]){
    order <- c('input', 'putative_reads', 'putative_reads.filtered', 'putative_reads.filtered.c_strand', 'putative_reads.filtered.g_strand', 'subtelo_pass', 'subtelo_pass.c_strand', 'subtelo_pass.g_strand', 'adaptor', 'adaptor.c_strand', 'adaptor.g_strand', 'telo_retained', 'telo_retained.c_strand', 'telo_retained.g_strand')
  }
  else{
    order <- c('input', 'putative_reads', 'putative_reads.filtered', 'subtelo_pass', 'adaptor', 'telo_retained')
  }
} 

if(args[3] == "filtered"){
  if (args[2]){
    order <- c('input', 'non_telomeric', '20_80_removed_reads', 'subtelo_fail', 'subtelo_fail.c_strand', 'subtelo_fail.g_strand', 
               'adaptor_filtered', 'adaptor_filtered.c_strand', 'adaptor_filtered.g_strand', 'no_telo_start', 'no_telo_start.c_strand', 'no_telo_start.g_strand',
               "below_telo_threshold", "below_telo_threshold.c_strand", "below_telo_threshold.g_strand")
  }
  else{
    order <- c("input", "non_telomeric", "20_80_removed_reads", "subtelo_fail", "adaptor_filtered", "no_telo_start", "below_telo_threshold")
  }
}


read_stats$file <- factor(read_stats$file, levels=order)
quality_stats$Sample <- factor(quality_stats$Sample, levels=order)

# boxplot length
length_box <- ggplot(data=read_stats) +
  geom_boxplot(mapping=aes(ymin=min_len,lower=Q1, middle=Q2, upper=Q3, ymax=max_len, x=file), stat='identity') +
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        axis.text = element_text(size=15),
        axis.title = element_text(size=20)) +
  coord_flip() + ylab("Read Length")

length_box_zoomed <- ggplot(data=read_stats) +
  geom_boxplot(mapping=aes(ymin=min_len,lower=Q1, middle=Q2, upper=Q3, ymax=max_len, x=file), stat='identity') +
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        axis.text = element_text(size=15),
        axis.title = element_text(size=20)) +
  ylim(0,50000) + coord_flip() + ylab("Read Length")

ggsave(paste(args[3], ".length_boxplot.pdf", sep=''), device="pdf", plot=length_box, width=12, height=10)
ggsave(paste(args[3], ".length_boxplot_zoomed.pdf", sep=''), device="pdf", plot=length_box_zoomed, width=12, height=10)

# quality score plot
qscore <- ggplot(data=quality_stats) +
  geom_boxplot(mapping=aes(x=Sample, ymin=Min,lower=Q1, middle=Q2, upper=Q3, ymax=Max),stat='identity') +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size=15),
        axis.title=element_text(size=20)) +
  coord_flip() + ylab("Quality Score")

ggsave(paste(args[3], ".quality_score.pdf", sep=''), device="pdf", plot=qscore, width=12, height=10)

# number of reads
read_counts <- ggplot(data=read_stats) +
  geom_bar(mapping=aes(x=file, y=num_seqs), stat='identity') +
  theme_minimal() +
  theme(axis.title = element_text(size=20),
        axis.title.x = element_blank(),
        axis.text = element_text(size=15)) +
  ylab("Number of Sequences")

read_counts_no_input <- ggplot(data=read_stats[!read_stats$file %in% c("input", 'non_telomeric'),]) +
  geom_bar(mapping=aes(x=file, y=num_seqs), stat='identity') +
  theme_minimal() +
  theme(axis.title = element_text(size=20),
        axis.title.x = element_blank(),
        axis.text = element_text(size=15)) +
  ylab("Number of Sequences")

ggsave(paste(args[3], ".read_counts.pdf", sep=''), device="pdf", plot=read_counts, width=12, height=10)
ggsave(paste(args[3], ".read_counts.no_input.pdf", sep=''), device="pdf", plot=read_counts_no_input, width=12, height=10)
