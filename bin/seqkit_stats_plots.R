#!/usr/bin/env Rscript
.libPaths("/home/ndeimler/R/x86_64-pc-linux-gnu-library/4.3")
library(ggplot2)
library(dplyr)
library(stringr)

args = commandArgs(trailingOnly=TRUE)

read_stats <- read.table(args[1], header=TRUE)

read_stats$file <- str_split_i(read_stats$file, "/", -1)
read_stats$file <- str_split_i(read_stats$file, ".fastq", 1)

if(args[3] == "telomeric"){
  if (args[2]){
    order <- NULL 
  }
  else{
    order <- c('input', 'putative_reads', 'putative_reads.c_g_filtered', 'subtelo', 'adaptor', 'telomeric')
  }
} 

print(order)

if(args[3] == "filtered"){
  if (args[2]){
    order <- NULL 
  }
  else{
    order <- c("input", "non_telomeric", "20_80_removed_reads", "subtelo_filtered", "adaptor_filtered", "below_telo_%_threshold")
  }
}

print(order)

read_stats <- read_stats[match(order, read_stats$file),]
read_stats$file <- factor(read_stats$file, levels=order)


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
qscore <- ggplot(data=read_stats) +
  geom_line(mapping=aes(x=file, y=AvgQual), group=1) +
  ylab("Average Quality Score") +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size=15),
        axis.title=element_text(size=20))

# ggsave(paste(args[3], ".quality_score.pdf", sep=''), device="pdf", plot=qscore, width=12, height=10)

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

# ggsave(paste(args[3], ".read_counts.pdf", sep=''), device="pdf", plot=read_counts, width=12, height=10)
# ggsave(paste(args[3], ".read_counts.no_input.pdf", sep=''), device="pdf", plot=read_counts_no_input, width=12, height=10)
