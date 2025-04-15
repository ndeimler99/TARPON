#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)
library(stringr)

COLORS=c('#0173b2', '#de8f05', '#029e73', '#d55e00', '#cc78bc', '#ca9161', '#fbafe4', '#949494', '#ece133', '#56b4e9')

args = commandArgs(trailingOnly=TRUE)


stats <- read.table("combined_stats.VRR.txt", header=TRUE)

stats$Sample_ID <- factor(stats$Sample_ID, levels=sort(stats$Sample_ID))
  # number of reads
read_count <- ggplot(data=stats) + 
    geom_bar(mapping=aes(x=Sample_ID, y=Number_of_Reads), stat='identity', fill=COLORS[1:length(stats$Sample_ID)]) +
    theme_minimal() + ylab("Number of Telomeric Reads") +
    theme(axis.title.x = element_blank(),
          axis.text=element_text(size=15),
          axis.title.y = element_text(size=20),       
          axis.text.x = element_text(size=15, angle=45, vjust=1, hjust=1))
  
ggsave("sample_comparison.read_count.pdf", device=pdf, plot=read_count, width=12, height=10)
  
  # telo length box plot differences
telo_length_box <- ggplot(data=stats) +
    geom_boxplot(mapping=aes(x=Sample_ID, ymin=Min_VRR_Telo_Length, lower=Q1, middle=Q2, upper=Q3, ymax=Max_VRR_Telo_Length), stat="identity", fill=COLORS[1:length(stats$Sample_ID)]) +
    theme_minimal() +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_text(size=15, angle=45, hjust=1, vjust=1),
          axis.title.y = element_text(size=20),
          axis.text.y = element_text(size=15)) +
    ylab("VRR Telomere Length")
  
  
ggsave("sample_comparison.vrr_telo_length.boxplot.pdf", device=pdf, plot=telo_length_box, width=12, height=10)
