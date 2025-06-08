#!/usr/bin/env Rscript

# load appropriate libraries
library(ggplot2)
library(dplyr)

# get args args[1] = telo_stats (file), args[2] = plot_telo_length (boolean), args[3] = plot_vrr_length (boolean), args[4] = strand_comparison (boolean)
args = commandArgs(trailingOnly=TRUE)
# open telo stats file into dataframe
telo_stats <- read.table(args[1], header=TRUE)

# of reads in clusters
# table of cluster summary stats
# telomere length by cluster boxplots sorted from shortest to longest with text saying how many telomere it is
# barplot showing nujmber of reads in each cluster sorted from shortest to longest by length