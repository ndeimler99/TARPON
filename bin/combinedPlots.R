#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)
library(viridis)
library(ggpointdensity)

# scripts serves the purposes to combine all multiplexed data into one set of easy to understand figures for easy sample comparison

args = commandArgs(trailingOnly=TRUE)
telo_stats <- read.table(args[1], header=TRUE, sep=",")

telo_lengths_for_binning <- c(10500, 9500, 8500, 7500, 6500, 5500, 4500, 3500, 2500, 1500, 500)

telo_stats$bin_telo_length <- unlist(lapply(telo_stats$telo_length, function(x) telo_lengths_for_binning[which.min(abs(telo_lengths_for_binning-x))]))
telo_stats$bin_telo_length <- factor(telo_stats$bin_telo_length, levels=c(10500, 9500, 8500, 7500, 6500, 5500, 4500, 3500, 2500, 1500, 500))

telo_bar_hist <- ggplot(data=telo_stats) +
  geom_bar(mapping=aes(x=sample, fill=bin_telo_length), position="fill") +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size=15),
        axis.title = element_text(size=20),
        axis.text.x = element_text(angle=45)) +
  ylab("Proportion of Telomere Reads") +
  guides(fill=guide_legend(title="Telomere Length (BP)")) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(breaks=c(10500, 9500, 8500, 7500, 6500, 5500, 4500, 3500, 2500, 1500, 500),
                    labels=c("10000+", "9000-9999", "8000-8999", "7000-7999", "6000-6999", "5000-5999", "4000-4999", "3000-3999", "2000-2999", "1000-1999", "0-999"),
                    values=c("#F8766D", "#E68613", "#ABA300", "#0CB702", "#00BE67", "#00BFC4", "#00A9FF", "#8494FF", "#C77CFF", "#FF61CC", "#FF68A1"))

ggsave("sampleComparison.telo_length_bar_plot.pdf", plot = telo_bar_hist, device="pdf", width=length(unique(telo_stats$sample))*4, height=10)

telo_box <- ggplot(data=telo_stats) +
  geom_boxplot(mapping=aes(x=sample,y=telo_length)) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle=45),
        axis.text = element_text(size=15),
        axis.title = element_text(size=20)) +
  ylab("Telomere Length")

ggsave("sampleComparison.telo_length_box_plot.pdf", plot = telo_box, device="pdf", width=length(unique(telo_stats$sample))*4, height=10)


if (args[2]) {
  # also plot vrr length
  
  telo_stats$bin_telo_length <- unlist(lapply(telo_stats$vrr_telo_length, function(x) telo_lengths_for_binning[which.min(abs(telo_lengths_for_binning-x))]))
  telo_stats$bin_telo_length <- factor(telo_stats$bin_telo_length, levels=c(10500, 9500, 8500, 7500, 6500, 5500, 4500, 3500, 2500, 1500, 500))
  
  
  telo_bar_hist <- ggplot(data=telo_stats) +
    geom_bar(mapping=aes(x=sample, fill=bin_telo_length), position="fill") +
    theme_minimal() +
    theme(axis.title.x = element_blank(),
          axis.text = element_text(size=15),
          axis.title = element_text(size=20),
          axis.text.x = element_text(angle=45)) +
    ylab("Proportion of Telomere Reads") +
    guides(fill=guide_legend(title="VRR Telomere Length (BP)")) +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(breaks=c(10500, 9500, 8500, 7500, 6500, 5500, 4500, 3500, 2500, 1500, 500),
                      labels=c("10000+", "9000-9999", "8000-8999", "7000-7999", "6000-6999", "5000-5999", "4000-4999", "3000-3999", "2000-2999", "1000-1999", "0-999"),
                      values=c("#F8766D", "#E68613", "#ABA300", "#0CB702", "#00BE67", "#00BFC4", "#00A9FF", "#8494FF", "#C77CFF", "#FF61CC", "#FF68A1"))
  
  ggsave("sampleComparison.VRR_telo_length_bar_plot.pdf", plot = telo_bar_hist, device="pdf", width=length(unique(telo_stats$sample))*4, height=10)
  
  telo_box <- ggplot(data=telo_stats) +
    geom_boxplot(mapping=aes(x=sample,y=vrr_telo_length)) +
    theme_minimal() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle=45),
          axis.text = element_text(size=15),
          axis.title = element_text(size=20)) +
    ylab("VRR Telomere Length")

  ggsave("sampleComparison.VRR_telo_length_box_plot.pdf", plot = telo_box, device="pdf", width=length(unique(telo_stats$sample))*4, height=10)
}