#!/usr/bin/env Rscript
library(ggplot2)
library(dplyr)
library(viridis)
library(ggpointdensity)

args = commandArgs(trailingOnly=TRUE)
telo_stats <- read.table(args[1], header=TRUE)
telo_lengths_for_binning <- c(10500, 9500, 8500, 7500, 6500, 5500, 4500, 3500, 2500, 1500, 500)
read_length_hist <- ggplot(data=telo_stats) +
  geom_histogram(mapping=aes(read_len), binwidth=50) +
  geom_vline(xintercept=mean(telo_stats$read_len), color='red') +
  xlab("Read Length") + ylab("Read Frequency") +
  theme_minimal() +
  theme(axis.text = element_text(size=15),
        axis.title = element_text(size=20))

ggsave("read_length_hist.pdf", device="pdf", plot=read_length_hist, width=10, height=10)

if (args[2]){
  # plot telo length
  telo_length_hist <- ggplot(data=telo_stats) +
    geom_histogram(mapping=aes(telo_length), binwidth=50) +
    geom_vline(xintercept=mean(telo_stats$telo_length), color='red') +
    xlab("Telomere Length (Base Pairs)") + ylab("Read Frequency") +
    theme_minimal() +
    theme(axis.text = element_text(size=15),
          axis.title = element_text(size=20))
  
  ggsave("telo_length_hist.pdf", device="pdf", plot=telo_length_hist, width=10, height=10)
  
  telo_length_scatter <- ggplot(data=telo_stats) +
    geom_pointdensity(mapping=aes(x=telo_length, y=read_len), adjust=4) +
    scale_color_viridis() +
    theme_minimal() +
    theme(axis.text = element_text(size=15),
          axis.title = element_text(size=20)) +
    xlab("Telomere Length (Base Pairs)") +
    ylab("Read Length (Base Pairs)")
  
  ggsave("telo_length_by_read_length.pdf", device="pdf", plot=telo_length_scatter, width=10, height=10)
  
  telo_stats$bin_telo_length <- unlist(lapply(telo_stats$telo_length, function(x) telo_lengths_for_binning[which.min(abs(telo_lengths_for_binning-x))]))
  telo_stats$bin_telo_length <- factor(telo_stats$bin_telo_length, levels=c(10500, 9500, 8500, 7500, 6500, 5500, 4500, 3500, 2500, 1500, 500))
  
  telo_bar_hist <- ggplot(data=telo_stats) +
    geom_bar(mapping=aes(x=1, fill=bin_telo_length), position="fill") +
    theme_minimal() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank()) +
    ylab("Proportion of Telomere Reads") +
    guides(fill=guide_legend(title="Telomere Length (BP)")) +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(breaks=c(10500, 9500, 8500, 7500, 6500, 5500, 4500, 3500, 2500, 1500, 500),
                      labels=c("10000+", "9000-9999", "8000-8999", "7000-7999", "6000-6999", "5000-5999", "4000-4999", "3000-3999", "2000-2999", "1000-1999", "0-999"),
                      values=c("#F8766D", "#E68613", "#ABA300", "#0CB702", "#00BE67", "#00BFC4", "#00A9FF", "#8494FF", "#C77CFF", "#FF61CC", "#FF68A1"))
  
  ggsave("telo_length_bar_plot.pdf", plot = telo_bar_hist, device="pdf", width=6, height=10)
}

if (args[3]){
  # plot VRR length
  vrr_length_hist <- ggplot(data=telo_stats) +
    geom_histogram(mapping=aes(vrr_telo_length), binwidth=50) +
    geom_vline(xintercept=mean(telo_stats$vrr_telo_length), color='red') +
    xlab("VRR Telomere Length (Base Pairs)") + ylab("Read Frequency") +
    theme_minimal() +
    theme(axis.text = element_text(size=15),
          axis.title = element_text(size=20))
  
  ggsave("vrr_length_hist.pdf", device="pdf", plot=vrr_length_hist, width=10, height=10)
  
  vrr_length_scatter <- ggplot(data=telo_stats) +
    geom_pointdensity(mapping=aes(x=vrr_telo_length, y=read_len), adjust=4) +
    scale_color_viridis() +
    theme_minimal() +
    theme(axis.text = element_text(size=15),
          axis.title = element_text(size=20)) +
    xlab("VRR Telomere Length (Base Pairs)") +
    ylab("Read Length (Base Pairs)")
  
  ggsave("vrr_length_by_read_length.pdf", device="pdf", plot=vrr_length_scatter, width=10, height=10)
  
  telo_stats$bin_vrr_length <- unlist(lapply(telo_stats$vrr_telo_length, function(x) telo_lengths_for_binning[which.min(abs(telo_lengths_for_binning-x))]))
  telo_stats$bin_vrr_length <- factor(telo_stats$bin_vrr_length, levels=c(10500, 9500, 8500, 7500, 6500, 5500, 4500, 3500, 2500, 1500, 500))
  
  vrr_bar_hist <- ggplot(data=telo_stats) +
    geom_bar(mapping=aes(x=1, fill=bin_vrr_length), position="fill") +
    theme_minimal() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank()) +
    ylab("Proportion of Telomere Reads") +
    guides(fill=guide_legend(title="VRR Telomere\nLength (BP)")) +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(breaks=c(10500, 9500, 8500, 7500, 6500, 5500, 4500, 3500, 2500, 1500, 500),
                      labels=c("10000+", "9000-9999", "8000-8999", "7000-7999", "6000-6999", "5000-5999", "4000-4999", "3000-3999", "2000-2999", "1000-1999", "0-999"),
                      values=c("#F8766D", "#E68613", "#ABA300", "#0CB702", "#00BE67", "#00BFC4", "#00A9FF", "#8494FF", "#C77CFF", "#FF61CC", "#FF68A1"))

  ggsave("vrr_length_bar_plot.pdf", plot = vrr_bar_hist, device="pdf", width=6, height=10)
}

if (args[4]) {
  # plot strand comparison
  read_length_by_strand <- ggplot(data=telo_stats) +
    geom_histogram(mapping=aes(read_len), binwidth=50) +
    facet_grid(strand ~ .) +
    geom_vline(data=filter(telo_stats, strand=="C"), aes(xintercept=mean(telo_stats[telo_stats$strand == "C",]$read_len)), color='red') +
    geom_vline(data=filter(telo_stats, strand=="G"), aes(xintercept=mean(telo_stats[telo_stats$strand == "G",]$read_len)), color='red') +
    xlab("Read Length") + ylab("Read Frequency") +
    theme_minimal() +
    theme(axis.text = element_text(size=15),
          axis.title = element_text(size=20),
          panel.border = element_rect(color="black", fill=NA))
  
  ggsave("C_G_COMPARISON/C_vs_G.Read_Length.pdf", device="pdf", width=10, height=12, plot=read_length_by_strand, create.dir = TRUE)

  if (args[2]){
    telo_length_by_strand <- ggplot(data=telo_stats) +
      geom_histogram(mapping=aes(telo_length), binwidth=50) +
      facet_grid(strand ~ .) +
      geom_vline(data=filter(telo_stats, strand=="C"), aes(xintercept=mean(telo_stats[telo_stats$strand == "C",]$telo_length)), color='red') +
      geom_vline(data=filter(telo_stats, strand=="G"), aes(xintercept=mean(telo_stats[telo_stats$strand == "G",]$telo_length)), color='red') +
      xlab("Telomere Length") + ylab("Telomere Frequency") +
      theme_minimal() +
      theme(axis.text = element_text(size=15),
            axis.title = element_text(size=20),
            panel.border = element_rect(color="black", fill=NA))
    
    ggsave("C_G_COMPARISON/C_vs_G.Telo_Length.pdf", device="pdf", width=10, height=12, plot=telo_length_by_strand)
  
    
    telo_stats$bin_telo_length <- unlist(lapply(telo_stats$telo_length, function(x) telo_lengths_for_binning[which.min(abs(telo_lengths_for_binning-x))]))
    telo_stats$bin_telo_length <- factor(telo_stats$bin_telo_length, levels=c(10500, 9500, 8500, 7500, 6500, 5500, 4500, 3500, 2500, 1500, 500))
    
    telo_bar_hist <- ggplot(data=telo_stats) +
      geom_bar(mapping=aes(x=strand, fill=bin_telo_length), position="fill") +
      theme_minimal() +
      theme(axis.title.x = element_text(size=20),
            axis.text.x = element_text(size=15)) +
      ylab("Proportion of Telomere Reads") + xlab("Strand") +
      guides(fill=guide_legend(title="Telomere Length (BP)")) +
      scale_y_continuous(labels = scales::percent) +
      scale_fill_manual(breaks=c(10500, 9500, 8500, 7500, 6500, 5500, 4500, 3500, 2500, 1500, 500),
                        labels=c("10000+", "9000-9999", "8000-8999", "7000-7999", "6000-6999", "5000-5999", "4000-4999", "3000-3999", "2000-2999", "1000-1999", "0-999"),
                        values=c("#F8766D", "#E68613", "#ABA300", "#0CB702", "#00BE67", "#00BFC4", "#00A9FF", "#8494FF", "#C77CFF", "#FF61CC", "#FF68A1"))

    ggsave("C_G_COMPARISON/C_vs_G.telo_length_bar_plot.pdf", plot = telo_bar_hist, device="pdf", width=9, height=10, create.dir = TRUE)
  
  }

  if (args[3]){
    
    vrr_length_by_strand <- ggplot(data=telo_stats) +
      geom_histogram(mapping=aes(vrr_telo_length), binwidth=50) +
      facet_grid(strand ~ .) +
      geom_vline(data=filter(telo_stats, strand=="C"), aes(xintercept=mean(telo_stats[telo_stats$strand == "C",]$vrr_telo_length)), color='red') +
      geom_vline(data=filter(telo_stats, strand=="G"), aes(xintercept=mean(telo_stats[telo_stats$strand == "G",]$vrr_telo_length)), color='red') +
      xlab("VRR Telomere Length") + ylab("Telomere Frequency") +
      theme_minimal() +
      theme(axis.text = element_text(size=15),
            axis.title = element_text(size=20),
            panel.border = element_rect(color="black", fill=NA))
    
    ggsave("C_G_COMPARISON/C_vs_G.VRR_Telo_Length.pdf", device="pdf", width=10, height=12, plot=vrr_length_by_strand, create.dir = TRUE)
  
    telo_stats$bin_vrr_length <- unlist(lapply(telo_stats$vrr_telo_length, function(x) telo_lengths_for_binning[which.min(abs(telo_lengths_for_binning-x))]))
    telo_stats$bin_vrr_length <- factor(telo_stats$bin_vrr_length, levels=c(10500, 9500, 8500, 7500, 6500, 5500, 4500, 3500, 2500, 1500, 500))
    
    vrr_bar_hist <- ggplot(data=telo_stats) +
      geom_bar(mapping=aes(x=strand, fill=bin_vrr_length), position="fill") +
      theme_minimal() +
      theme(axis.title.x = element_text(size=20),
            axis.text.x = element_text(size=15)) +
      ylab("Proportion of Telomere Reads") + xlab("Strand") +
      guides(fill=guide_legend(title="VRR Telomere\nLength (BP)")) +
      scale_y_continuous(labels = scales::percent) +
      scale_fill_manual(breaks=c(10500, 9500, 8500, 7500, 6500, 5500, 4500, 3500, 2500, 1500, 500),
                        labels=c("10000+", "9000-9999", "8000-8999", "7000-7999", "6000-6999", "5000-5999", "4000-4999", "3000-3999", "2000-2999", "1000-1999", "0-999"),
                        values=c("#F8766D", "#E68613", "#ABA300", "#0CB702", "#00BE67", "#00BFC4", "#00A9FF", "#8494FF", "#C77CFF", "#FF61CC", "#FF68A1"))
    
    ggsave("C_G_COMPARISON/C_vs_G.vrr_telo_length_bar_plot.pdf", plot = vrr_bar_hist, device="pdf", width=9, height=10, create.dir = TRUE)
    
  }
}
  

