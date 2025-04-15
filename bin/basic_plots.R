#!/usr/bin/env Rscript

# load appropriate libraries
library(ggplot2)
library(dplyr)
library(viridis)
library(ggpointdensity)

# get args args[1] = telo_stats (file), args[2] = plot_telo_length (boolean), args[3] = plot_vrr_length (boolean), args[4] = strand_comparison (boolean)
args = commandArgs(trailingOnly=TRUE)
# open telo stats file into dataframe
telo_stats <- read.table(args[1], header=TRUE)
telo_lengths_for_binning <- c(10500, 9500, 8500, 7500, 6500, 5500, 4500, 3500, 2500, 1500, 500)

# plot read length histogram
read_length_hist <- ggplot(data=telo_stats) +
  geom_histogram(mapping=aes(read_len), binwidth=50) +
  geom_vline(xintercept=mean(telo_stats$read_len), color='red', linetype="dashed") +
  xlab("Read Length") + ylab("Read Frequency") +
  theme_minimal() +
  theme(axis.text = element_text(size=15),
        axis.title = element_text(size=20))

ggsave("read_length_hist.pdf", device="pdf", plot=read_length_hist, width=10, height=10)


  # plot VRR length histogram
vrr_length_hist <- ggplot(data=telo_stats) +
    geom_histogram(mapping=aes(vrr_telo_length), binwidth=50) +
    geom_vline(xintercept=mean(telo_stats$vrr_telo_length), color='red') +
    xlab("VRR Telomere Length (Base Pairs)") + ylab("Read Frequency") +
    theme_minimal() +
    theme(axis.text = element_text(size=15),
          axis.title = element_text(size=20))
  
ggsave("vrr_length_hist.pdf", device="pdf", plot=vrr_length_hist, width=10, height=10)
  
vrr_length_boxplot <- ggplot(data=telo_stats) +
    geom_boxplot(mapping=aes(y=vrr_telo_length)) +
    ylab("VRR Telomere Length") +
    theme_minimal() + theme(axis.text.x = element_blank(),
                            axis.title.y = element_text(size=20),
                            axis.text = element_text(size=15))
  
ggsave("vrr_length_boxplot.pdf", device="pdf", plot=vrr_length_boxplot, width=5, height=10)
  
  # plot vrr length scatterplot by read length
vrr_length_scatter <- ggplot(data=telo_stats) +
    geom_pointdensity(mapping=aes(x=vrr_telo_length, y=read_len), adjust=4) +
    scale_color_viridis() +
    theme_minimal() +
    theme(axis.text = element_text(size=15),
          axis.title = element_text(size=20)) +
    xlab("VRR Telomere Length (Base Pairs)") +
    ylab("Read Length (Base Pairs)")
  
ggsave("vrr_length_by_read_length.pdf", device="pdf", plot=vrr_length_scatter, width=10, height=10)
  
  # plot vrr length grouped bar plot
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


#   # plot vrr boxplot comparing telo length to vrr length
#   telo_vrr_box <- ggplot(data = telo_stats) +
#     geom_boxplot(mapping=aes(x="Telo Length", y=telo_length)) +
#     geom_boxplot(mapping=aes(x="VRR Length", y=vrr_telo_length)) +
#     theme_minimal() +
#     theme(axis.title.x = element_blank(),
#           axis.text.x = element_text(angle=45),
#           axis.text = element_text(size=15),
#           axis.title = element_text(size=20))
#   
#   ggsave("vrr.telo_comparison.box.pdf", plot=telo_vrr_box, device="pdf", width=6, height=10)
#   
#   # plot vrr llength compared to telo length binned box plot
#   telo_stats$bin_telo_length <- unlist(lapply(telo_stats$telo_length, function(x) telo_lengths_for_binning[which.min(abs(telo_lengths_for_binning-x))]))
#   telo_stats$bin_telo_length <- factor(telo_stats$bin_telo_length, levels=c(10500, 9500, 8500, 7500, 6500, 5500, 4500, 3500, 2500, 1500, 500))
#   
#   telo_vrr_bar <- ggplot(data=telo_stats) +
#     geom_bar(mapping=aes(x="Telo Length", fill=bin_telo_length), position="fill") +
#     geom_bar(mapping=aes(x="VRR Length", fill=bin_vrr_length), position="fill") +
#     theme_minimal() +
#     theme(axis.title.x = element_blank(),
#           axis.text.x = element_text(angle=45),
#           axis.text = element_text(size=15),
#           axis.title = element_text(size=20))
#   ggsave("vrr.telo_comparison.bar.pdf", plot=telo_vrr_bar, device="pdf", width=6, height=10)
#   
# }


# if strand comparison - basically plot everything above again but now showing strand differences
if (args[2]) {
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
  
  ggsave("STRAND_COMPARISON/Read_Length.pdf", device="pdf", width=10, height=12, plot=read_length_by_strand, create.dir = TRUE)

    
  vrr_length_by_strand <- ggplot(data=telo_stats) +
      geom_histogram(mapping=aes(vrr_telo_length), binwidth=50) +
      facet_grid(strand ~ .) +
      geom_vline(data=filter(telo_stats, strand=="C"), aes(xintercept=mean(telo_stats[telo_stats$strand == "C",]$vrr_telo_length)), color='red', linetype="dashed") +
      geom_vline(data=filter(telo_stats, strand=="G"), aes(xintercept=mean(telo_stats[telo_stats$strand == "G",]$vrr_telo_length)), color='red', linetype="dashed") +
      xlab("VRR Telomere Length") + ylab("Telomere Frequency") +
      theme_minimal() +
      theme(axis.text = element_text(size=15),
            axis.title = element_text(size=20),
            panel.border = element_rect(color="black", fill=NA))
    
  ggsave("STRAND_COMPARISON/VRR_Telo_Length.pdf", device="pdf", width=10, height=12, plot=vrr_length_by_strand, create.dir = TRUE)
  
    
  vrr_length_by_strand_boxplot <- ggplot(data=telo_stats) +
      geom_boxplot(mapping=aes(x=strand,y=vrr_telo_length)) + theme_minimal() +
      theme(axis.title = element_text(size=20),
            axis.text=element_text(size=15)) +
      xlab("Strand") + ylab("VRR Telomere Length")
    
  ggsave("STRAND_COMPARISON/vrr_length_boxplot.pdf", device="pdf", width=8, height=12, plot=vrr_length_by_strand_boxplot)
    
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
    
  ggsave("STRAND_COMPARISON/vrr_telo_length_bar_plot.pdf", plot = vrr_bar_hist, device="pdf", width=9, height=10, create.dir = TRUE)
    
}

  

