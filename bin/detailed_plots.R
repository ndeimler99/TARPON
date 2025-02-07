#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)
library(viridis)
library(ggpointdensity)
library(ggpubr)

# script aims to create detailed figures from telomeric statistics

args = commandArgs(trailingOnly=TRUE)
telo_stats <- read.table(args[1], header=TRUE)
TELO_LENGTHS_FOR_BINNING <- c(10500, 9500, 8500, 7500, 6500, 5500, 4500, 3500, 2500, 1500, 500)


# read quality histogram
read_quality_histogram <- ggplot() +
  geom_histogram(mapping = aes(telo_stats$read_qual), binwidth = 1) +
  geom_vline(xintercept=mean(telo_stats$read_qual), color="red", linetype="dashed") +
  theme_minimal() + ylab("Frequency") + xlab("Read Quality (Phred Score)") +
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=15))

ggsave("read_quality_histogram.pdf", device="pdf", width=12, height=10, plot=read_quality_histogram)

# telo quality histogram
telo_quality_histogram <- ggplot() +
  geom_histogram(mapping = aes(telo_stats$telo_qual), binwidth = 1) +
  geom_vline(xintercept=mean(telo_stats$telo_qual), color="red", linetype="dashed") +
  theme_minimal() + ylab("Frequency") + xlab("Telomere Quality (Phred Score)") +
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=15))

ggsave("vrr_telo_quality_histogram.pdf", device="pdf", width=12, height=10, plot=telo_quality_histogram)

# telo vs read quality boxplots
quality_boxplots <- ggplot() +
  geom_boxplot(mapping=aes(y="Read Quality", x=telo_stats$read_qual)) +
  geom_boxplot(mapping=aes(y="Telomere Quality", x=telo_stats$telo_qual)) +
  xlab("Quality Score") + theme_minimal() +
  theme(axis.text=element_text(size=15),
        axis.title.x = element_text(size=20),
        axis.title.y = element_blank())

ggsave("read_vs_vrr_telo_quality_score_boxplots.pdf", device="pdf", plot=quality_boxplots, width=12, height=6, create.dir = TRUE)

# telo vs read quality scatterplot
read_vs_telo_quality <- ggplot() +
  geom_pointdensity(mapping=aes(x=telo_stats$read_qual, y=telo_stats$telo_qual), adjust=4) +
  scale_color_viridis() +
  theme_minimal() +
  xlab("Read Quality Score") + 
  ylab("Telomere Quality Score") +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20))

ggsave("read_vs_vrr_telo_quality_scatter.pdf", width=12, height=10, plot=read_vs_telo_quality, device="pdf", create.dir = TRUE)

# read quality by telo length (and vrr length)
# telo quality by telo length (and vrr length)
if (args[2]){
  # quality by length
  telo_length_by_telo_quality <- ggplot() +
    geom_pointdensity(mapping=aes(x=telo_stats$telo_length, y=telo_stats$telo_qual), adjust=4) +
    scale_color_viridis() +
    theme_minimal() +
    xlab("Telomere Length") + 
    ylab("Telomere Quality") +
    theme(axis.text=element_text(size=15),
          axis.title=element_text(size=20))
  telo_length_by_telo_quality
  ggsave("telomere_length_vs_telomere_quality.pdf", width=12, height=10, plot=telo_length_by_telo_quality, device="pdf", create.dir = TRUE)
  
  telo_length_by_read_quality <- ggplot() +
    geom_pointdensity(mapping=aes(x=telo_stats$telo_length, y=telo_stats$read_qual), adjust=4) +
    scale_color_viridis() +
    theme_minimal() +
    xlab("Telomere Length") + 
    ylab("Read Quality") +
    theme(axis.text=element_text(size=15),
          axis.title=element_text(size=20))
  ggsave("telomere_length_vs_read_quality.pdf", width=12, height=10, plot=telo_length_by_read_quality, device="pdf", create.dir = TRUE)  
}

if (args[3]){
  # quality by length
  vrr_length_by_telo_quality <- ggplot() +
    geom_pointdensity(mapping=aes(x=telo_stats$vrr_telo_length, y=telo_stats$telo_qual), adjust=4) +
    scale_color_viridis() +
    theme_minimal() +
    xlab("VRR Telomere Length") + 
    ylab("Telomere Quality") +
    theme(axis.text=element_text(size=15),
          axis.title=element_text(size=20))
  
  ggsave("vrr_telomere_length_vs_telomere_quality.pdf", width=12, height=10, plot=vrr_length_by_telo_quality, device="pdf", create.dir = TRUE)
  
  vrr_length_by_read_quality <- ggplot() +
    geom_pointdensity(mapping=aes(x=telo_stats$vrr_telo_length, y=telo_stats$read_qual), adjust=4) +
    scale_color_viridis() +
    theme_minimal() +
    xlab("VRR Telomere Length") + 
    ylab("Read Quality") +
    theme(axis.text=element_text(size=15),
          axis.title=element_text(size=20))
  ggsave("vrr_telomere_length_vs_read_quality.pdf", width=12, height=10, plot=vrr_length_by_read_quality, device="pdf", create.dir = TRUE)  
  
}

if (args[4]){
  # read quality by strand
  quality_boxplot <- ggplot() +
    geom_boxplot(mapping=aes(y="C", x=telo_stats[telo_stats$strand=="C",]$telo_qual)) +
    geom_boxplot(mapping=aes(y="G", x=telo_stats[telo_stats$strand=="G",]$telo_qual)) +
    theme_minimal() + ylab("Strand") + xlab("Telomere Quality") +
    theme(axis.text = element_text(size=15),
          axis.title = element_text(size=20))
  ggsave("telomere_quality_boxplot_by_strand.pdf", device="pdf", plot=quality_boxplot, width=12, height=6, create.dir = TRUE)
  
  # telo quality by strand
  
  read_quality_boxplot <- ggplot() +
    geom_boxplot(mapping=aes(y="C", x=telo_stats[telo_stats$strand=="C",]$read_qual)) +
    geom_boxplot(mapping=aes(y="G", x=telo_stats[telo_stats$strand=="G",]$read_qual)) +
    theme_minimal() + ylab("Strand") + xlab("Read Quality") +
    theme(axis.text = element_text(size=15),
          axis.title = element_text(size=20))
  ggsave("read_quality_boxplot_by_strand.pdf", device="pdf", plot=read_quality_boxplot, width=12, height=6, create.dir = TRUE)
}

# distribution of perfect telomeric repeats
perfect_repeat_hist <- ggplot() +
  geom_histogram(mapping=aes(telo_stats$wt_composition), binwidth = 5) + 
  theme_minimal() +
  ylab("Number of Reads") + xlab("% Telomeric (Perfect Repeats)") +
  theme(axis.text = element_text(size=15),
        axis.title = element_text(size=20))

ggsave("telo_repeat_frequency.pdf", plot=perfect_repeat_hist, device="pdf", width=12, height=10)

# boxplot of perfect telomeric repeats
percentage_repeat_box <- ggplot() +
  geom_boxplot(mapping=aes(y="Perfect Repeats", x=telo_stats$wt_composition)) +
  geom_boxplot(mapping=aes(y="One Nucleotide\nSubstitutions", x=telo_stats$one_nucl_variant_composition)) +
  theme_minimal() +
  xlab("Percentage of VRR Region") +
  theme(axis.title.y= element_blank(),
        axis.text = element_text(size=15),
        axis.title.x = element_text(size=20))

ggsave("telo_repeat_frequency_boxplot.pdf", width=6, height=10, device="pdf", plot=percentage_repeat_box, create.dir = TRUE)

# scatterplot telomere quality by perfect repeat composition
telo_quality_vs_perfect <- ggplot() +
  geom_pointdensity(mapping=aes(x=telo_stats$telo_qual, y=telo_stats$wt_composition), adjust=4) +
  scale_color_viridis() +
  theme_minimal() +
  xlab("Telomere Quality") + 
  ylab("% VRR of Perfect Repeats") +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20))

ggsave("telomere_quality_vs_percentage_repeats.pdf", width=12, height=10, plot=telo_quality_vs_perfect, device="pdf", create.dir = TRUE)


# distribution of telomer-like repeats
imperfect_repeat_hist <- ggplot() +
  geom_histogram(mapping=aes(telo_stats$perc_variant), binwidth = 5) + 
  theme_minimal() +
  ylab("Number of Reads") + xlab("% Telomeric (<= One Nucl. Substitutions Repeats)") +
  theme(axis.text = element_text(size=15),
        axis.title = element_text(size=20))

ggsave("telomere_like_repeat_frequency.pdf", width=12, height=10, plot=imperfect_repeat_hist, device="pdf", create.dir = TRUE)


# scatterplot telomere quality by telomer-like composition
telo_quality_vs_imperfect <- ggplot() +
  geom_pointdensity(mapping=aes(x=telo_stats$telo_qual, y=telo_stats$one_nucl_variant_composition), adjust=4) +
  scale_color_viridis() +
  theme_minimal() +
  xlab("Telomere Quality") + 
  ylab("% VRR of One Nucl. Substitution Repeats") +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20))

ggsave("telomere_quality_vs_percentage_telomere_like_repeats.pdf", width=12, height=10, plot=telo_quality_vs_imperfect, device="pdf", create.dir = TRUE)


# perfect vs telomer-like repeats scatterpot and boxplot

perfect_by_imperfect_scatter <- ggplot() +
  geom_pointdensity(mapping=aes(x=telo_stats$wt_composition, y=telo_stats$one_nucl_variant_composition), adjust=4) +
  scale_color_viridis() +
  theme_minimal() +
  xlab("% VRR Region Perfect Repeats") + 
  ylab("% VRR Region of One Nucl. Substitution Repeats") +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20))

ggsave("telomere_vs_telomere_like_repeat_frequency_scatter.pdf", width=12, height=10, plot=perfect_by_imperfect_scatter, device="pdf", create.dir = TRUE)

# telo length by vrr length scatter plot

vrr_vs_telo_length <- ggplot() +
  geom_pointdensity(mapping=aes(x=telo_stats$vrr_telo_length, y=telo_stats$telo_length), adjust=4) +
  scale_color_viridis() +
  theme_minimal() +
  xlab("VRR Telomere Length (BP)") + 
  ylab("Telomere Length (BP)") +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20))

ggsave("vrr_vs_telo_length_scatter.pdf", width=12, height=10, plot=vrr_vs_telo_length, device="pdf", create.dir = TRUE)



if (args[2]){
  # plot using telo lengths
  
  # scatterplot telo length by percentage
  telo_length_vs_perfect <- ggplot() +
    geom_pointdensity(mapping=aes(x=telo_stats$telo_length, y=telo_stats$wt_composition), adjust=4) +
    scale_color_viridis() +
    theme_minimal() +
    xlab("Telomere Length") + 
    ylab("% VRR of Perfect Repeats") +
    theme(axis.text=element_text(size=15),
          axis.title=element_text(size=20))
  telo_length_vs_perfect
  ggsave("telomere_length_vs_percentage_repeats_scatteer.pdf", width=12, height=10, plot=telo_length_vs_perfect, device="pdf", create.dir = TRUE)
  
  # scatterplot telo length percentage one subs
  telo_length_vs_imperfect <- ggplot() +
    geom_pointdensity(mapping=aes(x=telo_stats$telo_length, y=telo_stats$one_nucl_variant_composition), adjust=4) +
    scale_color_viridis() +
    theme_minimal() +
    xlab("Telomere Length") + 
    ylab("% VRR of One Nucl. Substitution Repeats") +
    theme(axis.text=element_text(size=15),
          axis.title=element_text(size=20))

    ggsave("telomere_length_vs_percentage_imperfect_scatter.pdf", width=12, height=10, plot=telo_length_vs_imperfect, device="pdf", create.dir = TRUE)
  
}

if (args[3]){
  # plot using VRR length

  # scatterplot telo length by percentage
  vrr_length_vs_perfect <- ggplot() +
    geom_pointdensity(mapping=aes(x=telo_stats$vrr_telo_length, y=telo_stats$wt_composition), adjust=4) +
    scale_color_viridis() +
    theme_minimal() +
    xlab("VRR Telomere Length") + 
    ylab("% VRR of Perfect Repeats") +
    theme(axis.text=element_text(size=15),
          axis.title=element_text(size=20))
  
  ggsave("vrr_telomere_length_vs_percentage_perfect_scatter.pdf", width=12, height=10, plot=vrr_length_vs_perfect, device="pdf", create.dir = TRUE)
  
  # scatterplot telo length percentage one subs
  vrr_length_vs_imperfect <- ggplot() +
    geom_pointdensity(mapping=aes(x=telo_stats$vrr_telo_length, y=telo_stats$one_nucl_variant_composition), adjust=4) +
    scale_color_viridis() +
    theme_minimal() +
    xlab("VRR Telomere Length") + 
    ylab("% VRR of One Nucl. Substitution Repeats") +
    theme(axis.text=element_text(size=15),
          axis.title=element_text(size=20))
  
  ggsave("vrr_telomere_length_vs_percentage_imperfect_scatter.pdf", width=12, height=10, plot=vrr_length_vs_imperfect, device="pdf", create.dir = TRUE)
  
}

if (args[4]) {
  # plot c-g strand comparison
  
  # % telo histogram
  telomere_perc_boxplot <- ggplot() +
    geom_boxplot(mapping=aes(y="C", x=telo_stats[telo_stats$strand=="C",]$wt_composition)) +
    geom_boxplot(mapping=aes(y="G", x=telo_stats[telo_stats$strand=="G",]$wt_composition)) +
    theme_minimal() + ylab("Strand") + xlab("Percent Perfect Telomeric") +
    theme(axis.text = element_text(size=15),
          axis.title = element_text(size=20))
  ggsave("repeat_composition_by_strand.pdf", device="pdf", plot=telomere_perc_boxplot, width=12, height=6, create.dir = TRUE)
  # % one subs histogram
  telomere_imperfect_perc_boxplot <- ggplot() +
    geom_boxplot(mapping=aes(y="C", x=telo_stats[telo_stats$strand=="C",]$one_nucl_variant_composition)) +
    geom_boxplot(mapping=aes(y="G", x=telo_stats[telo_stats$strand=="G",]$one_nucl_variant_composition)) +
    theme_minimal() + ylab("Strand") + xlab("Percent One Nucleotide Substitutions Telomeric") +
    theme(axis.text = element_text(size=15),
          axis.title = element_text(size=20))
  ggsave("telomere-like_repeat_composition_by_strand.pdf", device="pdf", plot=telomere_imperfect_perc_boxplot, width=12, height=6, create.dir = TRUE)
}


