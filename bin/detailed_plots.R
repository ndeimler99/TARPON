#!/usr/bin/env Rscript
library(ggplot2)
library(dplyr)
library(viridis)
library(ggpointdensity)
library(ggpubr)

args = commandArgs(trailingOnly=TRUE)
telo_stats <- read.table(args[1], header=TRUE)
telo_lengths_for_binning <- c(10500, 9500, 8500, 7500, 6500, 5500, 4500, 3500, 2500, 1500, 500)

# histogram perfect repeats
perfect_repeat_hist <- ggplot() +
  geom_histogram(mapping=aes(telo_stats$perc_perfect), binwidth = 5) + 
  theme_minimal() +
  ylab("Number of Reads") + xlab("% Telomeric (Perfect Repeats)") +
  theme(axis.text = element_text(size=15),
        axis.title = element_text(size=20))

# histogram one subs
imperfect_repeat_hist <- ggplot() +
  geom_histogram(mapping=aes(telo_stats$perc_variant), binwidth = 5) + 
  theme_minimal() +
  ylab("Number of Reads") + xlab("% Telomeric (<= One Nucl. Substitutions Repeats)") +
  theme(axis.text = element_text(size=15),
        axis.title = element_text(size=20))

repeat_histograms <- ggarrange(perfect_repeat_hist, imperfect_repeat_hist, nrow=2)
ggsave("DETAILED_STATS/percentage_telomere_histogram.pdf", width=12, height=10, plot=repeat_histograms, device="pdf", create.dir = TRUE)


percentage_repeat_box <- ggplot() +
  geom_boxplot(mapping=aes(y="Perfect Repeats", x=telo_stats$perc_perfect)) +
  geom_boxplot(mapping=aes(y="One Nucleotide\nSubstitutions", x=telo_stats$perc_variant)) +
  theme_minimal() +
  xlab("Percentage of VRR Region") +
  theme(axis.title.y= element_blank(),
        axis.text = element_text(size=15),
        axis.title.x = element_text(size=20))

ggsave("DETAILED_STATS/percentage_telomere_boxplot.pdf", width=12, height=6, device="pdf", plot=percentage_repeat_box, create.dir = TRUE)


perfect_by_imperfect_scatter <- ggplot() +
  geom_pointdensity(mapping=aes(x=telo_stats$perc_perfect, y=telo_stats$perc_variant), adjust=4) +
  scale_color_viridis() +
  theme_minimal() +
  xlab("% VRR Region Perfect Repeats") + 
  ylab("% VRR Region of One Nucl. Substitution Repeats") +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20))

ggsave("DETAILED_STATS/pecentage_telomeric_vs_one_nucl_subs.pdf", width=12, height=10, plot=perfect_by_imperfect_scatter, device="pdf", create.dir = TRUE)


read_quality_histogram <- ggplot() +
  geom_histogram(mapping = aes(telo_stats$read_quality), binwidth = 1) +
  geom_vline(xintercept=mean(telo_stats$read_quality), color="red") +
  theme_minimal() + ylab("Frequency") + xlab("Read Quality (Phred Score)") +
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=15))

telo_quality_histogram <- ggplot() +
  geom_histogram(mapping = aes(telo_stats$telo_quality), binwidth = 1) +
  geom_vline(xintercept=mean(telo_stats$telo_quality), color="red") +
  theme_minimal() + ylab("Frequency") + xlab("Telomere Quality (Phred Score)") +
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=15))

quality_histograms <- ggarrange(read_quality_histogram, telo_quality_histogram, nrow=2)

ggsave("DETAILED_STATS/quality_score_histograms.pdf", device="pdf", plot=quality_histograms, width=12, height=10, create.dir = TRUE)

quality_boxplots <- ggplot() +
  geom_boxplot(mapping=aes(y="Read Quality", x=telo_stats$read_quality)) +
  geom_boxplot(mapping=aes(y="Telomere Quality", x=telo_stats$telo_quality)) +
  xlab("Quality Score") + theme_minimal() +
  theme(axis.text=element_text(size=15),
        axis.title.x = element_text(size=20),
        axis.title.y = element_blank())

ggsave("DETAILED_STATS/quality_score_boxplots.pdf", device="pdf", plot=quality_boxplots, width=12, height=6, create.dir = TRUE)

read_vs_telo_quality <- ggplot() +
  geom_pointdensity(mapping=aes(x=telo_stats$read_quality, y=telo_stats$telo_quality), adjust=4) +
  scale_color_viridis() +
  theme_minimal() +
  xlab("Read Quality Score") + 
  ylab("Telomere Quality Score") +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20))

ggsave("DETAILED_STATS/read_vs_telo_quality_scatter.pdf", width=12, height=10, plot=read_vs_telo_quality, device="pdf", create.dir = TRUE)


vrr_vs_telo_length <- ggplot() +
  geom_pointdensity(mapping=aes(x=telo_stats$vrr_telo_length, y=telo_stats$telo_length), adjust=4) +
  scale_color_viridis() +
  theme_minimal() +
  xlab("VRR Telomere Length (BP)") + 
  ylab("Telomere Length (BP)") +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20))

ggsave("DETAILED_STATS/vrr_vs_telo_length_scatter.pdf", width=12, height=10, plot=vrr_vs_telo_length, device="pdf", create.dir = TRUE)

telo_quality_vs_perfect <- ggplot() +
  geom_pointdensity(mapping=aes(x=telo_stats$telo_quality, y=telo_stats$perc_perfect), adjust=4) +
  scale_color_viridis() +
  theme_minimal() +
  xlab("Telomere Quality") + 
  ylab("% VRR of Perfect Repeats") +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20))

ggsave("DETAILED_STATS/telomere_quality_vs_percentage_perfect.pdf", width=12, height=10, plot=telo_quality_vs_perfect, device="pdf", create.dir = TRUE)

telo_quality_vs_imperfect <- ggplot() +
  geom_pointdensity(mapping=aes(x=telo_stats$telo_quality, y=telo_stats$perc_variant), adjust=4) +
  scale_color_viridis() +
  theme_minimal() +
  xlab("Telomere Quality") + 
  ylab("% VRR of One Nucl. Substitution Repeats") +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20))

ggsave("DETAILED_STATS/telomere_quality_vs_percentage_imperfect.pdf", width=12, height=10, plot=telo_quality_vs_imperfect, device="pdf", create.dir = TRUE)

if (args[2]){
  # plot using telo lengths
  
  # scatterplot telo length by percentage
  telo_length_vs_perfect <- ggplot() +
    geom_pointdensity(mapping=aes(x=telo_stats$telo_length, y=telo_stats$perc_perfect), adjust=4) +
    scale_color_viridis() +
    theme_minimal() +
    xlab("Telomere Length") + 
    ylab("% VRR of Perfect Repeats") +
    theme(axis.text=element_text(size=15),
          axis.title=element_text(size=20))
  telo_length_vs_perfect
  ggsave("DETAILED_STATS/telomere_length_vs_percentage_perfect.pdf", width=12, height=10, plot=telo_length_vs_perfect, device="pdf", create.dir = TRUE)
  
  # scatterplot telo length percentage one subs
  telo_length_vs_imperfect <- ggplot() +
    geom_pointdensity(mapping=aes(x=telo_stats$telo_length, y=telo_stats$perc_variant), adjust=4) +
    scale_color_viridis() +
    theme_minimal() +
    xlab("Telomere Length") + 
    ylab("% VRR of One Nucl. Substitution Repeats") +
    theme(axis.text=element_text(size=15),
          axis.title=element_text(size=20))

    ggsave("DETAILED_STATS/telomere_length_vs_percentage_imperfect.pdf", width=12, height=10, plot=telo_length_vs_imperfect, device="pdf", create.dir = TRUE)
  
    # quality by length
    telo_length_by_telo_quality <- ggplot() +
      geom_pointdensity(mapping=aes(x=telo_stats$telo_length, y=telo_stats$telo_quality), adjust=4) +
      scale_color_viridis() +
      theme_minimal() +
      xlab("Telomere Length") + 
      ylab("Telomere Quality") +
      theme(axis.text=element_text(size=15),
            axis.title=element_text(size=20))
    telo_length_by_telo_quality
    ggsave("DETAILED_STATS/telomere_length_vs_telomere_quality.pdf", width=12, height=10, plot=telo_length_by_telo_quality, device="pdf", create.dir = TRUE)
  
    telo_length_by_read_quality <- ggplot() +
      geom_pointdensity(mapping=aes(x=telo_stats$telo_length, y=telo_stats$read_quality), adjust=4) +
      scale_color_viridis() +
      theme_minimal() +
      xlab("Telomere Length") + 
      ylab("Read Quality") +
      theme(axis.text=element_text(size=15),
            axis.title=element_text(size=20))
    ggsave("DETAILED_STATS/telomere_length_vs_read_quality.pdf", width=12, height=10, plot=telo_length_by_read_quality, device="pdf", create.dir = TRUE)  
}

if (args[3]){
  # plot using VRR length

  # scatterplot telo length by percentage
  vrr_length_vs_perfect <- ggplot() +
    geom_pointdensity(mapping=aes(x=telo_stats$vrr_telo_length, y=telo_stats$perc_perfect), adjust=4) +
    scale_color_viridis() +
    theme_minimal() +
    xlab("VRR Telomere Length") + 
    ylab("% VRR of Perfect Repeats") +
    theme(axis.text=element_text(size=15),
          axis.title=element_text(size=20))
  
  ggsave("DETAILED_STATS/vrr_telomere_length_vs_percentage_perfect.pdf", width=12, height=10, plot=vrr_length_vs_perfect, device="pdf", create.dir = TRUE)
  
  # scatterplot telo length percentage one subs
  vrr_length_vs_imperfect <- ggplot() +
    geom_pointdensity(mapping=aes(x=telo_stats$vrr_telo_length, y=telo_stats$perc_variant), adjust=4) +
    scale_color_viridis() +
    theme_minimal() +
    xlab("VRR Telomere Length") + 
    ylab("% VRR of One Nucl. Substitution Repeats") +
    theme(axis.text=element_text(size=15),
          axis.title=element_text(size=20))
  
  ggsave("DETAILED_STATS/vrr_telomere_length_vs_percentage_imperfect.pdf", width=12, height=10, plot=vrr_length_vs_imperfect, device="pdf", create.dir = TRUE)
  
  # quality by length
  vrr_length_by_telo_quality <- ggplot() +
    geom_pointdensity(mapping=aes(x=telo_stats$vrr_telo_length, y=telo_stats$telo_quality), adjust=4) +
    scale_color_viridis() +
    theme_minimal() +
    xlab("VRR Telomere Length") + 
    ylab("Telomere Quality") +
    theme(axis.text=element_text(size=15),
          axis.title=element_text(size=20))
  
  ggsave("DETAILED_STATS/vrr_telomere_length_vs_telomere_quality.pdf", width=12, height=10, plot=vrr_length_by_telo_quality, device="pdf", create.dir = TRUE)
  
  vrr_length_by_read_quality <- ggplot() +
    geom_pointdensity(mapping=aes(x=telo_stats$vrr_telo_length, y=telo_stats$read_quality), adjust=4) +
    scale_color_viridis() +
    theme_minimal() +
    xlab("VRR Telomere Length") + 
    ylab("Read Quality") +
    theme(axis.text=element_text(size=15),
          axis.title=element_text(size=20))
  ggsave("DETAILED_STATS/vrr_telomere_length_vs_read_quality.pdf", width=12, height=10, plot=vrr_length_by_read_quality, device="pdf", create.dir = TRUE)  
  
}

if (args[4]) {
  # plot c-g strand comparison
  
  # % telo histogram
  telomere_perc_boxplot <- ggplot() +
    geom_boxplot(mapping=aes(y="C", x=telo_stats[telo_stats$strand=="C",]$perc_perfect)) +
    geom_boxplot(mapping=aes(y="G", x=telo_stats[telo_stats$strand=="G",]$perc_perfect)) +
    theme_minimal() + ylab("Strand") + xlab("Percent Perfect Telomeric") +
    theme(axis.text = element_text(size=15),
          axis.title = element_text(size=20))
  ggsave("C_G_COMPARISON/telomere_percentage_perfect_repeats.pdf", device="pdf", plot=telomere_perc_boxplot, width=12, height=6, create.dir = TRUE)
  # % one subs histogram
  telomere_imperfect_perc_boxplot <- ggplot() +
    geom_boxplot(mapping=aes(y="C", x=telo_stats[telo_stats$strand=="C",]$perc_variant)) +
    geom_boxplot(mapping=aes(y="G", x=telo_stats[telo_stats$strand=="G",]$perc_variant)) +
    theme_minimal() + ylab("Strand") + xlab("Percent One Nucleotide Substitutions Telomeric") +
    theme(axis.text = element_text(size=15),
          axis.title = element_text(size=20))
  ggsave("C_G_COMPARISON/telomere_percentage_imperfect_repeats.pdf", device="pdf", plot=telomere_imperfect_perc_boxplot, width=12, height=6, create.dir = TRUE)
  # quality differences
  quality_boxplot <- ggplot() +
    geom_boxplot(mapping=aes(y="C", x=telo_stats[telo_stats$strand=="C",]$telo_quality)) +
    geom_boxplot(mapping=aes(y="G", x=telo_stats[telo_stats$strand=="G",]$telo_quality)) +
    theme_minimal() + ylab("Strand") + xlab("Telomere Quality") +
    theme(axis.text = element_text(size=15),
          axis.title = element_text(size=20))
  ggsave("C_G_COMPARISON/telomere_quality_boxplot.pdf", device="pdf", plot=quality_boxplot, width=12, height=6, create.dir = TRUE)
  
  read_quality_boxplot <- ggplot() +
    geom_boxplot(mapping=aes(y="C", x=telo_stats[telo_stats$strand=="C",]$read_quality)) +
    geom_boxplot(mapping=aes(y="G", x=telo_stats[telo_stats$strand=="G",]$read_quality)) +
    theme_minimal() + ylab("Strand") + xlab("Read Quality") +
    theme(axis.text = element_text(size=15),
          axis.title = element_text(size=20))
  ggsave("C_G_COMPARISON/telomere_quality_boxplot.pdf", device="pdf", plot=read_quality_boxplot, width=12, height=6, create.dir = TRUE)
  
}


