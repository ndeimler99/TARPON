#!/usr/bin/env Rscript

# load appropriate libraries
library(ggplot2)
library(dplyr)
library(viridis)
library(ggpointdensity)

# get args args[1] = telo_stats (file), args[2] = plot_telo_length (boolean), args[3] = plot_vrr_length (boolean), args[4] = strand_comparison (boolean)
args = commandArgs(trailingOnly=TRUE)
# open telo stats file into dataframe
aln <- read.table(args[1], sep="\t", header=TRUE)
stats_df <- read.table(args[3], sep="\t", header=TRUE)

stats_df <- left_join(stats_df, aln, by="read_id")

summary_stats <- stats_df %>% group_by(Cluster) %>% summarize(count=n(), mean=mean(vrr_telo_length), sd=sd(vrr_telo_length),
                                                              q1=quantile(vrr_telo_length, 0.25), median=quantile(vrr_telo_length, 0.5),
                                                              q3=quantile(vrr_telo_length, 0.75))

clustering_summary <- stats_df %>% group_by(Cluster) %>% summarize(count=n())
clustering_summary$proportion <- clustering_summary$count / sum(clustering_summary$count) * 100
na_reads <- clustering_summary[is.na(clustering_summary$Cluster),]
na_reads$proportion <- as.numeric(na_reads$proportion)
# percentage of reads clustered
proportion_of_reads <- ggplot(data=na_reads) +
  geom_bar(mapping=aes(x=1, y=100-proportion), stat="identity") +
  ylim(0,100) +
  theme_minimal() +
  theme(axis.title=element_text(size=20),
        axis.text=element_text(size=15),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x=element_blank()) +
  ylab("Percentage of Reads\nAssigned to Cluster")



clustering_summary <- clustering_summary[!is.na(clustering_summary$Cluster),]
clustering_summary$proportion <- clustering_summary$count / sum(clustering_summary$count) * 100

# cluster size distribution box plot
cluster_sizes <- ggplot(data=clustering_summary) +
  geom_boxplot(mapping=aes(y=proportion)) +
  theme_minimal() +
  ylab("Percentage of Reads\nAssigned to Each Cluster") +
  theme(axis.text = element_text(size=15),
        axis.title=element_text(size=20),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank()) +
  geom_hline(mapping=aes(yintercept=1.08), linetype="dashed", color="red")


# number of clusters

cluster_count <- ggplot() + 
  geom_bar(mapping=aes(x=1, y=length(clustering_summary$Cluster)), stat="identity") +
  theme_minimal() +
  theme(axis.text = element_text(size=15),
        axis.title=element_text(size=20),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank()) +
  geom_hline(mapping=aes(yintercept=92), linetype="dashed", color="red") +
  ylab("Number of Clusters")

# telomere length by cluster boxplot

cluster_order <- stats_df %>% group_by(Cluster) %>% summarize(mean=mean(vrr_telo_length)) %>% arrange(mean)
cluster_order <- cluster_order[!is.na(cluster_order$Cluster),]
cluster_levels=c(NA, cluster_order$Cluster)

stats_df$Cluster <- factor(stats_df$Cluster, levels=cluster_levels)

cluster_by_telo_length <- ggplot(data=stats_df) +
  geom_boxplot(mapping=aes(x=Cluster, y=vrr_telo_length)) +
  theme_minimal() +
  theme(axis.title=element_text(size=20),
        axis.text=element_text(size=15)) +
  ylab("VRR Telomere Length") 


if (grepl("cluster", args[1])){
  write.table(summary_stats, paste(args[2], ".clustering_summary_stats.txt", sep=""), sep="\t")
  ggsave("cluster_by_telo_length.pdf", width=20, height=10, plot=cluster_by_telo_length)  
  ggsave("cluster_counts.pdf", width=6, height=10, plot=cluster_count)
  ggsave("cluster_sizes_boxplot.pdf", width=6, height=10, plot=cluster_sizes)
  ggsave("reads_assigned_to_clusters.pdf", width=6, height=10, plot=proportion_of_reads)
} else {
  write.table(summary_stats, paste(args[2], ".chrom_arm_summary.txt", sep=""), sep="\t")
  ggsave("chrom_by_telo_length.pdf", width=20, height=10, plot=cluster_by_telo_length)  
  ggsave("chrom_counts.pdf", width=6, height=10, plot=cluster_count)
  ggsave("chrom_sizes_boxplot.pdf", width=6, height=10, plot=cluster_sizes)
  ggsave("reads_assigned_to_chrom.pdf", width=6, height=10, plot=proportion_of_reads)
}