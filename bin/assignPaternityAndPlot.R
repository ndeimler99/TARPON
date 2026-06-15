#!/usr/bin/env Rscript

# load appropriate libraries
library(ggplot2)
library(dplyr)
library(stringr)

# get args args[1] = telo_stats (file), args[2] = sample name
args = commandArgs(trailingOnly=TRUE)

# args <- c("/media/data_01/ndeimler/WORKFLOWS/TARPON/work/8a/91b3b18860da67cbaa79e22d9b228c/alignment_stats.txt",
#           "/media/data_01/ndeimler/WORKFLOWS/TARPON/work/7f/156751488fbb93978de1fe74d60739/alignment_stats.txt",
#           "/media/data_01/ndeimler/DKC_PROJECT/DKC10/146_DKC10_ROBOT/146_0.7.0_TARPON/sample/sample.telo_stats_with_cluster.txt")
#           
maternal_aln <- read.table(args[1], sep="\t", header=TRUE)
paternal_aln <- read.table(args[2], sep="\t", header=TRUE)

# open telo stats file into dataframe
offspring_stats <- read.table(args[3], sep="\t", header=TRUE)

maternal_bests <- maternal_aln %>% group_by(offspring_cluster) %>% slice_max(p_ident, n=1, with_ties = FALSE) %>% ungroup()
maternal_bests <- maternal_bests[maternal_bests$p_ident > 85,]
paternal_bests <- paternal_aln %>% group_by(offspring_cluster) %>% slice_max(p_ident, n=1, with_ties = FALSE) %>% ungroup()
paternal_bests <- paternal_bests[paternal_bests$p_ident > 85,]

merged_df <- left_join(maternal_bests, paternal_bests, by="offspring_cluster", suffix=c(".maternal", ".paternal"))

merged_df$paternity <- ifelse(is.na(merged_df$p_ident.maternal) & !is.na(merged_df$p_ident.paternal), "Paternal",
                               ifelse(is.na(merged_df$p_ident.paternal) & !is.na(merged_df$p_ident.maternal), "Maternal",
                                      ifelse(abs(merged_df$p_ident.maternal - merged_df$p_ident.paternal) <= 1, NA,
                                             ifelse(merged_df$p_ident.paternal > merged_df$p_ident.maternal, "Paternal", "Maternal"))))


write.table(merged_df, "cluster_assignment_stats.txt", sep="\t", quote=FALSE)

assignment_df <- merged_df[,c("offspring_cluster", "paternity")]
assignment_df$offspring_cluster <- str_split_i(assignment_df$offspring_cluster, "_", 2)
assignment_df$offspring_cluster <- as.integer(assignment_df$offspring_cluster)
offspring_stats <- left_join(offspring_stats, assignment_df, by=c("Cluster"="offspring_cluster"))

write.table(offspring_stats, "offspring_telo_stats.txt", sep="\t", quote=FALSE)


# telomere length by cluster sorted from left to right on increasing telomere length colored by inheritance

order_df <- offspring_stats %>% group_by(Cluster) %>% 
                                summarise(telo_median=median(vrr_telo_length), .groups="drop") %>%
                              arrange(telo_median) %>% mutate(Cluster=factor(Cluster, levels=Cluster))

offspring_stats$Cluster <- factor(offspring_stats$Cluster, levels=levels(order_df$Cluster))

plt <- ggplot(data=offspring_stats) +
  geom_boxplot(mapping=aes(x=Cluster, y=vrr_telo_length/1000, fill=paternity)) +
  theme_minimal() + xlab("Cluster") + ylab("Telomere Length [kbp]") +
  theme(axis.text.x=element_text(size=10, angle=90, hjust=1, vjust=0.5),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12)) +
  scale_fill_manual(breaks=c("Maternal", "Paternal", NA), values=c("#D81B60", "#1E88E5", "gray"), name="Paternity")

ggsave("telomere_length_by_cluster_labeled.pdf", plot=plt, width=30, height=10, units="cm")

# telomere length by cluster sorted from left to right on increasing telomere length segregated and colored by inheritance

order_df <- offspring_stats %>% group_by(Cluster, paternity) %>% 
  summarise(telo_median=median(vrr_telo_length), .groups="drop") %>%
  arrange(paternity, telo_median) %>% mutate(Cluster=factor(Cluster, levels=Cluster))

offspring_stats$Cluster <- factor(offspring_stats$Cluster, levels=levels(order_df$Cluster))

plt <- ggplot(data=offspring_stats) +
  geom_boxplot(mapping=aes(x=Cluster, y=vrr_telo_length/1000, fill=paternity)) +
  theme_minimal() + xlab("Cluster") + ylab("Telomere Length [kbp]") +
  theme(axis.text.x=element_text(size=10, angle=90, hjust=1, vjust=0.5),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12)) +
  scale_fill_manual(breaks=c("Maternal", "Paternal", NA), values=c("#D81B60", "#1E88E5", "gray"), name="Paternity")

ggsave("telomere_length_by_cluster_labeled_and_segregated.pdf", plot=plt, width=30, height=10, units="cm")

# boxplot summary of cluster median telomere lengths comparing maternal to paternal telomeres

summary_df <- offspring_stats %>% group_by(Cluster, paternity) %>% summarize(median_telo=median(vrr_telo_length))
plt <- ggplot(data=summary_df) +
  geom_boxplot(mapping=aes(x=paternity, y=median_telo/1000)) +
  geom_jitter(mapping=aes(x=paternity, y=median_telo/1000, color=paternity), width=0.2) +
  theme_minimal() + ylab("Cluster Median\nTelomere Length [kbp]") + 
  theme(axis.text.x=element_text(size=10, angle=45, hjust=1),
        axis.title.x=element_blank(),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12)) +
  scale_color_manual(breaks=c("Maternal", "Paternal", NA), values=c("#D81B60", "#1E88E5", "gray"), name="Paternity")
ggsave("cluster_median_telo_length.pdf", plot=plt, width=10, height=10, units="cm")

