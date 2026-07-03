#!/usr/bin/env Rscript

library(ggpubr)
library(dplyr)
library(ggplot2)
library(stringr)

args = commandArgs(trailingOnly=TRUE)
# args <- c("/media/data_01/ndeimler/WORKFLOWS/TARPON/test_data/concordance_sample_sheet.csv",
#           "/media/data_01/ndeimler/HUMAN/TARPON_LINEAGE_TRACKING/HG002.test.txt")

sample_df <- read.table(args[1], sep=",",header=TRUE)
# 
# for (i in seq_len(nrow(sample_df))) {
#   assign(sample_df$sample[i], read.table(sample_df$stats[i], sep="\t", header=TRUE))
# }


samples <- setNames(
  lapply(sample_df$stats, read.table, header = TRUE),
  sample_df$sample
)

concordance_df <- read.table(args[2], sep="\t", header=TRUE, check.names=FALSE)

results <- list()#vector("list", nrow(concordance_df))
for (concord_cluster in concordance_df$Concordance_Cluster) {
  results[[concord_cluster]] <- lapply(names(samples), function(samp) {
      cluster = str_split_i(concordance_df[concordance_df$Concordance_Cluster == concord_cluster, samp], "_", 2)
      subset(samples[[samp]], Cluster==cluster)
  })
 names(results[[concord_cluster]]) <- names(samples)
}


for (concordance_cluster in names(results)){
  
  plot_df <- do.call(rbind, lapply(names(results[[concordance_cluster]]),
                                   function(samp) {
                                     df <- results[[concordance_cluster]][[samp]]
                                     if (nrow(df) == 0) {
                                       return(NULL)
                                     }
                                     df$sample <- samp
                                     df
                                  }))
  
  comparisons <- combn(unique(plot_df$sample), 2, simplify = FALSE)
  plt <- ggplot(data=plot_df, mapping=aes(x=sample, y=vrr_telo_length/1000)) +
    geom_boxplot() +
    theme_minimal() + ylab("Telomere Length (kbp)") +
    stat_compare_means(comparisons=comparisons,method="wilcox.test", label="p.signif", size=5) +
    theme(axis.title.x=element_blank(),
          axis.title.y = element_text(size=11),
          axis.text = element_text(size=9),
          axis.text.x = element_text(angle=45, hjust=1))
  
  ggsave(paste(paste("BOXPLOTS/", concordance_cluster, sep=""), ".pdf", sep=""), plot=plt, width=8, height=6, units="cm", device="pdf")
 
  plt <- ggplot(data=plot_df, mapping=aes(x=sample, y=vrr_telo_length/1000)) +
    geom_violin() +
    theme_minimal() + ylab("Telomere Length (kbp)") +
    stat_compare_means(comparisons=comparisons,method="wilcox.test", label="p.signif", size=5) +
    theme(axis.title.x=element_blank(),
          axis.title.y = element_text(size=11),
          axis.text = element_text(size=9),
          axis.text.x = element_text(angle=45, hjust=1))
  ggsave(paste(paste("VIOLIN_PLOTS/", concordance_cluster, sep=""), ".pdf", sep=""), plot=plt, width=8, height=6, units="cm", device="pdf")
}


