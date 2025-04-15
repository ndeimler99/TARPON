#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)
library(stringr)

COLORS=c('#0173b2', '#de8f05', '#029e73', '#d55e00', '#cc78bc', '#ca9161', '#fbafe4', '#949494', '#ece133', '#56b4e9')
TELO_LENGTHS_FOR_BINNING <- c(10500, 9500, 8500, 7500, 6500, 5500, 4500, 3500, 2500, 1500, 500)

args = commandArgs(trailingOnly=TRUE)

vrr_master_df <- c()

for (arg in args[1:length(args)]){
  df <- read.table(arg, header=TRUE)
  df$binned_telo_length <- unlist(lapply(df$vrr_telo_length, function(x) TELO_LENGTHS_FOR_BINNING[which.min(abs(TELO_LENGTHS_FOR_BINNING-x))]))
  binned_telo_length <- as.data.frame(table(df$binned_telo_length))
  binned_telo_length$sample <- str_split_i(str_split_i(arg, "\\.", 1), "/", -1)
  binned_telo_length$Freq <- binned_telo_length$Freq / sum(binned_telo_length$Freq) * 100
  vrr_master_df <- rbind(vrr_master_df, binned_telo_length)
}

vrr_master_df$sample <- factor(vrr_master_df$sample, levels=unique(sort(vrr_master_df$sample)))
vrr_master_df$Var1 <- factor(vrr_master_df$Var1, levels=c(10500, 9500, 8500, 7500, 6500, 5500, 4500, 3500, 2500, 1500, 500))
  
telo_bar_hist <- ggplot(data=vrr_master_df) +
    geom_bar(mapping=aes(x=sample, y=Freq, fill=Var1), stat='identity') +
    theme_minimal() +
    theme(axis.title.x = element_blank(),
          axis.text = element_text(size=15),
          axis.title.y=element_text(size=20),
          axis.text.x = element_text(size=15, angle=45, hjust=1, vjust=1)) +
    ylab("Proportion of Telomere Reads") +
    guides(fill=guide_legend(title="Telomere Length (BP)")) +
    scale_fill_manual(breaks=c(10500, 9500, 8500, 7500, 6500, 5500, 4500, 3500, 2500, 1500, 500),
                      labels=c("10000+", "9000-9999", "8000-8999", "7000-7999", "6000-6999", "5000-5999", "4000-4999", "3000-3999", "2000-2999", "1000-1999", "0-999"),
                      values=c("#F8766D", "#E68613", "#ABA300", "#0CB702", "#00BE67", "#00BFC4", "#00A9FF", "#8494FF", "#C77CFF", "#FF61CC", "#FF68A1"))
  
ggsave("sample_comparison.vrr_telo_length.pdf", device="pdf", plot=telo_bar_hist, width=12, height=10)

