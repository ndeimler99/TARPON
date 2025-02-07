#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)
library(stringr)

COLORS=c('#0173b2', '#de8f05', '#029e73', '#d55e00', '#cc78bc', '#ca9161', '#fbafe4', '#949494', '#ece133', '#56b4e9')

args = commandArgs(trailingOnly=TRUE)

wt_processivity <- read.table(args[1], header=FALSE)
mt_processivity <- read.table(args[2], header=FALSE)
repeat_seq <- args[3]
mutant <- args[4]

# wt_processivity <- read.table("/media/data_01/ndeimler/WORKFLOWS/TArPON/work/af/cffac0af2dfcb373d66dc0ea90422b/sample.wt.processivity_stats.txt", header=FALSE)
# mt_processivity <- read.table("/media/data_01/ndeimler/WORKFLOWS/TArPON/work/af/cffac0af2dfcb373d66dc0ea90422b/sample.mt.processivity_stats.txt", header=FALSE)

colnames(wt_processivity) <- c("repeat_count", "occurences")
colnames(mt_processivity) <- c("repeat_count", "occurences")

wt_processivity$percentage_of_occurences <- wt_processivity$occurences / sum(wt_processivity$occurences)
wt_processivity$percentage_of_repeats <- (wt_processivity$occurences * wt_processivity$repeat_count) / (sum(wt_processivity$occurences * wt_processivity$repeat_count))

mt_processivity$percentage_of_occurences <- mt_processivity$occurences / sum(mt_processivity$occurences)
mt_processivity$percentage_of_repeats <- (mt_processivity$occurences * mt_processivity$repeat_count) / (sum(mt_processivity$occurences * mt_processivity$repeat_count))

wt_processivity$repeat_type <- "wild_type"
mt_processivity$repeat_type <- "mutant"

processivity_df <- rbind(wt_processivity, mt_processivity)

ggplot(data=processivity_df) +
  geom_bar(mapping=aes(x=repeat_count, y=percentage_of_occurences, fill=repeat_type), stat="identity", position=position_dodge2()) +
  theme_minimal() +
  xlab("Number of Consecutive Repeats") + ylab("Number of Occurences") +
  theme(axis.title=element_text(size=20),
        axis.text=element_text(size=15),
        legend.title=element_blank(),
        legend.text=element_text(size=15)) +
  scale_fill_manual(breaks=c("wild_type", "mutant"), values=c("#D81B60", "#FFC107"), labels=c(repeat_seq, mutant)) + 
  xlim(0,100)

ggplot(data=processivity_df) +
  geom_bar(mapping=aes(x=repeat_count, y=percentage_of_repeats, fill=repeat_type), stat="identity", position=position_dodge2()) +
  theme_minimal() +
  xlab("Number of Consecutive Repeats") + ylab("Number of Occurences") +
  theme(axis.title=element_text(size=20),
        axis.text=element_text(size=15),
        legend.title=element_blank(),
        legend.text=element_text(size=15)) +
  scale_fill_manual(breaks=c("wild_type", "mutant"), values=c("#D81B60", "#FFC107"), labels=c(repeat_seq, mutant)) + 
  xlim(0,100)

ggplot(data=processivity_df) +
  geom_point(mapping=aes(x=repeat_count, y=log(percentage_of_occurences), color=repeat_type), size=4) +
  theme_minimal() +
  xlab("Number of Consecutive Repeats") + ylab("Number of Occurences") +
  theme(axis.title=element_text(size=20),
        axis.text=element_text(size=15),
        legend.title=element_blank(),
        legend.text=element_text(size=15)) +
  scale_fill_manual(breaks=c("wild_type", "mutant"), values=c("#D81B60", "#FFC107"), labels=c(repeat_seq, mutant)) + 
  xlim(0,100)



