#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)
library(stringr)

COLORS=c('#0173b2', '#de8f05', '#029e73', '#d55e00', '#cc78bc', '#ca9161', '#fbafe4', '#949494', '#ece133', '#56b4e9')

args = commandArgs(trailingOnly=TRUE)

stats <- read.table(args[1], header=TRUE)
repeat_distribution <- read.table(args[2], header=FALSE)
repeat_seq <- args[3]
mutant <- args[4]

colnames(repeat_distribution) <- c("repeat", "pos", "replacement", "freq")

if (mutant != "false") {
  
  variant_histogram <- ggplot(data=stats) +
    geom_histogram(mapping=aes(one_nucl_variant_composition, fill="one_nucl_comp"), binwidth=1, alpha=0.5) +
    geom_histogram(mapping=aes(wt_composition, fill="wt_comp"), binwidth=1, alpha=0.5) +
    geom_histogram(mapping=aes(mutant_composition, fill="mt_comp"), binwidth=1, alpha=0.5) +
    theme_minimal() +
    theme(axis.title=element_text(size=20),
          axis.text=element_text(size=15),
          legend.title=element_blank(),
          legend.text=element_text(size=15)) +
    xlab("Telomere Composition (%)") + ylab("Number of Telomeric Reads") +
    scale_fill_manual(breaks=c("wt_comp", "mt_comp", "one_nucl_comp"), values=c("#D81B60", "#FFC107", "#1E88E5"), labels=c(repeat_seq, mutant, "One Nucleotide Substitutions"))
  
  ggsave("variant_repeat_composition_histogram.pdf", plot=variant_histogram, width=12, height=10, device="pdf")
  
  variant_wt_scatter <- ggplot(data=stats) +
    geom_point(mapping=aes(x=wt_composition, y=one_nucl_variant_composition, color="one_nucl")) +
    geom_point(mapping=aes(x=wt_composition, y=mutant_composition, color="mutant")) +
    theme_minimal() +
    xlab(paste(repeat_seq, "(%)", sep="")) + ylab("Percentage of Repeat") +
    theme(axis.text=element_text(size=15),
          axis.title=element_text(size=20)) +
    scale_color_manual(breaks=c("one_nucl", "mutant"), values=c("#D81B60", "#FFC107"), labels=c("One Nucleotide Substitutions", mutant))
  
  ggsave("variant_repeat_by_wt_composition.pdf", plot=variant_wt_scatter, width=12, height=10, device="pdf")
  
  variant_wt_scatter <- ggplot(data=stats) +
    geom_point(mapping=aes(x=wt_composition, y=vrr_telo_length, color="wt_comp")) +
    geom_point(mapping=aes(x=one_nucl_variant_composition, y=vrr_telo_length, color="one_nucl_comp")) +
    geom_point(mapping=aes(x=mutant_composition, y=vrr_telo_length, color="mt_comp")) +
    theme_minimal() +
    scale_color_manual(breaks=c("wt_comp", "mt_comp", "one_nucl_comp"), values=c("#D81B60", "#FFC107", "#1E88E5"), labels=c(repeat_seq, mutant, "One Nucleotide Substitutions")) +
    theme(axis.title=element_text(size=20),
          axis.text=element_text(size=15),
          legend.title=element_blank(),
          legend.text=element_text(size=15)) +
    xlab("Telomere Composition (%)") + ylab("VRR Telomere Length (bp)")
  
  ggsave("variant_repeat_composition_by_telo_length.pdf", plot=variant_wt_scatter, width=6, height=10, device="pdf")
  
  
  variant_wt_boxplot <-  ggplot(data=stats) +
    geom_boxplot(mapping=aes(x=repeat_seq, y=wt_composition, color="wt_comp")) +
    geom_boxplot(mapping=aes(x=mutant, y=mutant_composition, color="mt_comp")) +
    geom_boxplot(mapping=aes(x="One Nucleotide Substitutions", y=one_nucl_variant_composition, color="one_nucl_comp")) +
    theme_minimal() +
    scale_color_manual(breaks=c("wt_comp", "mt_comp", "one_nucl_comp"), values=c("#D81B60", "#FFC107", "#1E88E5"), labels=c(repeat_seq, mutant, "One Nucleotide Substitutions")) +
    theme(axis.title=element_text(size=20),
          axis.text=element_text(size=15),
          legend.title=element_blank(),
          legend.text=element_text(size=15),
          legend.position = "none",
          axis.title.x=element_blank()) +
    ylab("Telomere Composition (%)")
  
  ggsave("variant_repeat_composition_boxplot.pdf", plot=variant_wt_boxplot, width=6, height=10, device="pdf")
  
  
  stats <- stats %>%
    arrange(vrr_telo_length, desc(one_nucl_variant_composition))
  
  telo_lengths <- ggplot(data=stats) +
    geom_bar(mapping=aes(x=seq(1:length(vrr_telo_length)), y=vrr_telo_length), width=0.7, stat="identity") + 
    xlab("Telomeric Read Index") + ylab("VRR Telomere Length (BP)") +
    coord_flip() + theme_minimal() +
    theme(axis.title=element_text(size=20),
          axis.text=element_text(size=20)) 
  
  ggsave("vrr_telo_lengths_barplot.pdf", plot=telo_lengths, width=12, height=10, device="pdf")
  
  new_stats <- stats[,c("read_id", "vrr_telo_length", "wt_composition", "mutant_composition", "one_nucl_variant_composition")]
  new_stats$wt <- (new_stats$wt_composition/100) * new_stats$vrr_telo_length
  new_stats$variant <- (new_stats$one_nucl_variant_composition/100) * new_stats$vrr_telo_length 
  new_stats$mutant <- (new_stats$mutant_composition/100) * new_stats$vrr_telo_length
  new_stats$neither <- new_stats$vrr_telo_length * ((100-(new_stats$wt_composition+new_stats$one_nucl_variant_composition+new_stats$mutant_composition))/100)
  new_stats <- new_stats[,c("read_id", "wt", "mutant", "variant", "neither")]
  new_stats <- reshape2::melt(new_stats)
  new_stats$read_id <- factor(new_stats$read_id, levels=stats$read_id)
  
  telo_lengths_colored <- ggplot(data=new_stats) +
    geom_bar(mapping=aes(x=read_id, y=value, fill=variable), stat="identity") +
    coord_flip() +
    theme_minimal() +
    scale_fill_manual(breaks=c("wt", "mutant", "variant", "neither"), values=c("#D81B60","#FFC107", "#1E88E5", "#004D40"), labels=c(repeat_seq, mutant, "One Nucleotide Substitutions", "Other")) +
    theme(axis.text.y = element_blank(),
          axis.text=element_text(size=15),
          axis.title=element_text(size=20),
          legend.title=element_blank(),
          legend.text=element_text(size=15),
          panel.grid.major = element_blank()) +
    ylab("VRR Telomere Length (bp)") + xlab("Telomeric Read Index")
  
  ggsave("vrr_telo_lengths_barplot_by_composition.pdf", plot=telo_lengths_colored, width=12, height=10, device="pdf")
  
  if (sum(!unlist(str_split(repeat_seq, "")) == unlist(str_split(mutant, ""))) == 1){
    repeat_distribution$pos <- factor(repeat_distribution$pos, levels=c(sort(unique(repeat_distribution$pos), decreasing=TRUE)))
    repeat_distribution_plot <- ggplot(data=repeat_distribution) +
      geom_bar(mapping=aes(x=pos, y=freq, fill=replacement), stat='identity') +
      coord_flip() +
      ylab("Percentage of Repeats") + theme_minimal() +
      theme(axis.text=element_text(size=15),
            axis.title=element_text(size=20),
            axis.title.y = element_blank(),
            legend.title=element_blank(),
            legend.text=element_text(size=15)) +
      scale_fill_manual(values=c("#D81B60", "#1E88E5", "#FFC107", "#004D40")) +
      scale_x_discrete(breaks=seq(0,(nchar(repeat_seq)-1)), labels=unlist(str_split(repeat_seq, "")))
    
    ggsave("telomere_composition_one_nucleotide_variants.pdf", width=12, height=10, device="pdf", plot=repeat_distribution_plot)
  } else {
    repeat_distribution$pos <- factor(repeat_distribution$pos, levels=c(sort(unique(repeat_distribution$pos), decreasing=TRUE)))
    repeat_distribution_plot <- ggplot(data=repeat_distribution) +
      geom_bar(mapping=aes(x=pos, y=freq, fill=replacement), stat='identity') +
      coord_flip() +
      ylab("Percentage of Repeats") + theme_minimal() +
      theme(axis.text=element_text(size=15),
            axis.title=element_text(size=20),
            axis.title.y = element_blank(),
            legend.title=element_blank(),
            legend.text=element_text(size=15)) +
      scale_fill_manual(breaks=c("A", "T", "C", "G", "mutant"), values=c("#D81B60", "#1E88E5", "#FFC107", "#004D40", "red")) +
      scale_x_discrete(breaks=c(seq(0,(nchar(repeat_seq)-1)),"mutant"), labels=c(unlist(str_split(repeat_seq, "")),mutant))
    
    ggsave("telomere_composition_one_nucleotide_variants.pdf", width=12, height=10, device="pdf", plot=repeat_distribution_plot)
  }
  
} else {
  variant_histogram <- ggplot(data=stats) +
    geom_histogram(mapping=aes(one_nucl_variant_composition, fill="one_nucl_comp"), binwidth=1, alpha=0.5) +
    geom_histogram(mapping=aes(wt_composition, fill="wt_comp"), binwidth=1, alpha=0.5) +
    theme_minimal() +
    theme(axis.title=element_text(size=20),
          axis.text=element_text(size=15),
          legend.title=element_blank(),
          legend.text=element_text(size=15)) +
    xlab("Telomere Composition (%)") + ylab("Number of Telomeric Reads") +
    scale_fill_manual(breaks=c("wt_comp", "one_nucl_comp"), values=c("#D81B60", "#1E88E5"), labels=c(repeat_seq, "One Nucleotide Substitutions"))
  
  ggsave("variant_repeat_composition_histogram.pdf", plot=variant_histogram, width=12, height=10, device="pdf")
  
  variant_wt_scatter <- ggplot(data=stats) +
    geom_point(mapping=aes(x=wt_composition, y=one_nucl_variant_composition)) +
    theme_minimal() +
    xlab(paste(repeat_seq, "(%)", sep="")) + ylab("One Nucleotide Substitutions (%)") +
    theme(axis.text=element_text(size=15),
          axis.title=element_text(size=20))
  
  ggsave("variant_repeat_by_wt_composition.pdf", plot=variant_wt_scatter, width=12, height=10, device="pdf")
  
  variant_wt_scatter <- ggplot(data=stats) +
    geom_point(mapping=aes(x=wt_composition, y=vrr_telo_length, color="wt_comp")) +
    geom_point(mapping=aes(x=one_nucl_variant_composition, y=vrr_telo_length, color="one_nucl_comp")) +
    theme_minimal() +
    scale_color_manual(breaks=c("wt_comp", "one_nucl_comp"), values=c("#D81B60", "#1E88E5"), labels=c(repeat_seq, "One Nucleotide Substitutions")) +
    theme(axis.title=element_text(size=20),
          axis.text=element_text(size=15),
          legend.title=element_blank(),
          legend.text=element_text(size=15)) +
    xlab("Telomere Composition (%)") + ylab("VRR Telomere Length (bp)")
  
  ggsave("variant_repeat_composition_by_telo_length.pdf", plot=variant_wt_scatter, width=6, height=10, device="pdf")
  
  variant_wt_boxplot <-  ggplot(data=stats) +
    geom_boxplot(mapping=aes(x=repeat_seq, y=wt_composition, color="wt_comp")) +
    geom_boxplot(mapping=aes(x="One Nucleotide Substitutions", y=one_nucl_variant_composition, color="one_nucl_comp")) +
    theme_minimal() +
    scale_color_manual(breaks=c("wt_comp", "one_nucl_comp"), values=c("#D81B60", "#1E88E5"), labels=c(repeat_seq, "One Nucleotide Substitutions")) +
    theme(axis.title=element_text(size=20),
          axis.text=element_text(size=15),
          legend.title=element_blank(),
          legend.text=element_text(size=15),
          legend.position = "none",
          axis.title.x=element_blank()) +
    ylab("Telomere Composition (%)")
  
  ggsave("variant_repeat_composition_boxplot.pdf", plot=variant_wt_boxplot, width=6, height=10, device="pdf")
  
  
  stats <- stats %>%
              arrange(vrr_telo_length, desc(one_nucl_variant_composition))

  telo_lengths <- ggplot(data=stats) +
    geom_bar(mapping=aes(x=seq(1:length(vrr_telo_length)), y=vrr_telo_length), width=0.7, stat="identity") + 
    xlab("Telomeric Read Index") + ylab("VRR Telomere Length (BP)") +
    coord_flip() + theme_minimal() +
    theme(axis.title=element_text(size=20),
          axis.text=element_text(size=20)) 
  
  ggsave("vrr_telo_lengths_barplot.pdf", plot=telo_lengths, width=12, height=10, device="pdf")
 
  new_stats <- stats[,c("read_id", "vrr_telo_length", "wt_composition", "one_nucl_variant_composition")]
  new_stats$wt <- (new_stats$wt_composition/100) * new_stats$vrr_telo_length
  new_stats$variant <- (new_stats$one_nucl_variant_composition/100) * new_stats$vrr_telo_length 
  new_stats$neither <- new_stats$vrr_telo_length * ((100-(new_stats$wt_composition+new_stats$one_nucl_variant_composition))/100)
  new_stats <- new_stats[,c("read_id", "wt", "variant", "neither")]
  new_stats <- reshape2::melt(new_stats)
  new_stats$read_id <- factor(new_stats$read_id, levels=stats$read_id)
  
  telo_lengths_colored <- ggplot(data=new_stats) +
    geom_bar(mapping=aes(x=read_id, y=value, fill=variable), stat="identity") +
    coord_flip() +
    theme_minimal() +
    scale_fill_manual(breaks=c("wt", "variant", "neither"), values=c("#D81B60", "#1E88E5", "#004D40"), labels=c(repeat_seq, "One Nucleotide Substitutions", "Other")) +
    theme(axis.text.y = element_blank(),
          axis.text=element_text(size=15),
          axis.title=element_text(size=20),
          legend.title=element_blank(),
          legend.text=element_text(size=15),
          panel.grid.major = element_blank()) +
    ylab("VRR Telomere Length (bp)") + xlab("Telomeric Read Index")
  
  ggsave("vrr_telo_lengths_barplot_by_composition.pdf", plot=telo_lengths_colored, width=12, height=10, device="pdf")
  
  
  repeat_distribution$pos <- factor(repeat_distribution$pos, levels=c(sort(unique(repeat_distribution$pos), decreasing=TRUE)))
  repeat_distribution_plot <- ggplot(data=repeat_distribution) +
    geom_bar(mapping=aes(x=pos, y=freq, fill=replacement), stat='identity') +
    coord_flip() +
    ylab("Percentage of Repeats") + theme_minimal() +
    theme(axis.text=element_text(size=15),
          axis.title=element_text(size=20),
          axis.title.y = element_blank(),
          legend.title=element_blank(),
          legend.text=element_text(size=15)) +
  scale_fill_manual(values=c("#D81B60", "#1E88E5", "#FFC107", "#004D40")) +
  scale_x_discrete(breaks=seq(0,(nchar(repeat_seq)-1)), labels=unlist(str_split(repeat_seq, "")))
  
  ggsave("telomere_composition_one_nucleotide_variants.pdf", width=12, height=10, device="pdf", plot=repeat_distribution_plot)
}

