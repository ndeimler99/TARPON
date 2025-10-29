library(ggplot2)
library(stringr)
library(dplyr)


figure_directory <- "/media/data_01/ndeimler/paternity_test/PATERNITY_ASSIGNMENT/"
output_stats <- "/media/data_01/ndeimler/paternity_test/PATERNITY_ASSIGNMENT/offspring_paternity_stats.txt"

final_df <- read.table(output_stats, header=TRUE, sep="\t")

##### Inheritance Patterns ####

final_df$cluster <- factor(final_df$cluster, level=final_df %>% distinct(cluster, paternity) %>% arrange(paternity, cluster) %>% pull(cluster))

order_df <- final_df %>% group_by(cluster, paternity) %>% summarise(mean_val=mean(telo_length), .groups="drop") %>%
	arrange(paternity, mean_val) %>% mutate(cluster=factor(cluster, levels=cluster))

final_df$cluster <- factor(final_df$cluster, levels=levels(order_df$cluster)) 

plt <- ggplot(data=final_df) +
  geom_boxplot(mapping=aes(x=cluster, y=telo_length, fill=paternity)) + 
  theme_minimal() + ylab("Telomere Length (bp)") +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=18),
        legend.text =element_text(size=14),
        axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=10))

ggsave(paste(figure_directory, "/telomere_length_segregated.png", sep=""), width=16, height=10, plot=plt)


order_df <- final_df %>% group_by(cluster, paternity) %>% 
  summarise(mean_val=mean(telo_length), .groups="drop") %>%
  arrange(mean_val) %>%
  mutate(cluster=factor(cluster, levels=cluster))

final_df$cluster <- factor(final_df$cluster, levels=levels(order_df$cluster))

plt <- ggplot(data=final_df) +
  geom_boxplot(mapping=aes(x=cluster, y=telo_length, fill=paternity)) + 
  theme_minimal() + ylab("Telomere Length (bp)") +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=18),
        legend.text =element_text(size=14),
        axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=10))
plt
ggsave(paste(figure_directory, "/telomere_length_by_cluster.png", sep=""), width=16, height=10, plot=plt)
