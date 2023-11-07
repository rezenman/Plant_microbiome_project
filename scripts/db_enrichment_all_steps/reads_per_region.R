library(dplyr)
library(ggplot2)

## create a table for reads per regions
sample_name = strsplit(getwd(), "/")[[1]][9]

regions = read.table("mod_report_file.txt", sep = " ",  fill = NA) %>% select(-c("V1", "V3", "V5", "V6", "V7", "V8", "V10"))
colnames(regions) = c("primers", "primers_sequences", "num_reads", "num_reads_comp")
regions$num_reads_non_comp = regions$num_reads - regions$num_reads_comp
regions$sample_name = sample_name
  
regions_plot = ggplot(regions, aes(x = primers, y = num_reads)) + geom_col() + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90)) + 
  labs(x = "Primer combinations", y = "Number of sequences")

write.csv(file = "reads_per_region.csv", x = regions)
ggsave(filename = "regions_plot.jpg", regions_plot, units = "in", height = 12, width = 12)

#create a table for total reads before and after
total_reads_list = list.files(pattern = "length.*txt", full.names = T)
names(total_reads_list) = c("After_trimming", "After_demultiplex", "Raw_reads")
df_list = lapply(total_reads_list, read.table, header = T)
df = do.call("rbind", df_list)
df$file = rownames(df)
rownames(df) = 1:3
df$num_seqs = as.numeric(gsub(df$num_seqs, pattern = ",", replacement = ""))
df$file = factor(df$file, levels = c("Raw_reads", "After_demultiplex", "After_trimming"))
df$sample_name = sample_name

total_reads = ggplot(df, aes(x = file, y = num_seqs, label = num_seqs)) + geom_col() + 
  theme_bw() + 
  geom_label() + 
  labs(x = "Step", y = "Number of sequences")

write.csv(file = "total_reads.csv", x = df)
ggsave(filename = "total_reads.jpg", total_reads, units = "in", height = 12, width = 12)