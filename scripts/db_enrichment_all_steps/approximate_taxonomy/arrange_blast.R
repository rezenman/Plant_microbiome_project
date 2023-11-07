setwd("~/PacBio/All_raw_fastas/db_enrichment/approximate_taxonomy/")
rm(list = ls())


col_names = c("query", "target", "percent_id", "align_length", "num_mismatch", "num_gaps", "qstart", "qend", "tstart", "tend", "e_value", "bit_score")
df = read.table("All_new_seeds_after_splitted_blast_sorted.oufmt6", header = F, col.names = col_names)
head(df)
dim(df)

df_max_bit = df %>% dplyr::group_by(query) %>% filter(bit_score == max(bit_score)) %>% 
  filter(percent_id == max(percent_id))
head(df_max_bit)
dim(df_max_bit)


number_of_matches = as.data.frame(table(df_max_bit$query)) %>% arrange(desc(Freq)) %>% dplyr::rename("query" = "Var1")
distinct_matches = df_max_bit %>% distinct(query, .keep_all = T)
head(distinct_matches)


head_to_tax = read.csv("~/New_UMI_and_Swift_after_smurf/Arrange_taxonomy/SILVA_138.1_SSURef_NR99_tax_silva_trunc_HEADER.csv")
head(head_to_tax)
tail(head_to_tax)

colnames(head_to_tax)[1] = "query"

head_to_tax$target = sapply(strsplit(head_to_tax$taxon, " "), `[`, 1)


head(distinct_matches)
head(head_to_tax)

merged = merge(distinct_matches, head_to_tax  %>% select(-query), by = "target")
dim(merged)
dim(distinct_matches)
head(merged)
merged$blast_par = paste0(merged$bit_score, "_", merged$percent_id, "_", merged$align_length)
head(merged)
head(number_of_matches)
merged = merge(merged, number_of_matches, by = "query") 
head(merged)

write.csv(x = merged %>% select(query, taxon, blast_par, Freq) %>% dplyr::rename("num_headers" = "Freq"), file = "new_seeds_blast_matches.csv")
