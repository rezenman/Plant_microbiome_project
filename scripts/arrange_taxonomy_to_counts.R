#load packages
library(ShortRead)
library(dplyr)

##set the working directory to the root directory where all smurf groups are found in the sub-directories
setwd("/home/labs/bfreich/shaharr/New_UMI_and_Swift_after_smurf/NEW_UMI/Results_New_UMI_samples_16102023_based_UMI_BASED_database_augmented_Oct2023_RL_135_REGIONS_1_2_3_4_5_6/")

##read the csv file (produced by Noam) which match the nes headers in our database(running numbers from 1....) to silva original headers
head_to_tax = read.csv("/home/labs/bfreich/shaharr/New_UMI_and_Swift_after_smurf/Arrange_taxonomy/SILVA_138.1_SSURef_NR99_tax_silva_trunc_HEADER.csv")
head(head_to_tax)
colnames(head_to_tax)[1] = "id" #change the first column name to id for merging in following steps
head_to_tax$blast_par = NA 
head_to_tax$blast_matches = NA #create two new columns with NA values (those will be relevant only to new headers) 
View(head_to_tax)


##read the csv containing the approximated taxonomy for new headers
new_headers_blast_alignments = read.csv("/home/labs/bfreich/shaharr/PacBio/All_raw_fastas/db_enrichment/approximate_taxonomy/new_seeds_blast_matches.csv") %>% 
  select(-X) %>% dplyr::rename("blast_matches" = "num_headers", "id" = "query")
View(new_headers_blast_alignments)


#binding the two taxonomy tables to create one contatinig original taxonomy of headers from silva and approximated taxonomy for new headers
head_to_tax = rbind(head_to_tax, new_headers_blast_alignments)
head(head_to_tax)


##reading the counts data, this is a modified version of the smurf output (delete the first line and coulumns "domain":"species")
counts = read.table("/home/labs/bfreich/shaharr/New_UMI_and_Swift_after_smurf/NEW_UMI/Results_New_UMI_samples_16102023_based_UMI_BASED_database_augmented_Oct2023_RL_135_REGIONS_1_2_3_4_5_6/groups_arranged.txt",
                    header = T, check.names = F)
colnames(counts)
dim(counts)
head(counts)


##listing all groups produced by smurf
all_groups = list.files(path = "./Groups/", pattern = "^sample.*fasta$", recursive = T, full.names = T)
head(all_groups)
all_groups

##seperating the hash from the groups full name
only_hash = sapply(strsplit(all_groups, "group[0-9].*_|\\.fasta$"), `[`, 2)

##keeping only unique hashes (the same group with the same hash can appear in two samples but contains same sequences)
only_unique_hash = unique(sapply(strsplit(all_groups, "group[0-9].*_|\\.fasta$"), `[`, 2))


#finding only one of the matching groups per hash
group_per_hash = NULL
for(i in 1:length(only_unique_hash)){

  print(paste0("analyzing sample: ", i, " -------- ", length(only_unique_hash)))
  pattern = only_unique_hash[i]
  match = grep(pattern = pattern, all_groups, value = T)[1]
  group_per_hash[i] = match
  
}


##for each hash finding its first id number from the fasta file and the number of headers in a afasta file, and creating a data frame 
ids = NULL
num_headers = NULL
hash = NULL

for(j in 1:length(only_unique_hash)){
  print(paste0("analyzing sample: ", j, " -------- ", length(only_unique_hash)))

  hash_to_look = only_unique_hash[j]
  fasta_to_read = grep(hash_to_look, group_per_hash, value = T)
  fasta = readFasta(fasta_to_read)
  
  id = as.data.frame(fasta@id[1])[1,1]
  id = strsplit(id, "\\.")[[1]][1]
  len = dim(as.data.frame(fasta@sread))[1]
  
  ids[j] = id
  hash[j] = hash_to_look
  num_headers[j] = len
  
}

hash_to_id = data.frame(hash = hash, id = as.numeric(ids), num_headers = num_headers)
head(hash_to_id)
head(head_to_tax)


#merging has and id with taxon, arrangig and assigining the string "No_Taxonomy" to headers which had no match
id_to_taxon_to_hash = merge(hash_to_id, head_to_tax, by = "id", all.x= T)
id_to_taxon_to_hash$id = as.numeric(id_to_taxon_to_hash$id)
id_to_taxon_to_hash[is.na(id_to_taxon_to_hash$taxon), 4] = "No_Taxonomy"
colnames(id_to_taxon_to_hash)[2] = "groups"
head(id_to_taxon_to_hash)

##merging the count table of each sample with the taxonomy,id, hash table
counts_with_taxonomy = merge(counts, id_to_taxon_to_hash, by = "groups")

## Arranging the table to containg all metadata in the first coulumns and only then each sample and its counts
counts_with_taxonomy = cbind(counts_with_taxonomy %>% select(groups, taxon, id, num_headers, blast_par, blast_matches), 
                             counts_with_taxonomy %>% select(-c(groups, taxon, id, num_headers, blast_par, blast_matches)))
head(counts_with_taxonomy)



##optinal write the full table to csv file
write.csv(x = counts_with_taxonomy, file = "counts_with_taxonomy_fixed.csv")



# seperate each sample to its own csv file----------------------------------------------------

dir.create("seperate_samples")  #this is only to create a seperate folder to contain all seperate samples csv files

#for each sample writing a seperate csv file
for(i in 5:dim(counts_with_taxonomy)[2]){
  sample_name = colnames(counts_with_taxonomy)[i]
  
  print(paste0("Analyzing sample: ", sample_name))
  
  df_loop = counts_with_taxonomy[,c(1:6, i)] %>% filter(.data[[sample_name]] > 0 )
  save_name = paste0("seperate_samples/", sample_name, ".csv")
  write.csv(x = df_loop, file = save_name)
}

