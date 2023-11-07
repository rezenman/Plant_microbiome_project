#!/bin/bash
module load seqkit/2.3.1
module load seqtk/1.3-GCC-8.3.0


##Create an original header to new header tab dilimited file
grep -E "^>" cd_hit_2d_clustered_2.99.fasta | sed 's/>//g' | awk 'BEGIN{ FS = OFS = "\t" } { print $0, (NR+1999999) }' > original_headers_2.txt 

##create a new fasta file with new headers
seqkit replace -p '(m.*)$' -r '{kv}' -k original_headers_2.txt cd_hit_2d_clustered_2.99.fasta > Sequences_to_add_non_linear.fasta

##linearize sequences 
seqtk seq -l 0 Sequences_to_add_non_linear.fasta > Seqs_to_add_30082023.fasta
