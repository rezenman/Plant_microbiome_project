#!/bin/bash
source activate blast

bsub -q gsla-cpu -n 20 -R "rusage[mem=2000]" -R "span[hosts=1]" -J blastn -o blastn.out -e blastn.err \
blastn -query "/home/labs/bfreich/shaharr/PacBio/All_raw_fastas/db_enrichment/Seqs_to_add_30082023.fasta" \
 -db Silva_NR99_shortened_headers.fasta -outfmt 6 -evalue 1e-5 -out New_seeds_smurf.outfmt6  -num_threads 20