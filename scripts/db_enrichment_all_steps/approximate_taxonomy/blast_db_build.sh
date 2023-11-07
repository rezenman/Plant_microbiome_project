#!/bin/bash
source activate blast

db="/home/labs/bfreich/shaharr/PacBio/All_raw_fastas/db_enrichment/approximate_taxonomy/Silva_NR99_shortened_headers.fasta"

bsub -q gsla-cpu -n 1 -R "rusage[mem=32000]" -J silva_db_build -o silva_db_build.out -e silva_db_build.err \
makeblastdb -in $db -parse_seqids -title silva -dbtype nucl 