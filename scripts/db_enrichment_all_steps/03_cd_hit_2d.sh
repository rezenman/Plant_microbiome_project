#!/bin/bash

module load cd-hit/4.8.1
silva_db="/home/labs/bfreich/shaharr/PacBio/All_raw_fastas/db_enrichment/Silva_with_LNA.fasta"

bsub -q gsla-cpu -R "rusage[mem=2000]" -R "span[hosts=1]" -n 40 -J "cd_hit2d_97"$sample_name -o $sample_name"cd_hit2d_97.out" -e $sample_name"cd_hit2d_97.err" \
cd-hit-est-2d \
    -i $silva_db \
    -i2 all_clustered.99.fasta \
    -o cd_hit_2d_clustered_2.97.fasta -c 0.97 -n 10 -g 1 -d 0 -M 0 -r 0 -T 0
