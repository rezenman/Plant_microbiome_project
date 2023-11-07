#!/bin/bash
module load cd-hit/4.8.1

bsub -q gsla-cpu -R "rusage[mem=5000]" -R "span[hosts=1]" -n 10 -J cd_hit_"$sample_name" -o $sample_name"_cd_hit99.out" -e $sample_name"_cd_hit99.err" \
cd-hit-est -i All_ready_for_cd_hit.fasta -o all_clustered.99.fasta -c 0.99 -n 10 -d 0 -M 0 -r 0 -g 1
