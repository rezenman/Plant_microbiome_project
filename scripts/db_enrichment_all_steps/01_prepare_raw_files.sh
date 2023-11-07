#!/bin/bash
#steps to perform
#1. create a primer file
#2. change the path to the primer file and the r script below

module load cutadapt/4.2-GCCcore-11.3.0
module load seqkit/2.3.1
module load R/4.0.0

primers="/home/labs/bfreich/shaharr/PacBio/All_raw_fastas/raw_for_pipeline/primers_rc.fasta"
r_script="/home/labs/bfreich/shaharr/PacBio/All_raw_fastas/raw_for_pipeline/reads_per_region.R"

for i in $(find ~+ -maxdepth 1 -name "*fastq.gz" | grep -E "Israel.fastq.gz|France.fastq.gz")
do      
        name=$(echo $i | rev |  cut -d'/' -f 1 | rev | cut -d'.' -f 1)
        mkdir $name
        cd $name

        #2. Run cutadapt with the primer file using the following code
        cutadapt  \
                --no-indels \
                -g file:$primers \
                -o {name}.fastq \
                --action=retain \
                --rc \
                $i &>>"$name"_log.txt

        #3. reverse complement sequences which only contain the reverse primer

        for j in $(ls Rev*)
        do
                seqkit seq $j -g -r -p -o "RC_"$j
        done 

        #4. Add a tag to the header of sequences who contatains only one sequence 
        for k in $(ls Fwd*fastq | grep -Ev "Comp")
        do      
                sed -i 's/ccs/ccs_fwd_only/g' $k
        done
        
        for k in $(ls RC*)
        do      
                sed -i 's/ccs/ccs_reverse_only/g' $k
        done

        #5. concatenate all regions to one file 
        cat Fwd* RC* > All_sequences.fastq

        #6. Trasnform to fasta file using bbmap
        reformat.sh in=All_sequences.fastq out=All_sequences.fasta --nullifybrokenquality fastawrap=0

        #7. replace spaces to underscores
        cat All_sequences.fasta | sed 's/ /_/g' > All_sequences_linearized.fasta
        
        ##to get stats about length "cat All_sequences_linearized.fasta | seqkit seq | seqkit stats"

        #8.filter by length - bigger than 1100 and smaller than 4000
        ##using bbmap, loaded by default
        reformat.sh in=All_sequences_linearized.fasta out="$name"_All_sequences_length_filtered.fasta minlength=1100 maxlength=4000

        #9. print stats
        cat All_sequences_linearized.fasta | seqkit seq | seqkit stats > length_before_trimming.txt
        cat "$name"_All_sequences_length_filtered.fasta | seqkit seq | seqkit stats > length_after_trimming.txt
        zcat $i | seqkit seq | seqkit stats > length_raw_file.txt
        cat *log.txt | grep -E -A 2 "Adapter" > sequences_per_region.txt
        sed -E 's/=== | ===|--//g' sequences_per_region.txt | grep -Ev "^$" | sed 'N;s/\n/; /' | sed -E "s/ [t|T]rimmed| 5' trimmed: [0-9]+ times;| 3'| 5'| times|:|;//g" > mod_report_file.txt
        Rscript $r_script
        cd ..
done
