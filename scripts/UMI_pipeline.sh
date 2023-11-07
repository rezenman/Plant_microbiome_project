#!/bin/bash
module load cutadapt/4.2-GCCcore-11.3.0
module load cd-hit/4.8.1
module load BBMap/38.90-GCC-10.2.0
module load seqtk/1.2
module load seqkit/2.3.1
module load R/3.5.1

#run as following:
# the script will create all files in the current directory
# bsub -q gsla-cpu -n 20 -R "rusage[mem=5000]" -R "span[hosts=1]" -J new_pip -o new_pip.out -e new_pip.err \
#  ./UMI_pipeline.sh

# for_primers="/home/labs/bfreich/shaharr/new_microbiome_pipeline/forward.fasta"
# rev_primers="/home/labs/bfreich/shaharr/new_microbiome_pipeline/reverse.fasta"
# read_1="/home/labs/bfreich/shaharr/new_microbiome_pipeline/380S_7RXS_S16_R1_001.fastq.gz"
# read_2="/home/labs/bfreich/shaharr/new_microbiome_pipeline/380S_7RXS_S16_R2_001.fastq.gz"
# sample_name="Test"

for_primers=$1
rev_primers=$2
read_1=$3
read_2=$4
sample_name=$5
rscript="/home/labs/bfreich/shaharr/new_microbiome_pipeline/aux_dada_script.R"

##seperating to regions with proper V pairing (V1V2 forward can only be with V1V2 reverse etc.)
##This step also trasnform the adapter, umi and stepper to lower case letters, log with number of reads per region is "cutadapt_log.txt"
cutadapt \
    --no-indels \
    --action=lowercase \
    -g file:$for_primers \
    -G file:$rev_primers \
    --cores=0 \
    -o {name1}-{name2}.1.fastq.gz -p {name1}-{name2}.2.fastq.gz \
     $read_1 $read_2 &>> cutadapt_log.txt

##for every variable region, creating a seperate fasta file with only headers, and the primer,stepper and a umi

#find all combinations of regions found, keep only those with more than 10000 reads
var_regions=$(find . -name "V*-V*.1.fastq.gz" -exec zgrep -EHc "^@" {} \; | sed 's/:/ /g' | sort -k 2n | awk '$2 >= 10000 {print $1}')

echo "region" "num_read1" "num_read2" "num_seeds" "num_read1_final_after_length_filter" "num_read2_final_after_length_filter" > clustering_log.txt

for i in $var_regions
do  
    #extracting region name, finding the name files read 1 and read 2
    region=$(echo $i | sed 's/.1.fastq.gz//g' | cut -d'/' -f 2)
    r1=$(ls $region*.1.fastq.gz)
    r2=$(ls $region*.2.fastq.gz)

    echo "Analyzing region:" $region

    #Extracting only lowercase letters from sequences of read 1 and read2, and extracting only headers to another file
    zgrep -Eo "^[a-z]+" $r1 > $region'_sequences1.txt'
    zgrep -Eo "^[a-z]+" $r2 > $region'_sequences2.txt'
    zgrep -E "^@" $r1 | sed  "s/^@/>/g" > $region'_headers.txt'

    #Concatenating lowercase sequences from both reads (stepper umis and adapters from both sides)
    paste -d '' $region'_sequences1.txt' $region'_sequences2.txt' > $region'_MergedSequences.txt'
    #Concatenating the headers with the concatenated sequences to create a fasta file to continue to work on
    paste -d '\n' $region'_headers.txt' $region'_MergedSequences.txt' > $region'_adapter_umi_stepper.fasta'

    #printing an histogram of umi,stepper,adapter lengths (most suppose to be between 46-56 bp length)
    readlength.sh in=$region'_adapter_umi_stepper.fasta' out=$region'_histogram.txt'

    #Clustering the umi,stepper,adapter sequences to remove redundant sequences and keep only one seed per group
    dedupe.sh in=$region'_adapter_umi_stepper.fasta' out=$region'_adapter_umi_stepper_97.fasta' outd=$region'_duplicates.fasta'

    #Extracting the seed headers, and filtering the original fastq files (following cutadapt) to contain only one sequence per cluster
    grep "^>" $region'_adapter_umi_stepper_97.fasta' | sed 's/^>//g' > $region'_headers_to_keep_after_clustering.txt'
    filterbyname.sh \
        in=$r1 in2=$r2 \
        out=$region'_filtered.1.fastq.gz' out2=$region'_filtered.2.fastq.gz' \
        names=$region'_headers_to_keep_after_clustering.txt' include=t overwrite=true

    #a sanity check to remove duplicated sequences,s since the previous command tens to create duplicates by mistake
    seqkit rmdup -n $region'_filtered.1.fastq.gz' -o $region'_filtered_dedup.1.fastq.gz'
    seqkit rmdup -n $region'_filtered.2.fastq.gz' -o $region'_filtered_dedup.2.fastq.gz'

    #running cutadapt again this time removing the lowercase letters and their quality scores, for each region seperately
    primer_for_2nd_cut=$(echo $region | cut -d'-' -f 1 | grep -Ef - $for_primers -A1 | tail -n 1)
    primer_rev_2nd_cut=$(echo $region | cut -d'-' -f 2 | grep -Ef - $rev_primers -A1 | tail -n 1)

    cutadapt \
        -e 0.15 \
        --pair-adapters \
        --no-indels \
        --minimum-length 100 \
        --action=trim \
        -g $primer_for_2nd_cut \
        -G $primer_rev_2nd_cut \
        --cores=0 \
        -o $region'_for_dada2.1.fastq.gz' -p $region'_for_dada2.2.fastq.gz' \
        $region'_filtered_dedup.1.fastq.gz' $region'_filtered_dedup.2.fastq.gz'


    find . -name "*unknown*" -exec rm {} \;
    
    # counting number of reads before and after and printing to report
    num_read1=$(zgrep -Ec "^@" $r1)
    num_read2=$(zgrep -Ec "^@" $r2)
    num_seeds=$(grep -Ec "^>" $region'_adapter_umi_stepper_97.fasta')
    num_read1_final=$(zgrep -Ec "^@" $region'_for_dada2.1.fastq.gz')
    num_read2_final=$(zgrep -Ec "^@" $region'_for_dada2.2.fastq.gz')

    echo $region $num_read1 $num_read2 $num_seeds $num_read1_final $num_read2_final >> clustering_log.txt
done

#Creating directories for fwd and rev and moving all relevant files there to continue to the dada pipline

mkdir FWD
mkdir REV

find . -maxdepth 1 -name "*_for_dada2.1.fastq.gz" -exec mv {} "./FWD" \;
find . -maxdepth 1 -name "*_for_dada2.2.fastq.gz" -exec mv {} "./REV" \;

Rscript $rscript

##Arranging data output fo smurf
# replacing all spaces to new lines to create fasta file
for i in $(ls Samp*.fasta); do sed 's/ /\n/g' $i > "$i"_mod ; done

# adding back primers for smurf analysis
#resolving ambigous base-pairs: M, R, W, V, H -> A
#                               Y, B, N, S -> C
#                               K, D -> G

for fa in $(ls *_for.fasta_mod)
  do
    fwd=$(echo $fa | sed -E 's/Samp|_for.*//g' | cut -d'-' -f 1)
    seq_for=$(grep -A1 $fwd $for_primers | grep -v $fwd | sed -E 's/M|R|W|V|H/A/g' | sed -E 's/Y|B|N|S/C/g' | sed -E 's/K|D/G/g')
    sed "2~2 s/^/$seq_for/" $fa > $fa"_primers_added"  
  done

for fa in $(ls *_rev.fasta_mod)
  do
    rev=$(echo $fa | sed -E 's/Samp|_for.*//g' | cut -d'-' -f 2)
    seq_rev=$(grep -A1 $rev $rev_primers | grep -v $rev | sed -E 's/M|R|W|V|H/A/g' | sed -E 's/Y|B|N|S/C/g' | sed -E 's/K|D/G/g')
    sed "2~2 s/^/$seq_rev/" $fa > $fa"_primers_added"  
  done

# reformating fasta files as pseudo fastq files with fake quality score of 40
for i in $(ls *primers_added); do reformat.sh in=$i out="$i".fastq qfake=40 ; done

##moving files to folders
mkdir fastq_with_primers
find . -maxdepth 1 -name "*primers_added.fastq" -exec mv {} "./fastq_with_primers" \; 

mkdir fasta_files
find . -maxdepth 1 -name "*.fasta*" -exec mv {} "./fasta_files" \; 

mkdir text_files
find . -maxdepth 1 -name "*.txt" -exec mv {} "./text_files" \; 

mkdir compressed_fastq_files
find . -maxdepth 1 -name "*.fastq.gz" -exec mv {} "./compressed_fastq_files" \; 



## Concatenating all regions to create the final files
cd fastq_with_primers

#locating files that are not perfect pairing of for and rev primers (e.g, V1V2_f and V2V3_r) and changin their names
for i in $(ls *fastq)
do
  reg_1=$(echo $i | grep -Eo 'V[0-9]V[0-9]|ITS[0-9]' | head -n 1)
  reg_2=$(echo $i | grep -Eo 'V[0-9]V[0-9]|ITS[0-9]' | tail -n 1)
  echo $reg_1 $reg_2
  if [ $reg_1 == $reg_2 ]; then echo "Good"; else mv $i 'Diff_regions_'$i && echo "changed name"; fi
done

name_1=$sample_name"_L001_R1_001.fastq.gz"
name_2=$sample_name"_L001_R2_001.fastq.gz"

#concatenating all fastq files with proper primer pairing (e.g, only V1V2_f with V1V2_r, V2V3_f with V2V3_r etc.) and gzipping
cat Samp*_for.fasta*fastq | gzip -c > $name_1 
cat Samp*_rev.fasta*fastq | gzip -c > $name_2 

echo "Finished All"
##final files to transfer to Noam for SMURF analysis called: $sample_name"_L001_R1_001.fastq.gz" and $sample_name"_L001_R2_001.fastq.gz"