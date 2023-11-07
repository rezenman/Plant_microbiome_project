#!/bin/bash

## this script will receive a fasta/fastq file from initial pacbio analysis, will remove redundant sequences and cluster against a chosen database to find sequences
## which does not match the database

##to use: ./enorich_db.sh -I "fasta/fastq input file" -D "database to compare to in fasta format" -N "sample name" -c "OPTIONAL - clustering identity to use with cd-hit"
## the only package needed is cd-hit, here I use the already installed version on the Weizmann institute wexac cluster:
module load cd-hit/4.8.1

clustering_id=0.99
## enlisting all arguments required to run the script
while getopts 'c:hI:D:N:' OPTION; do
  case "$OPTION" in
    c)
      Cvalue="$OPTARG"
      clustering_id=$OPTARG
      ;;
    h)
      echo "To run script please provide all required arguments and paths" >&2
      echo "script usage: $0 [-I path to long read file] [-D path to database] [-N sample name] [-c (optional) clustering identity, default is 0.99]" >&2
      exit 1
      ;;
    I)
      Ivalue="$OPTARG"
      input=$OPTARG
      ;;
    D)
      Dvalue="$OPTARG"
      db=$OPTARG
      ;;
    N)
      Nvalue="$OPTARG"
      sample_name=$OPTARG
      ;;
  esac
done

if [ ! "$Ivalue" ] || [ ! "$Dvalue" ] || [ ! "$Nvalue" ]
then
    echo "To run script please provide all required arguments and paths" >&2
    echo "script usage: $0 [-I path to long read file] [-D path to database] [-N sample name] [-c (optional) clustering identity, default is 0.99]" >&2
    exit 1
fi

echo -e "-----------Creating sequences to enrich database-----------\n" > $sample_name"_log.txt"

echo "Analyzing sample: $input" >> $sample_name"_log.txt"
echo "Reference database: $db" >> $sample_name"_log.txt"
echo "Sample name: $sample_name" >> $sample_name"_log.txt"
echo "Clustering identity: $clustering_id" >> $sample_name"_log.txt"

##1 - Clustering each sample to remove redundant sequences

echo -e "\n-----------Clustering sample to remove redundancy-----------\n" >> $sample_name"_log.txt"

out_file=$sample_name"_clustered.99.clst"
cd-hit-est -i $input -o $out_file -c $clustering_id -n 10 -d 0 -M 0 -r 0 &>> $sample_name"_log.txt"

if [ $? -eq 0 ] ; then echo -e "Finished clustering: $sample_name" >> $sample_name"_log.txt" ; else echo "Failed clustering $sample_name, check log" >>$sample_name"_log.txt" && exit  ; fi

#reformating cluster table to txt fromat
clstr2txt.pl $sample_name"_clustered.99.clst.clstr" > $sample_name"_clustered.99.clst.clstr.txt" 



##2 - Clustering non-redundant sequences to enrich database with non matching sequences
echo -e "\n-----------Clustering non-redundant sequences against the database-----------\n" >> $sample_name"_log.txt"
out_file_2d=$sample_name"_clustered2d_silva.fasta"
cd-hit-est-2d \
    -i $db \
    -i2 $out_file \
    -o $out_file_2d -c $clustering_id -n 10 -g 1 -d 0 -M 0 -r 0 -T 0 &>> $sample_name"_log.txt"

if [ $? -eq 0 ] ; then echo -e "Finished comparing: $sample_name to the database" >> $sample_name"_log.txt" ; else echo "Failed comparing $sample_name to the database, check log" >>$sample_name"_log.txt" && exit  ; fi
clstr2txt.pl $sample_name"_clustered.99.clst.clstr" > $sample_name"_clustered2d_silva.fasta.clst.clstr.txt"


echo -e "All non-redundant sequences not clustered with the database are found in the file: $out_file_2d" >> $sample_name"_log.txt"
