# Plant_microbiome_project
![Untitled](https://user-images.githubusercontent.com/67236735/225245574-23cfd9dd-1831-4b0b-b2e0-aa2ec11fa570.png)  
This page lists all scripts used to analyze data for the plant microbiome project managed by the Reich lab in the Weizmann institute of science in Rehovot, Israel.  
In the folder scripts you will find all kinds of scripts used for different analysis steps, here I will list each script purpose and an example of how to run it properly, please read carefully before starting.

## Enriching the database with long reads
All Relevant scripts are in the scripts/db_enrichment_all_steps folder

### First script to run
Use the: __01_prepare_raw_files.sh__ script

This script will filter the long reads fastq file and will only keep reads which match the primers that are in the primer file you will list, 
it will also reverse complement the reads so in the end all reads will have the same directionality

steps to perform
- create a primer file
- change the path to the primer file and the r script called reads_per_region.R

Run this script from the root folder where all long read files are found

### Second script to run    
Use the: __enrich_db.sh__ script  
Two files are required:  
- Long read file in a fasta/fastq format(can be compressed)  
- Database to enrich  

Packages needed: cd-hit version 4.8.1, if working in the weizmann wexac cluster just use: 
```
module load cd-hit/4.8.1
```
To use the script: 
```
./enrich_db.sh -I "fasta/fastq input file" -D "database to compare to in fasta format" -N "sample name" -c "OPTIONAL - clustering identity to use with cd-hit"
```

__NOTE__ - the three arguments -I, -D, -N are mandatory, otherwise the script won't run  
__NOTE__ - make sure to allocate enough memory, the larger the database the more memory needed

Outputs:
- First output is a fasta file containing non-redundant sequences that does not match the database  
- Second output is a log.txt file containing all verbose


## Running the UMI pipeline  
Use the: __UMI_pipeline.sh__ script  
Five inputs are required:
- Forward primers to use for region demultiplexing (fasta format, see in sample data)
- Reverse primers to use for region demultiplexing (fasta format, see in sample data)
- Read1 fastq file (can be compressed)
- Read2 fastq file (can be compressed)
- Sample name 

Packages needed: if working in the weizmann wexac cluster just use: 
```
module load cutadapt/4.2-GCCcore-11.3.0
module load cd-hit/4.8.1
module load BBMap/38.90-GCC-10.2.0
module load seqtk/1.2
module load seqkit/2.3.1
module load R/3.5.1
```
To use the script: 
```
./UMI_pipeline.sh "for_primers" "reverse primers" "read1 fastq" "read2 fastq" "sample_name"
``` 
__NOTE__ - make sure to allocate enough memory and cores, the cd-hit process can take some time, more cores are much faster

Outputs:
- Two fastq files, R1 and R2, called $sample_name"_L001_R1_001.fastq.gz" and $sample_name"_L001_R2_001.fastq.gz"  
- Log files, cutadapt_log.txt, clustering_df.txt - here you can see reads per region in the raw data and in the final files
