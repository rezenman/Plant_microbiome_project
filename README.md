# Plant_microbiome_project

This page lists all scripts used to analyze data for the plant microbiome project managed by the Reich lab in the Weizmann institute of science in Rehovot, Israel.  
In the folder scripts you will find all kinds of scripts used for different analysis steps, here I will list each script purpose and an example of how to run it properly, please read carefully before starting.

## Enriching the database with long reads   
Use the: __enrich_db.sh__ script  
Two files are required:  
- Long read file in a fasta/fastq format(can be compressed)  
- Database to enrich  

To use the script: 
```
./enrich_db.sh -I "fasta/fastq input file" -D "database to compare to in fasta format" -N "sample name" -c "OPTIONAL - clustering identity to use with cd-hit"
```

__NOTE__ - the three arguments -I, -D, -N are mandatory, otherwise the script won't run  
__NOTE__ - make sure to allocate enough memory, the larger the database the more memory needed

