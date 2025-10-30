#!/bin/bash

#SBATCH --job-name=ChIP_seq.sh     # name for job
#SBATCH -N 1                    # number of nodes (always 1)
#SBATCH -n 1                    # number of jobs / tasks (always 1)
#SBATCH -c 4                    # number of cores (1-4
#SBATCH -p general        # SLURM partition (always mcbstudent)
#SBATCH --qos=general        # SLURM Quality of service (always mcbstudent)
#SBATCH --mem=1G                # RAM (memory) requested
#SBATCH --mail-type=ALL
#SBATCH --mail-user=henry.schober@uconn.edu
#SBATCH -o ChIP_seq.sh_%j.out
#SBATCH -e ChIP_seq.sh_%j.err



#subsetting the data for testing

gunzip chipJ.fastq.gz 
head -200000 chipJ.fastq > chipJ_head.fastq
gunzip inputJ.fastq.gz
head -200000 inputJ.fastq > inputJ_head.fastq

for x in $fileList do; echo ${x}; done
#create an empty array to hold prefixes
prefixes=()

#for loop for script, starts off with creating a prefix for better naming

#for example cp $x > $prefix.fastqc should rename it with the new prefix
fileList=("chipJ_head.fastq" "inputJ_head.fastq")
for file in "${fileList[@]}"
    do 
        prefix=$(echo "$file" | cut -d "." -f 1)
        prefixes+=("$prefix")
    done


################## VERY IMPORTANT SO I DONT FORGET ##########
#the prefix is not linked to the original data files, so what I want to do is use the prefix for naming and making directories and file, never for fastqc or calling upon bt2


#TGCTTGGACTACATATGGTTGAGGGTTGTATGGAATTCTCGGGTGCCAAGG

#load all of the modules you would need
module load fastqc
module load bowtie2
module load cutadapt

#create a variable to easily call the bowtie files
bt2_hg38="/home/FCAM/mcb5430/MCB5430/genomes/hg38/hg38_base_bt2/"


for x in SRR412199_head.fastq SRR412199_tail.fastq
    do
        fastqc $x #run fastqc on origin file
        echo "finished fastqc of ${x}" #progress report of finishing fastqc
        cutadapt -j 2 -a TGCTTGGACTACATATGGTTGAGGGTTGTATGGAATTCTCGGGTGCCAAGG $x > ${x}_clipped #clipping out this sequence from the original data file"
        echo "finished clipping ${x}" #progress report of the clipping
        fastqc ${x}_clipped #run fastqc on the clipped portion
        echo "finished fastqc of ${x}_clipped" #progress report of fastqc on clipped data
        cutadapt -j 2 -q 30 -m 15 $x > ${x}_trimmed #trimming the data so that the quality scores are above 30 with a minimum length of 15
        echo "finished trimming ${x}" #progress report of trimming
        fastqc ${x}_trimmed #run fastqc on trimmed
        echo "finished fastqc of ${x}_trimmed" #progress report on running fastqc on trimmed data

        dataList=("${x}" "${x}_clipped" "${x}_trimmed") #added all of the previous fastqc files to an array so that I can make a for loop
        for y in ${dataList[@]} #for loop that calls every item from the array and runs bowtie2 on these items
            do
                bowtie2 -p2 -x ${bt2_Dros} -U ${y} -S ${y}.sam 2>&1 | tee bowtie_${y}_log.txt #added a log file to show the bowtie output
                echo "${y} is mapped and sam file is created"
            done
   done
