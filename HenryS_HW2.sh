#!/bin/bash

#SBATCH --job-name=HenryS_HW2.sh     # name for job
#SBATCH -N 1                    # number of nodes (always 1)
#SBATCH -n 1                    # number of jobs / tasks (always 1)
#SBATCH -c 1                    # number of cores (1-4
#SBATCH -p general        # SLURM partition (always mcbstudent)
#SBATCH --qos=general        # SLURM Quality of service (always mcbstudent)
#SBATCH --mem=1G                # RAM (memory) requested
#SBATCH --mail-type=ALL
#SBATCH --mail-user=henry.schober@uconn.edu
#SBATCH -o HenryS_HW2.sh_%j.out
#SBATCH -e HenryS_HW2.sh_%j.err




#Homework 2

# 1.

cp /home/FCAM/mcb5430/MCB5430/Homework/HW2/SRR412199.fastq .

wc -l SRR412199.fastq

# 19904076 lines in fast q file. You divide this by 4 in order to get the amount of sequence reads as there are four lines per sequence read.
# 4976019 sequences

# 2.

#first 2 million

head -2000000 SRR412199.fastq > SRR412199_head.fastq

#last 2 million

tail -2000000 SRR412199.fastq > SRR412199_tail.fastq

#mv SRR412199_head.fastq ../data/.
#mv SRR412199_tail.fastq ../data/.

#TGCTTGGACTACATATGGTTGAGGGTTGTATGGAATTCTCGGGTGCCAAGG

#load all of the modules you would need
module load fastqc
module load bowtie2
module load cutadapt

#create a variable to easily call the bowtie files
bt2_Dros="/isg/shared/databases/alignerIndex/animal/drosophila_melanogaster/current/bowtie2/Drosophila_melanogaster"

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
## Question 5

#For tail_original:  1.38% alignment rate, only 1179 aligned once, 5711 aligned more than once
#For tail_clipped: 2.7% alignment rate, only 1234 aligned once, 12280 aligned more than once
#For tail_trimmed: 2.87% alignment rate, only 1334 aligned once, 12926 aligned more than once

#Question 6

#After looking at the html files, I believe that all of the quality scores are all around 40. The issue with this is that from the per base sequence quality, it looks like we are not getting
#enough per base sequence quality reads. The major issue with the clipped data is that highest percentage sequence is blank(since we clipped the largest sequence), while both trimmed and original have a sequence that dominates at around 65%
#One thing that we could do for the trimmed data is to actually increase the minimum length as there is not much difference in the original data and the clipped data, there seems to be a lot of sequence length of >15 in the original data
# I would say the best way to deal with this data again is to use cutadapt on the clipped data, and get rid of sequence lengths <15 as this would get rid of the blank sequence
# Overall, I think we didn't trim enough and clipping the data just removes the sequence but replaces it with blank space.




