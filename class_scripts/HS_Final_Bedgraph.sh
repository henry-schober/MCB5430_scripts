#!/bin/bash

#SBATCH --job-name=RNA_Bedgraph.sh     # name for job
#SBATCH -N 1                    # number of nodes (always 1)
#SBATCH -n 1                    # number of jobs / tasks (always 1)
#SBATCH -c 4                    # number of cores (1-4
#SBATCH -p general        # SLURM partition (always mcbstudent)
#SBATCH --qos=general        # SLURM Quality of service (always mcbstudent)
#SBATCH --mem=24G                # RAM (memory) requested
#SBATCH --mail-type=ALL
#SBATCH --mail-user=henry.schober@uconn.edu
#SBATCH -o RNA_Bedgraph.sh_%j.out
#SBATCH -e RNA_Bedgraph.sh_%j.err

#Bedgraph Script

module load bedtools

#directory for genome
hg19_gen="/home/FCAM/mcb5430/MCB5430/MCB5430_final_2024/genomeInfo/hg19_chromInfo.txt"

# using genome coverage to make a bedgraph

for x in $(cat mergedlist.txt)
    do
        #getting prefixes
        prefix3=$(echo "$x" | cut -d "." -f 1 )
        prefix4=$(echo "$prefix3" | rev | cut -d "/" -f 1 | rev)
        
        #creating the bedgraph files
        bedtools genomecov -bga -split -i ${prefix3}_sorted.bam -g ${hg19_gen} > ${prefix3}.bedgraph

        #adding tracklines using an if else statement
        if [ $prefix4 == E2* ] 
            then
                awk -v var="${prefix4}" 'BEGIN {  print "browser position chr11:5,289,521-5,291,937"
                    print "track type=bedGraph name=\"" var "\" description=\"" var "\" visibility=full autoScale=on alwaysZero=on color=100,125,0 windowingFunction=maximum" } 
                    {  print $0}' ${prefix3}.bedgraph > ${prefix3}_header.bedgraph
        else
                awk -v var="${prefix4}" 'BEGIN {  print "browser position chr11:5,289,521-5,291,937"
                    print "track type=bedGraph name=\"" var "\" description=\"" var "\" visibility=full autoScale=on alwaysZero=on color=50,0,125 windowingFunction=maximum" } 
                    {  print $0}' ${prefix3}.bedgraph > ${prefix3}_header.bedgraph
        fi
done


