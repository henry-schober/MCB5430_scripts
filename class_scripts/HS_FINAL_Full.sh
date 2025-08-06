#!/bin/bash

#SBATCH --job-name=Full_final.sh     # name for job
#SBATCH -N 1                    # number of nodes (always 1)
#SBATCH -n 1                    # number of jobs / tasks (always 1)
#SBATCH -c 4                    # number of cores (1-4
#SBATCH -p general        # SLURM partition (always mcbstudent)
#SBATCH --qos=general        # SLURM Quality of service (always mcbstudent)
#SBATCH --mem=24G                # RAM (memory) requested
#SBATCH --mail-type=ALL
#SBATCH --mail-user=henry.schober@uconn.edu
#SBATCH -o full_final.sh_%j.out
#SBATCH -e full_final.sh_%j.err

#This is my bash code for the final, which includes everything from RNA seq to ChIP seq

#Script for RNA Seq

#modules
module load fastqc
module load cutadapt
module load hisat2
module load samtools
module load stringtie

# Part 1 - RNA Sequencing

#directory for RAW reads
fastqRAW="/home/FCAM/mcb5430/usr10/final/fastq_in/"
#directory for Processed Reads
fastqPROC="/home/FCAM/mcb5430/usr10/final/fastq_out/"
#directory for Hisat2 index for human genome
hg19_HISAT="/isg/shared/databases/alignerIndex/animal/hg19_ucsc/hg19_HISAT2/hg19"
#directory for HG19 gtf file
hg19_gtf="/isg/shared/databases/alignerIndex/animal/hg19_ucsc/hg19_HISAT2/hg19.gtf"
#directory for gtf file we want to use
hg19_gencode="/home/FCAM/mcb5430/MCB5430/MCB5430_final_2024/annotations/hg19_gencode_v19.gtf"
#given merged and annotated directory
merged_gtf="/home/FCAM/mcb5430/MCB5430/MCB5430_final_2024/DE_tables/gtf/E2_merged.annotated.gtf"


# make out directory if we don't have one already
if [ ! -d "$fastqPROC" ] 
    then mkdir -p "$fastqPROC"
fi

#clears merged list for stringtie
rm mergedlist.txt

# 1.1 Create a Pipeline to process the RNA-seq data

fileList=("E2_rep1.fastq" "E2_rep2.fastq" "untr_rep1.fastq" "untr_rep2.fastq")
for file in "${fileList[@]}"
    do 
        #getting prefixes for easier naming

        prefix=$(echo "$file" | cut -d "." -f 1)
        
        #Fastqc on RAW data

        fastqc ${fastqRAW}${prefix}.fastq
        echo "Finished original fastqc"

        #trim reads by removing bases with a quality score lower then 30
        cutadapt -j 2 -q 30 ${fastqRAW}${prefix}.fastq > ${fastqPROC}${prefix}_trimmed.fastq
        echo "Finished trimming the reads"

        #run fastqc on the trimmed reads
        fastqc ${fastqPROC}${prefix}_trimmed.fastq
        echo "Finished fastqc on trimmed reads"

        dataList=("${fastqPROC}${prefix}_trimmed.fastq" "${fastqRAW}${prefix}.fastq")
            for y in ${dataList[@]}
                do
                    #make another prefix to clean up file names

                    prefix2=$(echo "$y" | cut -d "." -f 1)

                    #hisat2
                    hisat2 -p 4 -u 15000000 --rna-strandness F -x ${hg19_HISAT} -U ${prefix2}.fastq -S ${prefix2}.sam

                    #sam to bam and sort
                    samtools sort -@ 4 -o ${prefix2}_sorted.bam ${prefix2}.sam

                    #GTF annnotation
                    stringtie -p 4 -G ${hg19_gencode} -o ${prefix2}.gtf -l {$prefix2} ${prefix2}_sorted.bam
                    
                    #adding these file names to a text file
                    echo "${prefix2}.gtf" >> mergedlist.txt
            done
done


#merging the transcripts
stringtie --merge -p 8 -G ${hg19_gencode} -o stringtie_merged.gtf mergedlist.txt

# quanitfying using stringtie
for x in $(cat mergedlist.txt)
    do
        prefix3=$(echo "$x" | cut -d "." -f 1 | rev | cut -d "/" -f 1 | rev)
        stringtie -e -B -p 8 -G ${merged_gtf} -o ./tables/${prefix3}/${prefix3}_est.gtf -A ${prefix3}_abun.tab ${prefix2}_sorted.bam
done

# Generating Bedgraph Files

########################################
###########  Bedgraph Script  ##########
########################################

module load bedtools/2.26.0

#directory for genome
hg19_gen="/home/FCAM/mcb5430/MCB5430/MCB5430_final_2024/genomeInfo/hg19_chromInfo.txt"

#create clean genome file
awk '{print $1 "\t" $2}' ${hg19_gen} > hg19_chromInfo_clean.txt


# using genome coverage to make a bedgraph

for x in $(cat mergedlist.txt)
    do
        #getting prefixes
        prefix3=$(echo "$x" | cut -d "." -f 1 )
        prefix4=$(echo "$prefix3" | rev | cut -d "/" -f 1 | rev)
        
        #creating the bedgraph files
        bedtools genomecov -bga -split -ibam ${prefix3}_sorted.bam -g hg19_chromInfo_clean.txt > ${prefix3}.bedgraph

        #adding tracklines using an if else statement
        if [[ $prefix4 == E2* ]] 
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

###############################################
##########    ChIP Seq script    ##############
###############################################

module load bedtools

#Creating the TSS files

#bed file for TSS file
hg19_bed="/home/FCAM/mcb5430/MCB5430/MCB5430_final_2024/annotations/gencode_hg19.gtf.bed"
#Summits bed file
summit_E2="/home/FCAM/mcb5430/MCB5430/MCB5430_final_2024/ER_peaks/ER_E2_summits.bed"

#create a unique TSS bed file depending on the strand type
awk '{OFS="\t"} { if (NR<2) print $0 ; 
else if ($6 == "+") print $1, $2, $2, $4, $5, $6 ; 
else print $1, $3, $3, $4, $5, $6}' < ${hg19_bed} > hg19_geneTSS.bed

#sort this bed file for bedtools closest command
sort -k1,1 -k2,2n hg19_geneTSS.bed > hg19_TSS_sorted.bed

#bedtools closest command
bedtools closest -D b -a ${summit_E2} -b hg19_TSS_sorted.bed > E2_summits_closestTSS.txt


