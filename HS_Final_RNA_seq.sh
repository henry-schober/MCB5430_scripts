#!/bin/bash

#SBATCH --job-name=RNA_seq.sh     # name for job
#SBATCH -N 1                    # number of nodes (always 1)
#SBATCH -n 1                    # number of jobs / tasks (always 1)
#SBATCH -c 4                    # number of cores (1-4
#SBATCH -p general        # SLURM partition (always mcbstudent)
#SBATCH --qos=general        # SLURM Quality of service (always mcbstudent)
#SBATCH --mem=24G                # RAM (memory) requested
#SBATCH --mail-type=ALL
#SBATCH --mail-user=henry.schober@uconn.edu
#SBATCH -o RNA_seq.sh_%j.out
#SBATCH -e RNA_seq.sh_%j.err


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

# Generate Bedgraph files







