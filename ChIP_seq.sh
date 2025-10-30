#!/bin/bash

#SBATCH --job-name=ChIP_seq.sh     # name for job
#SBATCH -N 1                    # number of nodes (always 1)
#SBATCH -n 1                    # number of jobs / tasks (always 1)
#SBATCH -c 4                    # number of cores (1-4
#SBATCH -p general        # SLURM partition (always mcbstudent)
#SBATCH --qos=general        # SLURM Quality of service (always mcbstudent)
#SBATCH --mem=24G                # RAM (memory) requested
#SBATCH --mail-type=ALL
#SBATCH --mail-user=henry.schober@uconn.edu
#SBATCH -o ChIP_seq.sh_%j.out
#SBATCH -e ChIP_seq.sh_%j.err



#subsetting the data for testing

gunzip chipJ.fastq.gz 
head -200000 chipJ.fastq > chipJ_head.fastq
gunzip inputJ.fastq.gz
head -200000 inputJ.fastq > inputJ_head.fastq

#setting up directorys and variables
#input path for raw data
inPATH="/home/FCAM/mcb5430/usr10/midterm/"
#output path for aligned reads
outPATH="/home/FCAM/mcb5430/usr10/midterm/aligned/"
#genome variable
bt2_hg38="/home/FCAM/mcb5430/MCB5430/genomes/hg38/hg38_base_bt2/hg38_base"
#chromInfo
chromInfo="/home/FCAM/mcb5430/MCB5430/genomes/hg38/hg38_chromInfo_base.txt"
#peaks
peaks="/home/FCAM/mcb5430/usr10/midterm/peaks/"
#fasta files
hg38_fasta="/home/FCAM/mcb5430/MCB5430/genomes/hg38/fasta/hg38.fa"
#background files
background="/home/FCAM/mcb5430/MCB5430/genomes/hg38/hg38_markov_bkgrnd.txt"
#fasta file for specifically chr10
chr10="/home/FCAM/mcb5430/MCB5430/genomes/hg38/fasta/ind_chr/fasta/chr10.fa"
#JASPER data base
database="/home/FCAM/mcb5430/MCB5430/TFdb/JASPAR2018_CORE_non-redundant_pfms.meme"

#for loop for script, starts off with creating a prefix for better naming

#for example cp $x > $prefix.fastqc should rename it with the new prefix

module load fastqc
module load bowtie2
module load cutadapt
module load samtools
module load bedtools
module load macs
module load meme/4.12.0


#for loop for pipeline
fileList=("chipJ.fastq" "inputJ.fastq")
for file in "${fileList[@]}"
    do 
        #getting prefixes for easier naming

        prefix=$(echo "$file" | cut -d "." -f 1)
        
        #Fastqc on data

        fastqc ${inPATH}${prefix}.fastq
        echo "Finished original fastqc"

        #clipping ChIP adapter sequence

        cutadapt -j 2 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC ${inPATH}${prefix}.fastq > ${outPATH}${prefix}_clipped.fastq
        echo "Finished clipping out the sequence"
        fastqc ${outPATH}${prefix}_clipped.fastq
        echo "Finished fastqc on clipped reads"

        #trimming the reads to a min of 30bp and qs of 30

        cutadapt -j 2 -q 30 -m 30 ${outPATH}${prefix}_clipped.fastq > ${outPATH}${prefix}_clipped_trimmed.fastq
        echo "Finished trimming the reads"
        fastqc ${outPATH}${prefix}_clipped_trimmed.fastq
        echo "Finished fastqc on trimmed and clipped reads"

        dataList=("${outPATH}${prefix}_clipped_trimmed.fastq" "${inPATH}${prefix}.fastq") #added all of the previous fastqc files to an array so that I can make a for loop
            for y in ${dataList[@]} #for loop that calls every item from the array and runs bowtie2 on these items
                do
                    bowtie2 -x ${bt2_hg38} -U ${y} -S ${y}.sam 2>&1 | tee ${y}_bowtie_log.txt #added a log file to show the bowtie output
                    echo "${y} is mapped and sam file is created"

                    samtools view -S -b ${y}.sam > ${y}.bam

                    bedtools bamtobed -i ${y}.bam > ${y}.bed

                    sortBed -i ${y}.bed > ${y}_sorted.bed

                    bedtools genomecov -bg -i ${y}_sorted.bed -g ${chromInfo} > ${y}.bedgraph
                done
    done
#adding Track lines
        
awk 'BEGIN {  print "browser position chr11:5,289,521-5,291,937"
    print "track type=bedGraph name=\"chipJ_processed_bedgraph\" description=\"chipJ_processed_bedgraph\"  visibility=full autoScale=on alwaysZero=on color=0,125,0 windowingFunction=maximum"} 
    {  print $0}' ${outPATH}chipJ_clipped_trimmed.fastq.bedgraph > ${outPATH}chipJ_processed_header.bedgraph

        #When adding in the genome browser, remember to change it hg38
    
awk 'BEGIN {  print "browser position chr11:5,289,521-5,291,937"
    print "track type=bedGraph name=\"chipJ_raw_bedgraph\" description=\"chipJ_raw_bedgraph\"  visibility=full autoScale=on alwaysZero=on color=0,125,0 windowingFunction=maximum"}
    {  print $0}' ${inPATH}chipJ.fastq.bedgraph > ${inPATH}chipJ_raw_header.bedgraph

awk 'BEGIN {  print "browser position chr11:5,289,521-5,291,937"
    print "track type=bedGraph name=\"inputJ_processed_bedgraph\" description=\"inputJ_processed_bedgraph\"  visibility=full autoScale=on alwaysZero=on color=100,50,0 windowingFunction=maximum"}
    {  print $0}' ${outPATH}inputJ_clipped_trimmed.fastq.bedgraph > ${outPATH}inputJ_processed_header.bedgraph

awk 'BEGIN {  print "browser position chr11:5,289,521-5,291,937"
    print "track type=bedGraph name=\"inputJ_raw_bedgraph\" description=\"inputJ_raw_bedgraph\"  visibility=full autoScale=on alwaysZero=on color=100,50,0 windowingFunction=maximum"}
    {  print $0}' ${inPATH}inputJ.fastq.bedgraph > ${inPATH}inputJ_raw_header.bedgraph

################## VERY IMPORTANT SO I DONT FORGET ##########
#the prefix is not linked to the original data files, so what I want to do is use the prefix for naming and making directories and file, never for fastqc or calling upon bt2


#adding tracklines to summits and peaks bed files

awk 'BEGIN {  print "browser position chr11:5,289,521-5,291,937"
    print "track type=bed name=\"summits_bed\" description=\"summits_bed\" visibility=squish autoScale=on colorByStrand=\"255,0,0 0,0,255\""} 
    {  print $0}' ${peaks}summits_J.bed > ${peaks}summits_J_header.bed

awk 'BEGIN {  print "browser position chr11:5,289,521-5,291,937"
    print "track type=bed name=\"peaks_bed\" description=\"peaks_bed\"  visibility=squish autoScale=on colorByStrand=\"255,0,0 0,0,255\""} 
    {  print $0}' ${peaks}peaks_J.bed > ${peaks}peaks_J_header.bed

#MACS

#processed
macs -f "SAM" -t ${outPATH}chipJ_clipped_trimmed.fastq.sam -c ${outPATH}inputJ_clipped_trimmed.fastq.sam -n ChIPJ_processed -g hs

#raw
macs -f "SAM" -t ${inPATH}chipJ.fastq.sam -c ${inPATH}inputJ.fastq.sam -n ChIPJ_raw -g hs

#adding tracklines
peaksList=("${inPATH}ChIPJ_processed_peaks.bed" "${inPATH}ChIPJ_raw_peaks.bed" "${inPATH}ChIPJ_processed_summits.bed" "${inPATH}ChIPJ_raw_summits.bed")

for x in ${peaksList[@]}
    do

        awk -v var="${x}" 'BEGIN {  print "browser position chr1:1,505,048-1,516,424"
              print "track type=bed name=\"" var "\" description=\"" var "\" visibility=squish autoScale=on useScore=1"}
           {  print $0}' ${x} > ${x}_header.bed
done
#MEME


#bedtools and sort to get 100 bp around top 500 binding sites
bedtools slop -i ${peaks}summits_J.bed -g ${chromInfo} -b 50 > ${peaks}summits_J_100bp.bed
sort -n -r -k5 ${peaks}summits_J_100bp.bed | head -500 > ${peaks}summits_J_100bp_top500.bed
#convert to fasta format for meme analysis
bedtools getfasta -name -fi ${hg38_fasta} -bed ${peaks}summits_J_100bp_top500.bed -fo ${peaks}summits_J_100bp_top500.fa

#meme
meme ${peaks}summits_J_100bp_top500.fa -oc meme_out -dna -nmotifs 2 -minw 6 -maxw 12 -revcomp -mod zoops

#mast
mast meme_out/meme.txt ${peaks}summits_J_100bp_top500.fa -oc mast_out



#for 3.2

#get fasta file for ALL_summits_100bp
bedtools getfasta -name -fi ${hg38_fasta} -bed ${peaks}summits_J_100bp.bed -fo ${peaks}All_summits_J_100.fa

mast -hit_list meme_out/meme.txt  ${peaks}All_summits_J_100.fa -oc mast_out > mast_out/mast_summits_hits.txt
#gets us how many motifs we have
awk '{print $3}' mast_out/mast_summits_hits.txt | sort | uniq -c > uniq_motifs.txt


#fimo

fimo --oc fimo_out --bgfile ${background} meme_out/meme.txt ${chr10}

#converting fimo out to a bed file
awk 'BEGIN {OFS="\t"} NR>1 {print $3, $4, $5, $2, $8, $6}' fimo_out/fimo.txt > fimo_out/fimo.bed

# Intersecting fimo bed file with peaks file

bedtools intersect -wo -a fimo_out/fimo.bed -b ${peaks}peaks_J.bed -g ${chromInfo} > hg38_chr10_motifs_in_peaks.txt

#How many motifs in peaks?
wc -l hg38_chr10_motifs_in_peaks.txt #1170

#how many motifs out of peaks?
wc -l fimo_out/fimo.bed #40,058, so if we do 40,058 - 1170 we get 38,888 motifs outside the peaks

#tomtom

tomtom -eps -oc tomtom_out meme_out/meme.txt ${database}