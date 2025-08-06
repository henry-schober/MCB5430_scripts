#!/bin/bash

#SBATCH --job-name=post_ChIP.sh     # name for job
#SBATCH -N 1                    # number of nodes (always 1)
#SBATCH -n 1                    # number of jobs / tasks (always 1)
#SBATCH -c 4                    # number of cores (1-4
#SBATCH -p general        # SLURM partition (always mcbstudent)
#SBATCH --qos=general        # SLURM Quality of service (always mcbstudent)
#SBATCH --mem=24G                # RAM (memory) requested
#SBATCH --mail-type=ALL
#SBATCH --mail-user=henry.schober@uconn.edu
#SBATCH -o post_ChIP.sh_%j.out
#SBATCH -e post_ChIP.sh_%j.err


module load macs
module load meme/4.12.0
module load bedtools
module load mast
module load fimo

#setting up directorys and variables
#input path for raw data
inPATH="/home/FCAM/mcb5430/usr10/midterm/"
#output path for aligned reads
outPATH="/home/FCAM/mcb5430/usr10/midterm/aligned/"
#genome variable
bt2_hg38="/home/FCAM/mcb5430/MCB5430/genomes/hg38/hg38_base_bt2/hg38_base"
#chromInfo
chromInfo="/home/FCAM/mcb5430/MCB5430/genomes/hg38/hg38_chromInfo_base.txt"
#calling peaks path
peaks="/home/FCAM/mcb5430/usr10/midterm/peaks/"

hg38_fasta="/home/FCAM/mcb5430/MCB5430/genomes/hg38/fasta/hg38.fa"
background="/home/FCAM/mcb5430/MCB5430/genomes/hg38/hg38_markov_bkgrnd.txt"
chr10="/home/FCAM/mcb5430/MCB5430/genomes/hg38/fasta/ind_chr/fasta/chr10.fa"

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


#####################




#bedtools and sort to get 100 bp around top 500 binding sites
bedtools slop -i ${peaks}summits_J.bed -g ${chromInfo} -b 50 > ${peaks}summits_J_100bp.bed
sort -n -r -k5 ${peaks}summits_J_100bp.bed | head -500 > ${peaks}summits_J_100bp_top500.bed
#convert to fasta format for meme analysis
bedtools getfasta -name -fi ${hg38_fasta} -bed ${peaks}summits_J_100bp_top500.bed -fo ${peaks}summits_J_100bp_top500.fa

#meme
meme ${peaks}summits_J_100bp_top500.fa -oc meme_out -dna -nmotifs 2 -minw 6 -maxw 12 -revcomp -mod zoops

#mast
mast meme_out/meme.txt ${peaks}summits_J_100bp_top500.fa -oc mast_out

#for finding how many of our peaks have that motif, we can search through the whole file selecting only the peaks that have X motif, then can take number of motifs in peaks divided by the total number of peaks
#I apparently need to make fasta files for all peaks, unsure what the use of this is for yet

#for 3.2

bedtools getfasta -name -fi ${hg38_fasta} -bed ${peaks}summits_J_100bp.bed -fo ${peaks}All_summits_J_100.fa

mast -hit_list meme_out/meme.txt  ${peaks}All_summits_J_100.fa -oc mast_out > mast_out/mast_summits_hits.txt

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

database="/home/FCAM/mcb5430/MCB5430/TFdb/JASPAR2018_CORE_non-redundant_pfms.meme"

tomtom -eps -oc tomtom_out meme_out/meme.txt ${database}
