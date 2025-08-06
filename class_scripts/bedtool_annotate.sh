#!/bin/bash

#SBATCH --job-name=bed_anno     # name for job
#SBATCH -N 1                    # number of nodes (always 1)
#SBATCH -n 1                    # number of jobs / tasks (always 1)
#SBATCH -c 4                    # number of cores (1-4
#SBATCH -p general        # SLURM partition (always mcbstudent)
#SBATCH --qos=general        # SLURM Quality of service (always mcbstudent)
#SBATCH --mem=24G                # RAM (memory) requested
#SBATCH --mail-type=ALL
#SBATCH --mail-user=henry.schober@uconn.edu
#SBATCH -o bed_anno_%j.out
#SBATCH -e bed_anno_%j.err

module load bedtools

bedtools annotate -i ATF1_summits_chr1.bed -files hg38_promoters.bed hg38_geneBody.bed hg38_intergenic.bed -names promoter gene intergenic > ATF1_summits_overlap.txt

bedtools annotate -i fimo_ATF1_chr1.bed -files hg38_promoters.bed hg38_geneBody.bed hg38_intergenic.bed -names promoter gene intergenic > ATF1_motifs_overlap.txt
