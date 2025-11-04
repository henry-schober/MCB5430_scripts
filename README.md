# Example Scripts from MCB5430

--------------------------------------------------
This is a repository I created in order to showcase some of the scripts developed during my time with the MCB5430 course at UConn. Two main pipelines are included here: [ChIP-seq](./ChIP_seq.sh) and [RNA-seq](./HS_Final_RNA_seq.sh). Other pipelines and scripts developed are inclused in the [Class_Scripts](./Class_Scripts/) folder.


## RNA Seq pipeline processes

The example RNA-seq pipeline I built uses the following tools:
- FastQC
    -  Quality control of raw reads 
- Cutadapt
    -  Adapter trimming and removing bases with a quality score lower than 30
- HISAT2
    -  Aligning reads to the human genome (hg19) index
- Samtools
    -  Converting SAM files to BAM files, sorting and indexing BAM files
- StringTie
    - Used for GTF annotation and quantification of transcripts

### Bedgraphs were created in the [HS_Final_Bedgraph.sh](./class_scripts/HS_Final_Bedgraph.sh) script using the following tools:
- Bedtools
    -  Creating bedgraph files from sorted BAM files


## ChIP Seq pipeline processes
The example ChIP-seq pipeline I built uses the following tools:
- FastQC
    -  Quality control of raw reads
- Cutadapt
    - Adapter trimming of exampe adapeter sequence and removing bases with a quality score lower than 30
- Bowtie2
    - Aligning reads to the human genome (hg38) index
- Samtools
    - Creating SAM files, later converting to BAM files
- Bedtools
    - Creating bedgraph files from sorting bed files
- MACS2
    - Peak calling from sorted BAM files
- MEME
    - Motif calling
- MAST
    - Determine Motif presence in sequence data
- FIMO
    - Determine location of all occurences of motifs in sequences
- TomTom
    - Comparison of discovered motifs to known motifs in databases



### All of these visualizations can be found in [Visualizations](./Visualizations/Graphs.md).

