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

            #Running bowtie on the preprocessed and raw data with log files

        bowtie2 -x ${bt2_hg38} -U ${outPATH}${prefix}_clipped_trimmed.fastq -S ${outPATH}${prefix}_clipped_trimmed.sam 2>&1 | tee ${outPATH}${prefix}_clipped_trimmed.log

        bowtie2 -x ${bt2_hg38} -U ${inPATH}${prefix}.fastq -S ${outPATH}${prefix}.sam 2>&1 | tee ${outPATH}${prefix}.log

        # Converting sam file to bed file and sorting it
        #Processed data
        samtools view -S -b ${outPATH}${prefix}_clipped_trimmed.sam > ${outPATH}${prefix}_clipped_trimmed.bam
        #Raw Data
        samtools view -S -b ${outPATH}${prefix}.sam > ${outPATH}${prefix}.bam

        #converting to a bed file
        #processed data
        bedtools bamtobed -i ${outPATH}${prefix}_clipped_trimmed.bam > ${outPATH}${prefix}_clipped_trimmed.bed
        #Raw Data
        bedtools bamtobed -i ${outPATH}${prefix}.bam > ${outPATH}${prefix}.bed

        #sorting the bed files
        #processed data
        sortBed -i ${outPATH}${prefix}_clipped_trimmed.bed > ${outPATH}${prefix}_processed_sorted.bed
        #raw data
        sortBed -i ${outPATH}${prefix}.bed > ${outPATH}${prefix}_raw_sorted.bed

        #Creating the Bedgraph
        bedtools genomecov -bg -i ${outPATH}${prefix}_processed_sorted.bed -g ${chromInfo} > ${outPATH}${prefix}_processed.bedgraph

        #Raw data
        bedtools genomecov -bg -i ${outPATH}${prefix}_raw_sorted.bed -g ${chromInfo} > ${outPATH}${prefix}_raw.bedgraph

   done


