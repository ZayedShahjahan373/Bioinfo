#produce quality report for reads stored in reads.fastq.gz pipe into analysis directory 
#store under preprocessing folder in the analysis directory
fastqc reads.fastq.gz > $ANADIR/qc_preprocessing

#Next step is to use trimmomatic to filter out paired-end reads whose mean base quality
#is below 20
java -jar trimmomatic-0.32.jar PE -phred64 reads1.fastq.gz paired1.fq.gz unpaired1.fq.gz paired2.fq.gz paired2.fq.gz AVGQUAL:20

#For paired-end reads trim bases from 3' end when the base quality is below 20; also filter
#out reads which are shorter than 50 bases after trimming 
#NOTE:: order of the filtering steps is defined by the command
java -jar trimmomatic-0.32.jar PE -phred64 reads1.fastq.gz paired1.fq.gz unpaired1.fq.gz paired2.fq.gz paired2.fq.gz AVGQUAL:20 TRAILING:20 MINLEN:50

#Trimming using a sliding window for trimmomatic. first numeric tells you the base length of your window
#second numeric tells you the quality threshold
java -jar trimmomatic-0.32.jar PE -phred64 reads1.fastq.gz paired1.fq.gz unpaired1.fq.gz paired2.fq.gz paired2.fq.gz SLIDINGWINDOW:3:20 MINLEN:50

#Remove Illumina TruSeq2 adapters from paired-end reads:
