###################
# RNASeq Pipeline #
###################

## Genome
GENOME=.../GRCh37/Sequence/STAR
GENOMEGTF=.../GRCh37/Annotation/Genes/genes.gtf

## Softwre
FASTQC=/home/software/FastQC/fastqc -t 12 --noextract
TRIM=/home/software/anaconda2/bin/java -Xmx150g -jar /home/software/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 12 -phred33

STAR=/home/software/STAR/bin/Linux_x86_64/STAR -- readFilesCommand zcat --twopassMode Basic --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 0 --outFilterMismatchNmax 2 --outReadsUnmapped None --chimSegmentMin 12 --chimJunctionOverhangMin 12 --alignSJDBoverhangMin 10--alignMatesGapMax 100000 --alignIntronMax 100000 --chimSegmentReadGapMax parameter 3 --alignSJstitchMismatchNmax 5 -1 5 5 --runThreadN 12 --limitBAMsortRAM 31532137230 --sjdbGTFfile /storage01/data/iGenomes/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf --quantMode GeneCounts --outSAMtype BAM SortedByCoordinate

# SAMTOOLS
SAMVIEW=/home/software/anaconda2/bin/samtools view -@ 12
SAMSORT=/home/software/anaconda2/bin/samtools sort -@ 12
SAMRMDUP=/home/software/anaconda2/bin/samtools rmdup
SAMINDEX=/home/software/anaconda2/bin/samtools index
MPILEUP=/home/software/anaconda2/bin/samtools mpileup -B -C50 -f $(GENOME) -q 1
STATS=/home/software/anaconda2/bin/samtools stats

# PICARD INSERT SIZE ANALYSIS
PICARDINSERT="/home/software/anaconda2/bin/picard collectInsertSizeMetrics "

DATADIR=/$(MYFILE)/fastq
MYDIR=/mappingOnGenes/$(MYFILE)/

makedir : 
	if [ ! -d "$(MYDIR)" ]; then	mkdir $(MYDIR); fi

# FASTQC 1

$(MYDIR)$(MYFILE)_R1_fastqc.zip : $(DATADIR)/$(MYFILE)_R1.fastq.gz
	$(FASTQC) $(DATADIR)/$(MYFILE)_R1.fastq.gz -o $(MYDIR)
	
$(MYDIR)$(MYFILE)_R2_fastqc.zip : $(DATADIR)/$(MYFILE)_R2.fastq.gz
	$(FASTQC) $(DATADIR)/$(MYFILE)_R2.fastq.gz -o $(MYDIR)
	
# Trimming

$(MYDIR)$(MYFILE)_output1_trimmed.fastq.gz $(MYDIR)$(MYFILE)_output2_trimmed.fastq.gz: $(DATADIR)/$(MYFILE)_R1.fastq.gz $(DATADIR)/$(MYFILE)_R2.fastq.gz
	$(TRIM) $(DATADIR)/$(MYFILE)_R1.fastq.gz $(DATADIR)/$(MYFILE)_R2.fastq.gz $(MYDIR)$(MYFILE)_output1_trimmed.fastq.gz $(MYDIR)$(MYFILE)_unpaired_output1_trimmed.fastq.gz $(MYDIR)$(MYFILE)_output2_trimmed.fastq.gz $(MYDIR)$(MYFILE)_unpaired_output2_trimmed.fastq.gz ILLUMINACLIP:/home/software/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 HEADCROP:3 TRAILING:20 MINLEN:35
	
# FASTQC 2

$(MYDIR)$(MYFILE)_output1_trimmed.fastqc.zip : $(MYDIR)$(MYFILE)_output1_trimmed.fastq.gz
	$(FASTQC) $(MYDIR)$(MYFILE)_output1_trimmed.fastq.gz -o $(MYDIR)

$(MYDIR)$(MYFILE)_output2_trimmed.fastqc.zip : $(MYDIR)$(MYFILE)_output2_trimmed.fastq.gz
	$(FASTQC) $(MYDIR)$(MYFILE)_output2_trimmed.fastq.gz -o $(MYDIR)
	
# Alignment

$(MYDIR)$(MYFILE)_ReadsPerGene.out.tab : $(MYDIR)$(MYFILE)_output1_trimmed.fastq.gz $(MYDIR)$(MYFILE)_output2_trimmed.fastq.gz
	$(STAR) --genomeDir $(GENOME) --readFilesIn $(MYDIR)$(MYFILE)_output1_trimmed.fastq.gz $(MYDIR)$(MYFILE)_output2_trimmed.fastq.gz --outFileNamePrefix $(MYDIR)$(MYFILE)


all: makedir $(MYDIR)$(MYFILE)_R1_fastqc.zip $(MYDIR)$(MYFILE)_R2_fastqc.zip $(MYDIR)$(MYFILE)_output1_trimmed.fastq.gz $(MYDIR)$(MYFILE)_output2_trimmed.fastq.gz $(MYDIR)$(MYFILE)_output1_trimmed.fastqc.zip $(MYDIR)$(MYFILE)_output2_trimmed.fastqc.zip $(MYDIR)$(MYFILE)_ReadsPerGene.out.tab 