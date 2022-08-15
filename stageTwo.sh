#!/usr/bin/bash
######reproduction of the tutorial
##create directories 
 mkdir dc_workshop
 cd dc_workshop/
 mkdir -p data/ref_genome
####download Ecoli from NCBI and unzip the file
curl -L -o data/ref_genome/ecoli_rel606.fasta.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/017/985/GCA_000017985.1_ASM1798v1/GCA_000017985.1_ASM1798v1_genomic.fna.gz
gunzip data/ref_genome/ecoli_rel606.fasta.gz
###The real name of the genome is written in the first line of the fasta file. So we will try to print the first line
head data/ref_genome/ecoli_rel606.fasta 
######Download trimmed data to work with them. untar them 
curl -L -o sub.tar.gz https://ndownloader.figshare.com/files/14418248
tar xvf sub.tar.gz
#####create another directory and move the data untared in that directory
mkdir data/trimmed_fastq_small
mv sub/ ~/dc_workshop/data/trimmed_fastq_small
#Create in on commend line new directories for the results
mkdir -p results/sam results/bam results/bcf results/vcf
###index the reference genome
bwa index data/ref_genome/ecoli_rel606.fasta
####aligning the reads from just on sample
 bwa mem data/ref_genome/ecoli_rel606.fasta data/trimmed_fastq_small/sub/SRR2584866_1.trim.sub.fastq data/trimmed_fastq_small/sub/SRR2584866_2.trim.sub.fastq > results/sam/SRR2584866.aligned.sam
##########convert sam file to bam file
samtools view -S -b results/sam/SRR2584866.aligned.sam > results/bam/SRR2584866.aligned.bam
###sort BAM file by coordinates
samtools sort -o results/bam/SRR2584866.aligned.sorted.bam results/bam/SRR2584866.aligned.bam 
#####the following command help us lear more about the bam file             
samtools flagstat results/bam/SRR2584866.aligned.sorted.bam
######calculate the read coverage
###convert to BAM files
samtools view -S -b results/sam/SRR2584866.aligned.sam > results/bam/SRR2584866.aligned.bam
###sort BAM file by coordinates
samtools sort -o results/bam/SRR2584866.aligned.sorted.bam results/bam/SRR2584866.aligned.bam 
#####the following command help us lear more about the bam file             
samtools flagstat results/bam/SRR2584866.aligned.sorted.bam
######calculate the read coverage
bcftools mpileup -O b -o results/bcf/SRR2584866_raw.bcf -f data/ref_genome/ecoli_rel606.fasta results/bam/SRR2584866.al>
#######detect SNVs
bcftools call --ploidy 1 -m -v -o results/vcf/SRR2584866_variants.vcf results/bcf/SRR2584866_raw.bcf
##### the number of variants in the file
grep -v "#" results/vcf/SRR2584866_final_variants.vcf | wc -l
#### index the bam file
samtools index results/bam/SRR2584866.aligned.sorted.bam
###########visiualize our mapped read with tview
samtools tview results/bam/SRR2584866.aligned.sorted.bam data/ref_genome/ecoli_rel606.fasta



#reproduce the tutorial with the new data
#####download the files
wget https://zenodo.org/record/2582555/files/SLGFSK-N_231335_r1_chr5_12_17.fastq.gz
wget https://zenodo.org/record/2582555/files/SLGFSK-N_231335_r2_chr5_12_17.fastq.gz
wget https://zenodo.org/record/2582555/files/SLGFSK-T_231336_r1_chr5_12_17.fastq.gz
wget https://zenodo.org/record/2582555/files/SLGFSK-T_231336_r2_chr5_12_17.fastq.gz
wget https://zenodo.org/record/2582555/files/hg19.chr5_12_17.fa.gz
####unzip the files 
gunzip SLGFSK-N_231335_r1_chr5_12_17.fastq.gz 
gunzip SLGFSK-N_231335_r2_chr5_12_17.fastq.gz
gunzip SLGFSK-T_231336_r1_chr5_12_17.fastq.gz
gunzip SLGFSK-T_231336_r2_chr5_12_17.fastq.gz
gunzip hg19.chr5_12_17.fa.gz
#create a new directory and move the files there
mkdir  -p data/ref_genome
mkdir  -p data/trimmed_fastq_small
mv SLGFSK-N_231335_r1_chr5_12_17.fastq  SLGFSK-N_231335_r2_chr5_12_17.fastq  SLGFSK-T_231336_r1_chr5_12_17.fastq  SLGFSK-T_231336_r2_chr5_12_17.fastq data/trimmed_fastq_small/
mv hg19.chr5_12_17.fa data/ref_genome

####create directories for results
mkdir -p results/sam results/bam results/bcf results/vcf
#####remove the NNNN noise
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' < data/ref_genome/hg19.chr5_12_17.fa | sed 's/N//g' | tr "\t" "\n" > hg19.chr5_12_17_new.fa
####index the reference genome
 bwa index data/ref_genome/hg19.chr5_12_17_new.fa
#####################align reference genome
bwa mem -M data/ref_genome/hg19.chr5_12_17_new.fa <(cat data/trimmed_fastq_small/SLGFSK-N_231335_r1_chr5_12_17.fastq  data/trimmed_fastq_small/SLGFSK-T_231336_r1_chr5_12_17.fastq) <(cat data/trimmed_fastq_small/SLGFSK-N_231335_r2_chr5_12_17.fastq  data/trimmed_fastq_small/SLGFSK-T_231336_r2_chr5_12_17.fastq) > results/sam/reads.aligned.sam
###########################convert sam to bam
samtools view -S -b results/sam/reads.aligned.sam > results/bam/reads.aligned.bam
samtools sort -o results/bam/reads.aligned.sorted.bam results/bam/reads.aligned.bam 
samtools flagstat results/bam/reads.aligned.sorted.bam
####calculate the read coverage
bcftools mpileup -O b -o results/bcf/reads_raw.bcf -f data/ref_genome/hg19.chr5_12_17.fa  results/bam/reads.aligned.sorted.bam 

####Detect the single nucleotide variants
bcftools call --ploidy 1 -m -v -o results/vcf/reads_variants.vcf results/bcf/reads_raw.bcf 
#######################Filter and report the SNV variants in variant calling format 
vcfutils.pl varFilter results/vcf/reads_variants.vcf  > results/vcf/reads_final_variants.vcf
####display the file
less -S results/vcf/reads_final_variants.vcf 
#####viewing : indexing bam file
samtools index results/bam/reads.aligned.sorted.bam
##########visualize our mapped reads
samtools tview results/bam/reads.aligned.sorted.bam data/ref_genome/hg19.chr5_12_17_new.fa
