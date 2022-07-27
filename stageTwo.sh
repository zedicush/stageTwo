#!/usr/bin/bash
####number of sequence in DNA
wget https://raw.githubusercontent.com/HackBio-Internship/wale-home-tasks/main/DNA.fa
grep ">" DNA.fa |wc -l
##total A,T,C,G count
grep -v "^>" DNA.fa | wc -m
####setup conda
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm -rf ~/miniconda3/miniconda.sh
~/miniconda3/bin/conda init bash
~/miniconda3/bin/conda init zsh
###close and  open the terminal for the change to take effects
##setup en environment
conda create -y --name NGS
###following command allows me to activate env from bash file
source ~/miniconda3/etc/profile.d/conda.sh
conda activate NGS
###install fastqc, multiqc and fastp
conda install -c bioconda fastqc
conda install -c bioconda multiQC
conda install -c bioconda fastp
####download files
wget https://github.com/josoga2/yt-dataset/blob/main/dataset/raw_reads/SRR1552451.fastq.gz?raw=true
wget https://github.com/josoga2/yt-dataset/blob/main/dataset/raw_reads/SRR1552452.fastq.gz?raw=true
wget https://github.com/josoga2/yt-dataset/blob/main/dataset/raw_reads/Alsen_R1.fastq.gz?raw=true
wget https://github.com/josoga2/yt-dataset/blob/main/dataset/raw_reads/Alsen_R2.fastq.gz?raw=true
###rename files
mv SRR1552451.fastq.gz\?raw\=true SRR1552451.fastq.gz
mv SRR1552452.fastq.gz\?raw\=true SRR1552452.fastq.gz
mv Alsen_R1.fastq.gz\?raw\=true Alsen_R1.fastq.gz
mv Alsen_R2.fastq.gz\?raw\=true Alsen_R2.fastq.gz
### create output folder
mkdir output
####run fastqc on these files and save the results in the output file
conda activate NGS
fastqc SRR1552451.fastq.gz -O output
fastqc SRR1552452.fastq.gz -O output
fastqc Alsen_R1.fastq.gz -O output
fastqc Alsen_R2.fastq.gz -O output
###########Run multiQC
multiqc output/ -o output/
####run fastp: in1 is the input for the reverse and in2 the input for the forward. out* are the outputs and Alsen.html the html file related
fastp --in1 Alsen_R1.fastq.gz --in2 Alsen_R2.fastq.gz --out1 output/Alsen_R1.trimmed.fastq.gz --out2 output/Alsen_R2.trimmed.fastq.gz -h output/Alsen.html
