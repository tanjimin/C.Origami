#!/bin/bash
#SBATCH -J ctcf_pipeline
#SBATCH --output=./ctcf_pipeline_%j.log 
#SBATCH --mem=20G
#SBATCH --time=48:00:00
#SBATCH -N 8

module load seqtk/1.3
module load git
module load samtools/1.9

echo "Downloading"
## Download
# Create dir
mkdir -p fastq
cd fastq

mkdir -p ip
cd ip

# Download fastq
xargs -L 1 curl -s -O -J -L < ../../ctcf_ip.txt

cd ../
mkdir -p input
cd input

# Download fastq
xargs -L 1 curl -s -O -J -L < ../../ctcf_input.txt

cd ../../

echo "All files downloaded, merging"
# Merge fastq files
cat ./fastq/ip/* > ./fastq/ip/ip.fastq.gz
cat ./fastq/input/* > ./fastq/input/input.fastq.gz

echo "Files merged, subsampling"
# subsample fastq files
zcat ./fastq/ip/ip.fastq.gz | echo "IP reads: $((`wc -l`/4))"
seqtk sample -s 2021 ./fastq/ip/ip.fastq.gz 30000000 | gzip > ./fastq/ip/sub_ip_R1.fastq.gz
zcat ./fastq/ip/sub_ip_R1.fastq.gz | echo "sub IP reads: $((`wc -l`/4))"

zcat ./fastq/input/input.fastq.gz | echo "input reads: $((`wc -l`/4))"
seqtk sample -s 2021 ./fastq/input/input.fastq.gz 30000000 | gzip > ./fastq/input/sub_input_R1.fastq.gz
zcat ./fastq/input/sub_input_R1.fastq.gz | echo "sub input reads: $((`wc -l`/4))"

echo "Files merged, running pipeline"
# run sns pipeline
mkdir -p sns
cd sns
git clone --depth 1 https://github.com/igordot/sns
sns/generate-settings hg38
sns/gather-fastqs ../fastq
sns/run chip

echo "Pipeline job submitted"
