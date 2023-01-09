#!/bin/bash
#SBATCH -J atac
#SBATCH --output=./atac_pipeline_%j.log 
#SBATCH --mem=20G
#SBATCH --time=48:00:00
#SBATCH -N 8

module load seqtk/1.3
module load git
module load samtools/1.9

echo "Downloading"
## Download
# Create data dir
mkdir -p fastq
cd fastq

## Download All Files

# End 1
mkdir -p end1
cd end1

mkdir -p rep1
cd rep1
xargs -L 1 curl -s -O -J -L < ../../../rep1-end1.txt
cd ../

mkdir -p rep2
cd rep2
xargs -L 1 curl -s -O -J -L < ../../../rep2-end1.txt
cd ../

cd ../

# End 2
mkdir -p end2
cd end2

mkdir -p rep1
cd rep1
xargs -L 1 curl -s -O -J -L < ../../../rep1-end2.txt
cd ../

mkdir -p rep2
cd rep2
xargs -L 1 curl -s -O -J -L < ../../../rep2-end2.txt
cd ../

cd ../

cd ../

echo "All files downloaded, merging"
## Cat different replicates
cat ./fastq/end1/rep1/* ./fastq/end1/rep2/* > ./fastq/end1/merged1.fastq.gz
cat ./fastq/end2/rep1/* ./fastq/end2/rep2/* > ./fastq/end2/merged2.fastq.gz

echo "Files merged, subsampling end 1"
zcat ./fastq/end1/merged1.fastq.gz | echo "end 1 reads: $((`wc -l`/4))"
seqtk sample -s 2021 ./fastq/end1/merged1.fastq.gz 40000000 | gzip > ./fastq/sub_merged_R1.fastq.gz
zcat ./fastq/sub_merged_R1.fastq.gz | echo "sub end 1 reads: $((`wc -l`/4))"

echo "Files merged, subsampling end 2"
zcat ./fastq/end2/merged2.fastq.gz | echo "end 2 reads: $((`wc -l`/4))"
seqtk sample -s 2021 ./fastq/end2/merged2.fastq.gz 40000000 | gzip > ./fastq/sub_merged_R2.fastq.gz
zcat ./fastq/sub_merged_R2.fastq.gz | echo "sub end 2 reads: $((`wc -l`/4))"

echo "Files merged, running sns pipeline"
mkdir -p sns
cd sns
git clone --depth 1 https://github.com/igordot/sns
sns/generate-settings hg38
sns/gather-fastqs ../fastq
sns/run atac

echo "Pipeline job submitted"
