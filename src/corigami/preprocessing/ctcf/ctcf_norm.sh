#!/bin/bash
#SBATCH -J ctcf_norm
#SBATCH --output=./ctcf_norm_%j.log 
#SBATCH --mem=10G
#SBATCH --time=12:00:00
#SBATCH -c8

module load deeptools/3.2.1

bigwigCompare -b1 sns/BIGWIG/sub_ip.bin1.rpkm.bw -b2 sns/BIGWIG/sub_input.bin1.rpkm.bw -p 8 -o sns/BIGWIG/final.bigwig
