# C.Origami genomic features preprocessing

You will need a bigwig file of the corresponding CTCF ChIP-seq and ATAC-seq. We recommend using [Seq-N-Slide](https://igordot.github.io/sns/) pipeline for processing raw fastqs into bigwigs using the atac or chip route.

## CTCF ChIP-seq pipeline

CTCF pipeline (found under `./src/corigami/preprocessing/ctcf/ctcf_pipeline.sh`) can be run with the following steps:

1. Edit the `ctcf_ip.txt` and `ctcf_input.txt` files that links to immunoprecipitation(IP) and input fastq.gz file download links. Alternatively, you could download your own files and copy/softlink them to the corresponding location.
2. Run `ctcf_pipeline` and this will create directories to store the fastq.gz files, merge, subsample and submit a [sns-chip](https://igordot.github.io/sns/routes/chip.html) job. 
3. After the sns job is done, run `ctcf_norm.sh` to get the normalized CTCF ChIP-seq profile under name `./sns/BIGWIG/final.bigwig`. **This file is the CTCF bigwig file that we use as C.Origami input**.



## ATAC-seq pipeline 
ATAC pipeline (found under `./src/corigami/preprocessing/atac/atac_pipeline.sh`) can be run with the following steps:


1. Edit all four `rep{n}-end{m}.txt` (n, m = 1 or 2) files that indicates fastq files download links of two replicates and two read ends.
2. Run `atac_pipeline` to generate [sns-atac](https://igordot.github.io/sns/routes/atac.html) output. The target bigwig is named `./sns/BIGWIG/sub_merged.bin1.rpkm.bw`. **This file is the ATAC-seq bigwig file that we use as C.Origami input**.
