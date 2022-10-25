#!/bin/bash

python inference/editing.py \
    --chr "chr2" \
    --celltype "imr90" \
    --start 500000 \
    --model "../../data/model/epoch=159-step=95199.ckpt" \
    --seq "../../data/data/hg38/dna_sequence" \
    --ctcf "../../data/data/hg38/imr90/genomic_features/ctcf_log2fc.bw" \
    --atac "../../data/data/hg38/imr90/genomic_features/atac.bw" \
    --del-start 1500000 \
    --del-width 100000 \
