#!/bin/bash

python inference/screening.py \
    --chr "chr2" \
    --celltype "imr90" \
    --model "../../data/model/epoch=159-step=95199.ckpt" \
    --seq "../../data/data/hg38/dna_sequence" \
    --ctcf "../../data/data/hg38/imr90/genomic_features/ctcf_log2fc.bw" \
    --atac "../../data/data/hg38/imr90/genomic_features/atac.bw" \
    --screen-start 1250000 \
    --screen-end 2250000 \
    --perturb-width 1000 \
    --step-size 1000 \
    --plot-impact-score \
    --save-pred --save-perturbation --save-diff --save-bedgraph

