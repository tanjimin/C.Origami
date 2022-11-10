#!/bin/bash

python corigami/inference/screening.py \
    --chr "chr2" \
    --celltype "imr90" \
    --model "../corigami_data/model_weights/corigami_base.ckpt" \
    --seq "../corigami_data/data/hg38/dna_sequence" \
    --ctcf "../corigami_data/data/hg38/imr90/genomic_features/ctcf_log2fc.bw" \
    --atac "../corigami_data/data/hg38/imr90/genomic_features/atac.bw" \
    --screen-start 1250000 \
    --screen-end 2250000 \
    --perturb-width 1000 \
    --step-size 1000 \
    --plot-impact-score \
    --save-pred --save-perturbation --save-diff --save-bedgraph

