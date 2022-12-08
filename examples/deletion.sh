#!/bin/bash

python corigami/inference/editing.py \
    --chr "chr2" \
    --celltype "imr90" \
    --start 500000 \
    --model "../corigami_data/model_weights/corigami_base.ckpt" \
    --seq "../corigami_data/data/hg38/dna_sequence" \
    --ctcf "../corigami_data/data/hg38/imr90/genomic_features/ctcf_log2fc.bw" \
    --atac "../corigami_data/data/hg38/imr90/genomic_features/atac.bw" \
    --del-start 1500000 \
    --del-width 100000 \
