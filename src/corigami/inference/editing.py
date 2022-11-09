import os
import numpy as np
import pandas as pd
import sys
import torch

import corigami.inference.utils.inference_utils as infer
from corigami.inference.utils import plot_utils 

import argparse

def main():
    parser = argparse.ArgumentParser(description='C.Origami Editing Module.')
    
    # Output location
    parser.add_argument('--out', dest='output_path', 
                        default='outputs',
                        help='output path for storing results (default: %(default)s)')

    # Location related params
    parser.add_argument('--celltype', dest='celltype', 
                        help='Sample cell type for prediction, used for output separation', required=True)
    parser.add_argument('--chr', dest='chr_name', 
                        help='Chromosome for prediction', required=True)
    parser.add_argument('--start', dest='start', type=int,
                        help='Starting point for prediction (width is 2097152 bp which is the input window size)', required=True)
    parser.add_argument('--model', dest='model_path', 
                        help='Path to the model checkpoint', required=True)
    parser.add_argument('--seq', dest='seq_path', 
                        help='Path to the folder where the sequence .fa.gz files are stored', required=True)
    parser.add_argument('--ctcf', dest='ctcf_path', 
                        help='Path to the folder where the CTCF ChIP-seq .bw files are stored', required=True)
    parser.add_argument('--atac', dest='atac_path', 
                        help='Path to the folder where the ATAC-seq .bw files are stored', required=True)

    # Deletion related params
    parser.add_argument('--del-start', dest='deletion_start', type=int,
                        help='Starting point for deletion.', required=True)
    parser.add_argument('--del-width', dest='deletion_width', type=int,
                        help='Width for deletion.', required=True)
    parser.add_argument('--padding', dest='end_padding_type', 
                        default='zero',
                        help='Padding type, either zero or follow. Using zero: the missing region at the end will be padded with zero for ctcf and atac seq, while sequence will be padded with N (unknown necleotide). Using follow: the end will be padded with features in the following region (default: %(default)s)')
    parser.add_argument('--hide-line', dest='hide_deletion_line', 
                        action = 'store_true',
                        help='Remove the line showing deletion site (default: %(default)s)')

    args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])
    single_deletion(args.output_path, args.celltype, args.chr_name, args.start, 
                    args.deletion_start, args.deletion_width, 
                    args.model_path,
                    args.seq_path, args.ctcf_path, args.atac_path, 
                    show_deletion_line = not args.hide_deletion_line,
                    end_padding_type = args.end_padding_type)

def single_deletion(output_path, celltype, chr_name, start, deletion_start, deletion_width, model_path, seq_path, ctcf_path, atac_path, show_deletion_line = True, end_padding_type = 'zero'):

    # Define window which accomodates deletion
    window = 2097152 + deletion_width
    seq_region, ctcf_region, atac_region = infer.load_region(chr_name, 
            start, seq_path, ctcf_path, atac_path, window = window)
    # Delete inputs
    seq_region, ctcf_region, atac_region = deletion_with_padding(start, 
            deletion_start, deletion_width, seq_region, ctcf_region, 
            atac_region, end_padding_type)
    # Prediction
    pred = infer.prediction(seq_region, ctcf_region, atac_region, model_path)
    # Initialize plotting class
    plot = plot_utils.MatrixPlotDeletion(output_path, pred, 'deletion', 
            celltype, chr_name, start, deletion_start, deletion_width, 
            padding_type = end_padding_type,
            show_deletion_line = show_deletion_line)
    plot.plot()

def deletion_with_padding(start, deletion_start, deletion_width, seq_region, ctcf_region, atac_region, end_padding_type):
    ''' Delete all signals at a specfied location with corresponding padding at the end '''
    # end_padding_type takes values of either 'zero' or 'follow'
    if end_padding_type == 'zero':
        seq_region, ctcf_region, atac_region = zero_region(seq_region, 
                ctcf_region, atac_region)
    elif end_padding_type == 'follow': pass
    else: raise Exception('unknown padding')
    # Deletion
    seq_region, ctcf_region, atac_region = delete(deletion_start - start, 
            deletion_start - start + deletion_width, 
            seq_region, ctcf_region, atac_region)
    return seq_region, ctcf_region, atac_region

def zero_region(seq, ctcf, atac, window = 2097152):
    ''' Replace signal with zero. N for sequence and 0 for CTCF and ATAC '''
    seq[window:] = [0, 0, 0, 0, 1]
    ctcf[window:] = 0
    atac[window:] = 0
    return seq, ctcf, atac

def delete(start, end, seq, ctcf, atac, window = 2097152):
    seq = np.delete(seq, np.s_[start:end], axis = 0)
    ctcf = np.delete(ctcf, np.s_[start:end])
    atac = np.delete(atac, np.s_[start:end])
    return seq[:window], ctcf[:window], atac[:window]

if __name__ == '__main__':
    main()

