import os
import numpy as np
import pandas as pd
import sys
import torch

import corigami.inference.utils.inference_utils as infer
from corigami.inference import editing 
from corigami.inference.utils import plot_utils, model_utils

import argparse

def main():
    parser = argparse.ArgumentParser(description='C.Origami Screening Module.')
    
    # Output location
    parser.add_argument('--out', dest='output_path', 
                        default='outputs',
                        help='output path for storing results (default: %(default)s)')

    # Location related params
    parser.add_argument('--celltype', dest='celltype', 
                        help='Sample cell type for prediction, used for output separation', required=True)
    parser.add_argument('--chr', dest='chr_name', 
                        help='Chromosome for prediction', required=True)
    parser.add_argument('--model', dest='model_path', 
                        help='Path to the model checkpoint', required=True)
    parser.add_argument('--seq', dest='seq_path', 
                        help='Path to the folder where the sequence .fa.gz files are stored', required=True)
    parser.add_argument('--ctcf', dest='ctcf_path', 
                        help='Path to the folder where the CTCF ChIP-seq .bw files are stored', required=True)
    parser.add_argument('--atac', dest='atac_path', 
                        help='Path to the folder where the ATAC-seq .bw files are stored', required=True)

    # Screening related params
    parser.add_argument('--screen-start', dest='screen_start', type=int,
                        help='Starting point for screening.', required=True)
    parser.add_argument('--screen-end', dest='screen_end', type=int,
                        help='Ending point for screening.', required=True)
    parser.add_argument('--perturb-width', dest='perturb_width', type=int,
                        help='Width of perturbation used for screening.', required=True)
    parser.add_argument('--step-size', dest='step_size', type=int,
                        help='step size of perturbations in screening.', required=True)

    # Saving related params
    parser.add_argument('--plot-impact-score', dest='plot_impact_score', 
                        action = 'store_true',
                        help='Plot impact score and save png. (Not recommended for large scale screening (>10000 perturbations)')
    parser.add_argument('--save-pred', dest='save_pred', 
                        action = 'store_true',
                        help='Save prediction tensor (default: %(default)s)')
    parser.add_argument('--save-perturbation', dest='save_perturbation', 
                        action = 'store_true',
                        help='Save perturbed tensor (default: %(default)s)')
    parser.add_argument('--save-diff', dest='save_diff', 
                        action = 'store_true',
                        help='Save difference tensor (default: %(default)s)')
    parser.add_argument('--save-impact-score', dest='save_impact_score', 
                        action = 'store_true',
                        help='Save impact score array (default: %(default)s)')
    parser.add_argument('--save-bedgraph', dest='save_bedgraph', 
                        action = 'store_true',
                        help='Save bedgraph file for impact score (default: %(default)s)')
    parser.add_argument('--save-frames', dest='plot_frames', 
                        action = 'store_true',
                        help='Save each deletion instance with png and npy (Could be taxing on computation and screening, not recommended) (default: %(default)s)')

    # Perturbation related params
    parser.add_argument('--padding', dest='end_padding_type', 
                        default='zero',
                        help='Padding type, either zero or follow. Using zero: the missing region at the end will be padded with zero for ctcf and atac seq, while sequence will be padded with N (unknown necleotide). Using follow: the end will be padded with features in the following region (default: %(default)s)')

    args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])
    screening(args.output_path, args.celltype, args.chr_name, args.screen_start, 
              args.screen_end, args.perturb_width, args.step_size,
              args.model_path,
              args.seq_path, args.ctcf_path, args.atac_path, 
              args.save_pred, args.save_perturbation,
              args.save_diff, args.save_impact_score,
              args.save_bedgraph, args.plot_impact_score, args.plot_frames)

def screening(output_path, celltype, chr_name, screen_start, screen_end, perturb_width, step_size, model_path, seq_path, ctcf_path, atac_path, save_pred = False, save_deletion = False, save_diff = True, save_impact_score = True, save_bedgraph = True, plot_impact_score = True, plot_frames = False):
    # Store data and model in memory
    seq, ctcf, atac = infer.load_data_default(chr_name, seq_path, ctcf_path, atac_path)
    model = model_utils.load_default(model_path)
    # Generate pertubation windows
    # Windows are centered. Thus, both sides have enough margins
    windows = [w * step_size + screen_start for w in range(int((screen_end - screen_start) / step_size))]
    from tqdm import tqdm
    preds = np.empty((0, 256, 256))
    preds_deletion = np.empty((0, 256, 256))
    diff_maps = np.empty((0, 256, 256))
    perturb_starts = []
    perturb_ends = []
    print('Screening...')
    for w_start in tqdm(windows):
        pred_start = int(w_start + perturb_width / 2 - 2097152 / 2)
        pred, pred_deletion, diff_map = predict_difference(chr_name, pred_start, int(w_start), perturb_width, model, seq, ctcf, atac)
        if plot_frames:
            plot_combination(output_path, celltype, chr_name, pred_start, w_start, perturb_width, pred, pred_deletion, diff_map, 'screening')
        preds = np.append(preds, np.expand_dims(pred, 0), axis = 0)
        preds_deletion = np.append(preds_deletion, np.expand_dims(pred_deletion, 0), axis = 0)
        diff_maps = np.append(diff_maps, np.expand_dims(diff_map, 0), axis = 0)
        perturb_starts.append(w_start)
        perturb_ends.append(w_start + perturb_width)
    impact_scores = np.abs(diff_maps).mean(axis = (1, 2))
    plot = plot_utils.MatrixPlotScreen(output_path, perturb_starts, perturb_ends, impact_scores, diff_maps, preds, preds_deletion, 'screening', celltype, chr_name, screen_start, screen_end, perturb_width, step_size, plot_impact_score)
    figure = plot.plot()
    plot.save_data(figure, save_pred, save_deletion, save_diff, save_impact_score, save_bedgraph)

def predict_difference(chr_name, start, deletion_start, deletion_width, model, seq, ctcf, atac):
    # Define window which accomodates deletion
    end = start + 2097152 + deletion_width
    seq_region, ctcf_region, atac_region = infer.get_data_at_interval(chr_name, start, end, seq, ctcf, atac)
    # Unmodified inputs
    inputs = preprocess_prediction(chr_name, start, seq_region, ctcf_region, atac_region)
    pred = model(inputs)[0].detach().cpu().numpy() # Prediction
    # Inputs with deletion
    inputs_deletion = preprocess_deletion(chr_name, start, deletion_start, 
            deletion_width, seq_region, ctcf_region, atac_region) # Get data
    pred_deletion = model(inputs_deletion)[0].detach().cpu().numpy() # Prediction
    # Compare inputs:
    diff_map = pred_deletion - pred
    return pred, pred_deletion, diff_map

def plot_combination(output_path, celltype, chr_name, start, deletion_start, deletion_width, pred, pred_deletion, diff_map, plot_type = 'point_screening'):
    # Initialize plotting class
    # Plot predicted map
    plot = plot_utils.MatrixPlot(output_path, pred, plot_type, celltype, 
                                 chr_name, start)
    plot.plot()
    # Plot deleted map
    plot = plot_utils.MatrixPlotDeletion(output_path, pred_deletion, 
            plot_type, celltype, chr_name, start, deletion_start, 
            deletion_width, padding_type = 'zero', show_deletion_line = True)
    plot.plot()
    # Plot difference map
    plot = plot_utils.MatrixPlotPointScreen(output_path, diff_map, 
            plot_type, celltype, chr_name, start, deletion_start, 
            deletion_width, padding_type = 'zero', show_deletion_line = False)
    plot.plot()

def preprocess_prediction(chr_name, start, seq_region, ctcf_region, atac_region):
    # Delete inputs
    seq_region, ctcf_region, atac_region = trim(seq_region, ctcf_region, 
                                                atac_region)
    # Process inputs
    inputs = infer.preprocess_default(seq_region, ctcf_region, atac_region)
    return inputs

def preprocess_deletion(chr_name, start, deletion_start, deletion_width, seq_region, ctcf_region, atac_region):
    # Delete inputs
    seq_region, ctcf_region, atac_region = editing.deletion_with_padding(start, 
            deletion_start, deletion_width, seq_region, ctcf_region, 
            atac_region, end_padding_type = 'zero')
    # Process inputs
    inputs = infer.preprocess_default(seq_region, ctcf_region, atac_region)
    return inputs

def trim(seq, ctcf, atac, window = 2097152):
    return seq[:window], ctcf[:window], atac[:window]

if __name__ == '__main__':
    main()

