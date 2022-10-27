import corigami.inference.utils.inference_utils as infer
from corigami.inference.utils import plot_utils 
import argparse

def main():
    parser = argparse.ArgumentParser(description='C.Origami Prediction Module.')
    parser.add_argument('--out', dest='output_path', 
                        default='./outputs',
                        help='output path for storing results')
    parser.add_argument('--celltype', dest='celltype', 
                        help='Sample cell type for prediction, used for output separation')
    parser.add_argument('--chr', dest='chr_name', 
                        help='Chromosome for prediction')
    parser.add_argument('--start', dest='start', type=int,
                        help='Starting point for prediction (width defaults to 2097152 bp which is the input window size)')
    parser.add_argument('--model', dest='model_path', 
                        help='Path to the model checkpoint')
    parser.add_argument('--seq', dest='seq_path', 
                        help='Path to the folder where the sequence .fa.gz files are stored')
    parser.add_argument('--ctcf', dest='ctcf_path', 
                        help='Path to the folder where the CTCF ChIP-seq .bw files are stored')
    parser.add_argument('--atac', dest='atac_path', 
                        help='Path to the folder where the ATAC-seq .bw files are stored')

    args = parser.parse_args()
    single_prediction(args.output_path, args.celltype, 
                      args.chr_name, args.start,
                      args.model_path, 
                      args.seq_path, args.ctcf_path, args.atac_path)

def single_prediction(output_path, celltype, chr_name, start, model_path, seq_path, ctcf_path, atac_path):
    seq_region, ctcf_region, atac_region = infer.load_region(chr_name, 
            start, seq_path, ctcf_path, atac_path)
    pred = infer.prediction(seq_region, ctcf_region, atac_region, 
                                   model_path)
    plot = plot_utils.MatrixPlot(output_path, pred, 'prediction', celltype, 
                                 chr_name, start)
    plot.plot()

if __name__ == '__main__':
    main()
