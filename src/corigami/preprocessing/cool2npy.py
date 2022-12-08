#!/usr/bin/env python
import argparse
import numpy as np
from cooler import Cooler
from pathlib import Path

def main(path, save_path, resolution, window_size, balance=True):
    hic = Cooler(f'{path}::resolutions/{resolution}')
    data = hic.matrix(balance=balance, sparse=True)
    # main loop
    for chrom in hic.chromnames:
        mat = data.fetch(chrom)
        diags = compress_diag(mat, window_size)
        ucsc_chrom = f'{chrom}.npz' if chrom.startswith('chr') else f'chr{chrom}.npz'
        chrom_path = save_path / ucsc_chrom
        np.savez(chrom_path, **diags)

def compress_diag(mat, window):
    # NOTE: dict is probably suboptimal here. We could have a big list double the window_size
    diag_dict = {}
    for d in range(window):
        diag_dict[str(d)] = np.nan_to_num(mat.diagonal(d).astype(np.half))
        diag_dict[str(-d)] = np.nan_to_num(mat.diagonal(-d).astype(np.half))
    return diag_dict

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extract chromosome-matrix diagonals from mcool file')
    parser.add_argument('path', help='Path to mcool file')
    parser.add_argument('outdir', help='Directory to save files to. Will be created if need but not its parents')
    parser.add_argument('-r', '--resolution', type=int, default=10000,
                        help='Matrix resolution to use [default: 10000]')
    parser.add_argument('-w', '--window', type=int, default=256,
            help='Number of diagonals to extract [default: 256]')
    parser.add_argument('--no-balance', dest='balance', action='store_false', help='Do not use balanced matrix')
    argv = parser.parse_args()
    outdir = Path(argv.outdir)
    outdir.mkdir(exist_ok=True)
    main(argv.path, outdir, resolution=argv.resolution, window_size=argv.window, balance=argv.balance)
