# c-origami

[![LICENSE](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/tanjimin/C.Origami-release/blob/dev/LICENSE)

[Models](#Download-model-and-other-relevant-resource-files) |
[GitHub](https://github.com/tanjimin/C.Origami-release) |
[Publications](#list-of-papers)

C-Origami is a deep neural network model for predicting de novo cell type-specific chromatin architecture. By incorporating DNA sequence, CTCF binding, and chromatin accessibility profiles, C-Origami achieves accurate cell type-specific prediction.

Publications associated with the C. Origami project can be found
[at the end of this README](#list-of-papers).


## Documentation

### CTCF/ATAC/DNA data 
In order to use our pipeline we require the sequencing data to be pre-processed. The input for both the CTCF and ATAC data should be in the form of a bigwig (bw) file. The bigwig should be normalized to the total number of reads. Data quality can be inspected using an applications such as [IGV](https://igv.org).
C.Origami has been trained on the human IMR-90 cell line (hg38 assembly). You can download the training data here: 

TODO

For human chromosomes:
```bash
download.sh hg38
```
For mouse chromosomes:
```bash
download.sh mm10
```
## Dependencies and Installation

### Create and activate a new virtual environment
```bash
conda create -n corigami python=3.6
conda activate corigami
```
### Install the package and other requirements
```bash
pip install torch torchvision pybigwig pytorch_lightning pandas matplotlib
```

### Installing Directly from the Github

If you want to install directly from the GitHub source, clone the repository:
```
git clone https://github.com/tanjimin/C.Origami-release.git
```
## Download model and other relevant resource files

If you wish to use our pretrained model, you may download the model and other files needed for running C.Origami by running the command below. You can also train your own model instead of using ours. See the [training](#Training) section below.

Our model is trained on IMR-90 cell line. The training data can be accessed here:
```bash
wget -O data https://www.dropbox.com/s/oor01snnekyh4s5/imr90_data.tar.gz?dl=0
```

# Training

TODO

# Inference

C.Origami can perform de novo prediction of cell type-specific chromatin architecture using both DNA sequence features and cell type-specific genomic information.

For any inference application, download one of our pre-trained models or use your own model. Inference allows you to pick between 3 tasks: **predict**, **perturbation**, or **screening**. Examples for each one and the required parameters are under the `examples` folder. 

## Prediction

Prediction will produce both an image of the 2MB window as well as a numpy matrix for further downstream analysis. 

```docs

    Usage:
    python corigami/inference/prediction.py [options] 

    Options:
    -h --help       Show this screen.
    --out           Output path for storing results
    --celltype      Sample cell type for prediction, used for output separation
    --chr           Chromosome for prediction')
    --start         Starting point for prediction (width defaults to 2097152 bp which is the input window size)
    --model         Path to the model checkpoint')
    --seq           Path to the folder where the sequence .fa.gz files are stored
    --ctcf          Path to the folder where the CTCF ChIP-seq .bw files are stored
    --atac          Path to the folder where the ATAC-seq .bw files are stored
```

`An example of a C.Origami predicted (2MB window) Hi-C matrix for the IMR-90 cell line at chromosome 2 with start position 500,000:`
<p align="center">
  <img  src="https://github.com/tanjimin/C.Origami-release/blob/dev/src/corigami/examples/imgs/chr2_500000.png">
  </p>

## Editing/Perturbation

For now the only perturbation implemented is deletion. Specify the same parameters as before along with specific deletion parameters. If you want to do multiple deletions, you can specify in the config by creating additional start and end positions. 

```docs
  Usage:
    python corigami/inference/editing.py [options] 

    Options:
    -h --help       Show this screen.
    --out           Output path for storing results
    --celltype      Sample cell type for prediction, used for output separation
    --chr           Chromosome for prediction')
    --start         Starting point for prediction (width defaults to 2097152 bp which is the input window size)
    --model         Path to the model checkpoint')
    --seq           Path to the folder where the sequence .fa.gz files are stored
    --ctcf          Path to the folder where the CTCF ChIP-seq .bw files are stored
    --atac          Path to the folder where the ATAC-seq .bw files are stored

    --del-start     Starting point for deletion
    --del-width     Width for deletion
    --padding       Padding type, either zero or follow. Using zero: the missing region at the end will be padded with zero for ctcf and atac seq, while sequence will be padded with N (unknown necleotide). Using follow: the end will be padded with features in the following region
    --hide-line     Remove the line showing deletion site
```

`An example of a C.Origami predicted (2MB window) Hi-C matrix for the IMR-90 cell line at chromosome 2 with start position 500,000 and a deletion from 1.5MB to 1.6MB (100,000 basepairs deleted):`
<p align="center">
  <img  src="https://github.com/tanjimin/C.Origami-release/blob/dev/src/corigami/examples/imgs/chr2_500000_del_1500000_100000_padding_zero.png">
  </p>

## Screening

In silico genetic screening can be used to see what regions of perturbation lead to the greatest impact on the prediction. Running this task will result in a bedgraph file consisting of the chr number, start position, end position, and impact score. The more impact the perturbation had, the higher the impact score.

Screening can be done only for one chromosome at a time. The end position unless otherwise specified will be 2MB from the start position specified above it. The `perturb-width` is allows you to set the size of the deletion you want to make or in other words how many base pairs to remove. The `step-size` is how far each deletion is from the past deletion (start position) - please note it is fine for the deletions to overlap. 

```docs
  Usage:
    python corigami/inference/screening.py [options] 

    Options:
    -h --help       Show this screen.
    --out           Output path for storing results
    --celltype      Sample cell type for prediction, used for output separation
    --chr           Chromosome for prediction')
    --start         Starting point for prediction (width defaults to 2097152 bp which is the input window size)
    --model         Path to the model checkpoint')
    --seq           Path to the folder where the sequence .fa.gz files are stored
    --ctcf          Path to the folder where the CTCF ChIP-seq .bw files are stored
    --atac          Path to the folder where the ATAC-seq .bw files are stored
    
    --screen-start        Starting point for screening
    --screen-end          Ending point for screening
    --perturb-width       Width of perturbation used for screening
    --step-size           step size of perturbations in screening
    --plot-impact-score   Plot impact score and save png. (Not recommended for large scale screening (>10000 perturbations)')
    --save-pred           Save prediction tensor
    --save-perturbation   Save perturbed tensor
    --save-diff           Save difference tensor
    --save-impact-score   Save impact score array
    --save-bedgraph       Save bedgraph file for impact score
    --save-frames         Save each deletion instance with png and npy (Could be taxing on computation and screening, not recommended).')
    --padding             Padding type, either zero or follow. Using zero: the missing region at the end will be padded with zero for ctcf and atac seq, while sequence will be padded with N (unknown necleotide). Using follow: the end will be padded with features in the following region
```
**Please note that screening can be very computationally intensive especially when screening at a 1 Kb resolution or less. For instance, screening on chromosome 8, a medium-size chromosome which has a length of 146Mb, requires the model to make 146Mb / 1Kb * 2 predictions = 292,000 separate predictions.**

`An example of a barplot representing the impact score of each perturbation. C.Origami screened chromosome 2 from position 1.25 MB to 2.25 MB with a perturbation of 1000 basepairs (perturb-width) being made every 1000 basepairs (step-size):`
<p align="center">
  <img  src="https://github.com/tanjimin/C.Origami-release/blob/dev/src/corigami/examples/imgs/chr2_screen_1250000_2250000_width_1000_step_1000.png">
  </p>

## License

C-Origami is MIT licensed, as found in the [LICENSE file](https://github.com/tanjimin/C.Origami-release/blob/dev/LICENSE).

## Cite

If you use the C-Origami code in your project, please cite the bioRxiv paper:

```BibTeX
@inproceedings{tan2020,
    title={Cell type-specific prediction of 3D chromatin architecture},
    author={Jimin Tan and Javier Rodriguez-Hernaez and Theodore Sakellaropoulos and Francesco Boccalatte and Iannis Aifantis and Jane Skok and David Fenyö and Bo Xia and Aristotelis Tsirigos},
    journal = {bioRxiv e-prints},
    archivePrefix = "bioRxiv",
    doi = {10.1101/2022.03.05.483136},
    year={2022}
}
```


## List of Papers

The following lists titles of papers from the C-Origami project. 

Cell type-specific prediction of 3D chromatin architecture
Jimin Tan, Javier Rodriguez-Hernaez, Theodore Sakellaropoulos, Francesco Boccalatte, Iannis Aifantis, Jane Skok, David Fenyö, Bo Xia, Aristotelis Tsirigos
bioRxiv 2022.03.05.483136; doi: https://doi.org/10.1101/2022.03.05.483136
