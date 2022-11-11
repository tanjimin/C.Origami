# C.Origami

[Models](#Download-model-and-other-relevant-resource-files) |
[GitHub](https://github.com/tanjimin/C.Origami) |
[Publications](#list-of-papers)

C.Origami is a deep neural network model enables *de novo* cell type-specific chromatin architecture predictions. C.Origami originates from the Origami architecture that incorporates DNA sequence and cell type-specific features for downstream tasks. It can predict the effect of aberrant genome reorganization such as translocations. In addition, it can be used to perform high-throughput *in silico* genetic screening to identify chromatin related *trans*-acting factors. 

## Dependencies and Installation

### Create a new conda environment for C.Origami and install dependencies
```bash
conda create -n corigami python=3.9
conda activate corigami

pip install torch==1.12.0 torchvision==0.13.0 pandas==1.3.0 matplotlib==3.3.2 pybigwig==0.3.18 omegaconf==2.1.1 tqdm==4.64.0
```

## Download dataset files and pretrained model weights

| File name | Content | 
| - | - | 
| [corigami_data.tar.gz](https://zenodo.org/record/7226561/files/corigami_data.tar.gz?download=1) | DNA reference sequence, CTCF ChIP-seq(IMR-90), ATAC-seq(IMR-90), Hi-C matrix(IMR-90), pretrained model weights | 
| [corigami_data_gm12878_add_on.tar.gz](https://zenodo.org/record/7226561/files/corigami_data.tar.gz?download=1) | CTCF ChIP-seq(GM12878), ATAC-seq(GM12878), Hi-C matrix(GM12878) | 

To run inference or training, you may download [corigami_data.tar.gz](https://zenodo.org/record/7226561/files/corigami_data.tar.gz?download=1) which contains the training data from IMR-90 cell line, and pretrained model weights. 
To test performance on GM12878 *de novo* prediction, you need to additionally download the add-on data file [corigami_data_gm12878_add_on.tar.gz](https://zenodo.org/record/7226561/files/corigami_data.tar.gz?download=1) and unzip it under `corigami_data/data/hg38`.

### Prediction with DNA sequence, CTCF ChIP-seq, and ATAC-seq data 
In order to use our pipeline we require the sequencing data to be pre-processed. Reference DNA sequence is included in the corigami_data file provided above.The input for both the CTCF and ATAC data should be in the form of a bigwig (bw) file. The bigwig should be normalized to the total number of reads. Data quality can be inspected using an applications such as [IGV](https://igv.org).


# Inference

C.Origami can perform de novo prediction of cell type-specific chromatin architecture using both DNA sequence features and cell type-specific genomic information.

For any inference application, download one of our pre-trained models or use your own model. C.Origami is pre-trained on the human IMR-90 cell line (hg38 assembly). Before inference, please download dataset and change path according to the instruction.

Inference allows you to pick between 3 tasks: **predict**, **perturbation**, or **screening**. Examples for each one and the required parameters are under the `examples` folder. 

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
  <img  src="https://github.com/tanjimin/C.Origami/blob/dev/src/corigami/examples/imgs/chr2_500000.png">
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
  <img  src="https://github.com/tanjimin/C.Origami/blob/dev/src/corigami/examples/imgs/chr2_500000_del_1500000_100000_padding_zero.png">
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
  <img  src="https://github.com/tanjimin/C.Origami/blob/dev/src/corigami/examples/imgs/chr2_screen_1250000_2250000_width_1000_step_1000.png">
  </p>


# Training

Coming up soon!

## Citation

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
