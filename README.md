# C.Origami

[Models](#download-dataset-files-and-pretrained-model-weights) |
[GitHub](https://github.com/tanjimin/C.Origami) |
[Publications](#list-of-papers)

C.Origami is a deep neural network model enables *de novo* cell type-specific chromatin architecture predictions. C.Origami originates from the Origami architecture that incorporates DNA sequence and cell type-specific features for downstream tasks. It can predict the effect of aberrant genome reorganization such as translocations. In addition, it can be used to perform high-throughput *in silico* genetic screening to identify chromatin related *trans*-acting factors. 

## Dependencies and Installation

Create a new conda environment for C.Origami

```bash
conda create -n corigami python=3.9
conda activate corigami
```

First install PyTorch according to the instructions on the
[PyTorch Website](https://pytorch.org/get-started/) for your operating system
and CUDA setup.  

**For inference ONLY dependency use:**

```bash
pip install corigami
```

**For full dependency used for training use:**

```bash
 pip install corigami[training]
```

`pip` will handle all package dependencies. 

### Installing Directly from Source

If you want to install directly from the GitHub, git clone the repository and install all dependencies:

```bash
pip install torch==1.12.0 torchvision==0.13.0 pandas==1.3.0 matplotlib==3.3.2 pybigwig==0.3.18 omegaconf==2.1.1 tqdm==4.64.0
```
Then run:

```bash
pip install -e .
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
    corigami-predict [options] 

    Options:
    -h --help       Show this screen.
    --out           Output path for storing results
    --celltype      Sample cell type for prediction, used for output separation
    --chr           Chromosome for prediction
    --start         Starting point for prediction (width defaults to 2097152 bp which is the input window size)
    --model         Path to the model checkpoint
    --seq           Path to the folder where the sequence .fa.gz files are stored
    --ctcf          Path to the folder where the CTCF ChIP-seq .bw files are stored
    --atac          Path to the folder where the ATAC-seq .bw files are stored
```

`An example of a C.Origami predicted (2MB window) Hi-C matrix for the IMR-90 cell line at chromosome 2 with start position 500,000:`
<p align="center">
  <img  src="https://github.com/tanjimin/C.Origami/blob/dev/examples/imgs/chr2_500000.png">
  </p>

## Editing/Perturbation

For now the only perturbation implemented is deletion. Specify the same parameters as before along with specific deletion parameters. If you want to do multiple deletions, you can specify in the config by creating additional start and end positions. 

```docs
  
  Usage:
  corigami-edit [options] 

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
  <img  src="https://github.com/tanjimin/C.Origami/blob/dev/examples/imgs/chr2_500000_del_1500000_100000_padding_zero.png">
  </p>

## Screening

In silico genetic screening can be used to see what regions of perturbation lead to the greatest impact on the prediction. Running this task will result in a bedgraph file consisting of the chr number, start position, end position, and impact score. The more impact the perturbation had, the higher the impact score.

Screening can be done only for one chromosome at a time. The end position unless otherwise specified will be 2MB from the start position specified above it. The `perturb-width` is allows you to set the size of the deletion you want to make or in other words how many base pairs to remove. The `step-size` is how far each deletion is from the past deletion (start position) - please note it is fine for the deletions to overlap. 

```docs

  Usage:
    corigami-screen [options] 

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
  <img  src="https://github.com/tanjimin/C.Origami/blob/dev/examples/imgs/chr2_screen_1250000_2250000_width_1000_step_1000.png">
  </p>


# Training

You may train your own model on another human or cell mouse line. 

### Genomic Features

You will need a bigwig file of the corresponding **atac** and **ctcf chip** sequence peaks. We recommend using [Seq-N-Slide](https://igordot.github.io/sns/) pipeline for processing raw fastqs into bigwigs using the `atac` or `chip` route. 

### Hi-C data 

Experimental Hi-C matrices are needed for training. We recommend using [HiC-Pro](https://github.com/nservant/HiC-Pro) to process your experimental HiC data. C.origami accepts npz files for each chromosome - therefore, we have included a script under `src/corigami/preprocessing` to convert from mcool (output of HiC-Pro) to npz.

### Data directory

C.origami expects the input data to be structured in the following way:
```docs
root
└── hg38
    ├── centrotelo.bed
    ├── dna_sequence
    │   ├── chr10.fa.gz
    │   ├── chr11.fa.gz
    │   ├── chr12.fa.gz
    │   ├── chr13.fa.gz
    │   ├── chr14.fa.gz
    │   ├── chr15.fa.gz
    │   ├── chr16.fa.gz
    │   ├── chr17.fa.gz
    │   ├── chr18.fa.gz
    │   ├── chr19.fa.gz
    │   ├── chr1.fa.gz
    │   ├── chr20.fa.gz
    │   ├── chr21.fa.gz
    │   ├── chr22.fa.gz
    │   ├── chr2.fa.gz
    │   ├── chr3.fa.gz
    │   ├── chr4.fa.gz
    │   ├── chr5.fa.gz
    │   ├── chr6.fa.gz
    │   ├── chr7.fa.gz
    │   ├── chr8.fa.gz
    │   ├── chr9.fa.gz
    │   ├── chrX.fa.gz
    │   └── chrY.fa.gz
    └── IMR-90
        ├── genomic_features
        │   ├── atac.bw
        │   └── ctcf_log2fc.bw
        └── hic_matrix
            ├── chr10.npy
            ├── chr11.npy
            ├── chr12.npy
            ├── chr13.npy
            ├── chr14.npy
            ├── chr15.npy
            ├── chr16.npy
            ├── chr17.npy
            ├── chr18.npy
            ├── chr19.npy
            ├── chr1.npy
            ├── chr20.npy
            ├── chr21.npy
            ├── chr22.npy
            ├── chr2.npy
            ├── chr3.npy
            ├── chr4.npy
            ├── chr5.npy
            ├── chr6.npy
            ├── chr7.npy
            ├── chr8.npy
            ├── chr9.npy
            └── chrX.npy
```
**Note**: if you choose to download the data from [link above](https://zenodo.org/record/7226561/files/corigami_data.tar.gz?download=1) the data directory will automatically be structured in this way. Then when training set your `--data-root` option to the root directory as shown in the tree above. 

**Note**: if you wish to use another assembly (e.g. mm10) please make sure your data directory is structured as above with the assembly name --> cell type, centrotelo.bed (this is a bed file of any regions you wish to exclude ex. telomeres and centromeres), dna sequence directory. Under the each cell type you should have a folder called `genomic_features`
containing the atac and ctcf bigwigs (make sure to name your files the **exact** same!) and a `hic_matrix` containing a npz file per chromosome. There can be multiple cell types (and thus multiple atac/ctcf/hic files) but only one copy of the dna sequence and centroleo.bed is needed per assembly. 

  Usage:
    corigami-train [options]
    
    Options:
    -h --help       Show this screen.
    --seed          Random seed for training (defaults to 2077)
    --save_path     Path to the model checkpoint
    --data-root     Root path of training data
    --assembly      Genome assembly for training data
    --celltype      Sample cell type for prediction, used for output separation
    --model-type    Type of model architecture (defaults to CNN with Transformer)
    --patience      Epoches before early stopping
    --max-epochs    Max epochs
    --save-top-n    Top n models to save
    --num-gpu       Number of GPUs to use
    --batch-size    Batch size
    --ddp-disabled  Using ddp, adjust batch size
    --num-workers   Dataloader workers
```

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


[Models](#Download-model-and-other-relevant-resource-files) |
[GitHub](https://github.com/tanjimin/C.Origami) |
[Publications](#list-of-papers)

