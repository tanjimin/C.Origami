import pandas as pd
from torch.utils.data import Dataset
from corigami.data.chromosome_dataset import ChromosomeDataset
import corigami.data.data_feature as data_feature

class GenomeDataset(Dataset):
    '''
    Load all chromosomes
    '''
    def __init__(self, celltype_root, 
                       genome_assembly,
                       feat_dicts, 
                       mode = 'train', 
                       include_sequence = True,
                       include_genomic_features = True,
                       use_aug = True):
        self.data_root = celltype_root
        self.include_sequence = include_sequence
        self.include_genomic_features = include_genomic_features
        if not self.include_sequence: print('Not using sequence!')
        if not self.include_genomic_features: print('Not using genomic features!')
        self.use_aug = use_aug

        if mode != 'train': self.use_aug = False # Set augmentation

        # Assign train/val/test chromosomes
        self.chr_names = self.get_chr_names(genome_assembly)
        if mode == 'train':
            self.chr_names.remove('chr10')
            self.chr_names.remove('chr15')
            self.chr_names.remove('chrX') # chrX removed for consistency
        elif mode == 'val':
            self.chr_names = ['chr10']
        elif mode == 'test':
            self.chr_names = ['chr15']
        else:
            raise Exception(f'Unknown mode {mode}')

        # Load genomewide features
        self.genomic_features = self.load_features(f'{celltype_root}/genomic_features', feat_dicts)

        # Load regions to be ignored
        self.centrotelo_dict = self.proc_centrotelo(f'{celltype_root}/../centrotelo.bed')
        # Load chromsome_datasets as part of the dictionary
        self.chr_data, self.lengths = self.load_chrs(self.chr_names, self.genomic_features)
        # Build chrmosome lookup table from the genome
        self.ranges = self.get_ranges(self.lengths)

    def __getitem__(self, idx):
        # Query for chromosome name and where in the chromosome
        chr_name, chr_idx = self.get_chr_idx(idx)
        seq, features, mat, start, end = self.chr_data[chr_name][chr_idx]
        if self.include_sequence:
            if self.include_genomic_features: # Both
                outputs = seq, features, mat, start, end, chr_name, chr_idx
            else: # sequence only  
                outputs = seq, mat, start, end, chr_name, chr_idx
        else: 
            if self.include_genomic_features: # features only
                outputs = features, mat, start, end, chr_name, chr_idx
            else: raise Exception('Must have at least one of the sequence or features')
        return outputs

    def __len__(self):
        return sum(self.lengths)
        
    def load_chrs(self, chr_names, genomic_features):
        '''
        Load chromosome data into a dictionary
        '''
        print('Loading chromosome datasets...')
        chr_data_dict = {}
        lengths = []
        for chr_name in chr_names:
            omit_regions = self.centrotelo_dict[chr_name]
            chr_data_dict[chr_name] = ChromosomeDataset(self.data_root, chr_name, omit_regions, genomic_features, self.use_aug)
            lengths.append(len(chr_data_dict[chr_name]))
        print('Chromosome datasets loaded')
        return chr_data_dict, lengths

    def load_features(self, root_dir, feat_dicts):
        '''
        Args:
            features: a list of dicts with 
                1. file name
                2. norm status
        Returns:
            feature_list: a list of genomic features (bigwig files)
        '''
        feat_list = []
        for feat_item in list(feat_dicts.values()):
            file_name = feat_item['file_name']
            file_path = f'{root_dir}/{file_name}'
            norm = feat_item['norm']
            feat_list.append(data_feature.GenomicFeature(file_path, norm))
        return feat_list
        
    def get_chr_names(self, assembly):
        '''
        Get a list of all chr names. e.g. [chr1 , chr2, ...]
        '''
        print(f'Using Assembly: {assembly}')
        if assembly in ['hg38', 'hg19']:
            chrs = list(range(1, 23))
        elif assembly in ['mm10', 'mm9']:
            chrs = list(range(1, 20))
        else: raise Exception(f'Assembly {assembly} unknown')
        chrs.append('X')
        #chrs.append('Y')
        chr_names = []
        for chr_num in chrs:
            chr_names.append(f'chr{chr_num}')
        return chr_names

    def get_ranges(self, lengths):
        current_start = 0
        ranges = []
        for length in lengths:
            ranges.append([current_start, current_start + length - 1])
            current_start += length
        return ranges

    def get_chr_idx(self, idx):
        '''
        Check index and return chr_name and chr index
        '''
        for i, chr_range in enumerate(self.ranges):
            start, end = chr_range
            if start <= idx <= end:
                return self.chr_names[i], idx - start

    def proc_centrotelo(self, bed_dir):
        ''' Take a bed file indicating location, output a dictionary of items 
        by chromosome which contains a list of 2 value lists (range of loc)
        '''
        df = pd.read_csv(bed_dir , sep = '\t', names = ['chr', 'start', 'end'])
        chrs = df['chr'].unique()
        centrotelo_dict = {}
        for chr_name in chrs:
            sub_df = df[df['chr'] == chr_name]
            regions = sub_df.drop('chr', axis = 1).to_numpy()
            centrotelo_dict[chr_name] = regions
        return centrotelo_dict

if __name__ == '__main__':
    main()
