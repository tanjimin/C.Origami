import gzip
import numpy as np
import pyBigWig as pbw

class Feature():

    def __init__(self, **kwargs):
        self.load(**kwargs)
    
    def load(self):
        raise Exception('Not implemented')

    def get(self):
        raise Exception('Not implemented')

    def __len__(self):
        raise Exception('Not implemented')

class HiCFeature(Feature):

    def load(self, path = None):
        self.hic = self.load_hic(path)

    def get(self, start, window = 2097152, res = 10000):
        start_bin = int(start / res)
        range_bin = int(window / res)
        end_bin = start_bin + range_bin
        hic_mat = self.diag_to_mat(self.hic, start_bin, end_bin)
        return hic_mat

    def load_hic(self, path):
        print(f'Reading Hi-C: {path}')
        return dict(np.load(path))

    def diag_to_mat(self, ori_load, start, end):
        '''
        Only accessing 256 x 256 region max, two loops are okay
        '''
        square_len = end - start
        diag_load = {}
        for diag_i in range(square_len):
            diag_load[str(diag_i)] = ori_load[str(diag_i)][start : start + square_len - diag_i]
            diag_load[str(-diag_i)] = ori_load[str(-diag_i)][start : start + square_len - diag_i]
        start -= start
        end -= start

        diag_region = []
        for diag_i in range(square_len):
            diag_line = []
            for line_i in range(-1 * diag_i, -1 * diag_i + square_len):
                if line_i < 0:
                    diag_line.append(diag_load[str(line_i)][start + line_i + diag_i])
                else:
                    diag_line.append(diag_load[str(line_i)][start + diag_i])
            diag_region.append(diag_line)
        diag_region = np.array(diag_region).reshape(square_len, square_len)
        return diag_region

    def __len__(self):
        return len(self.hic['0'])

class GenomicFeatureSingleThread(Feature):

    def __init__(self, path, norm):
        self.path = path
        self.load(path)
        self.norm = norm
        print(f'Feature path: {path} \n Normalization status: {norm}')

    def load(self, path):
        self.feature = self.read_feature(path)

    def get(self, chr_name, start, end):
        feature = self.feature_to_npy(chr_name, start, end)
        feature = np.nan_to_num(feature, 0) # Important! replace nan with 0
        if self.norm == 'log':
            feature = np.log(feature + 1)
        elif self.norm is None:
            feature = feature
        else:
            raise Exception(f'Norm type {self.norm} undefined')
        return feature

    def read_feature(self, path):
        '''
        read bigwig file
        '''
        bw_file = pbw.open(path)
        return bw_file

    def feature_to_npy(self, chr_name, start, end):
        signals = self.feature.values(chr_name, start, end)
        return np.array(signals)

    def length(self, chr_name):
        return self.feature.chroms(chr_name)

class GenomicFeature(GenomicFeatureSingleThread):

    def __init__(self, path, norm):
        self.path = path
        self.norm = norm
        print(f'Feature path: {path} \n Normalization status: {norm}')

    def load(self, path):
        raise Exception('Left blank')

    def feature_to_npy(self, chr_name, start, end):
        with pbw.open(self.path) as bw_file:
            signals = bw_file.values(chr_name, int(start), int(end))
        return np.array(signals)

    def length(self, chr_name):
        with pbw.open(self.path) as bw_file:
            length = bw_file.chroms(chr_name)
        return length

class SequenceFeature(Feature):

    def load(self, path = None):
        self.seq = self.read_seq(path)

    def get(self, start, end):
        seq = self.seq_to_npy(self.seq, start, end)
        onehot_seq = self.onehot_encode(seq)
        return onehot_seq

    def __len__(self):
        return len(self.seq)

    def read_seq(self, dna_dir):
        '''
        Transform fasta data to numpy array
        
        Args:
            dna_dir (str): Directory to DNA .fa path

        Returns:
            array: A numpy char array that contains DNA for a chromosome
        '''
        print(f'Reading sequence: {dna_dir}')
        with gzip.open(dna_dir, 'r') as f:
            seq = f.read().decode("utf-8")
        seq = seq[seq.find('\n'):]
        seq = seq.replace('\n', '').lower()
        return seq
        
    def seq_to_npy(self, seq, start, end):
        '''
        Transform fasta data to integer numpy array
        
        Args:
            dna_dir (str): Directory to DNA .fa path

        Returns:
            array: A numpy char array that contains DNA for a chromosome
        '''
        seq = seq[start : end]
        en_dict = {'a' : 0, 't' : 1, 'c' : 2, 'g' : 3, 'n' : 4}
        en_seq = [en_dict[ch] for ch in seq]
        np_seq = np.array(en_seq, dtype = int)
        return np_seq

    def onehot_encode(self, seq):
        ''' 
        encode integer dna array to onehot (n x 5)
        Args:
            seq (arr): Numpy array (n x 1) of dna encoded as 0-4 integers

        Returns:
            array: A numpy matrix (n x 5)
        '''
        seq_emb = np.zeros((len(seq), 5))
        seq_emb[np.arange(len(seq)), seq] = 1
        return seq_emb

if __name__ == '__main__':
    main()
