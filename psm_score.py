import argparse
import yaml
import random
import itertools as itt
import numpy as np

parser = argparse.ArgumentParser(description='Evaluate PSM score from PeptideShaker reports.')
parser.add_argument('file_paths', metavar='file_paths', type=str, nargs='+',
                    help='space separated list of absolute filepaths to the PeptideShaker reports.')
args = parser.parse_args()


class Spectrum:
    '''
    Class to parse psm annotated spectra files from PeptideShaker
    '''
    def __init__(self, file_path):
        self.file_path = file_path

    def parse_report(self):
        annotated_spectrum = []
        with open(self.file_path, 'r') as f:
            for n,l in enumerate(f):
                ls = l.strip().split('\t')
                if n==0:
                    self.pep_seq = ls[1]
                    continue
                if n==1:
                    keys = ['peak_n',]+ls
                    print(keys)
                    continue
                line = dict(zip(keys, ls))
                line['Intensity'] = float(line['Intensity'])
                annotated_spectrum.append(line)
        return annotated_spectrum


class Score:
    '''
    Class to provide score calculators for parsed spectra
    '''

    def hyperscore(self, peptide, spectrum):
        hyp = dict()
        hyp['Ny'] = sum([1 for peak in spectrum if 'y' in peak['Type']])
        hyp['Nb'] = sum([1 for peak in spectrum if 'b' in peak['Type']])
        hyp['Iy'] = sum([peak['Intensity'] for peak in spectrum if 'y' in peak['Type']])
        hyp['Ib'] = sum([peak['Intensity'] for peak in spectrum if 'b' in peak['Type']])
        for k in hyp.keys():
            if hyp[k]==0:
                hyp[k]=1
        hscore = np.log(np.math.factorial(hyp['Nb'])*np.math.factorial(hyp['Ny'])*hyp['Iy']*hyp['Ib'])

        return hscore

    def mvhscore(self, peptide, spectrum):

        pass

    def get_random_peptides(self, pep_seq, n_random=10000):
        pep_aas = list(pep_seq)
        rand_pep_seqs = set()
        if len(pep_seq) < 8 or len(set(pep_seq)) < 7:
            return list(set([''.join(x) for x in itt.permutations(pep_seq)]))
        while len(rand_pep_seqs)<n_random:
            random.shuffle(pep_aas)
            rand_pep = ''.join(pep_aas)
            if rand_pep != pep_seq:
                rand_pep_seqs.add(rand_pep)
        return rand_pep_seqs

class TheoreticalSpectrum:
    def __init__(self, pep_seq):
        pep_len = len(pep_seq)
        with open('psm_score_calculator/aa_weight.yml', 'r') as f:
            self.aa_weights = yaml.load(f)

        b_ions, y_ions = np.zeros(pep_len), np.zeros(pep_len)
        for i in range(pep_len):
            b_ions[i]           = sum(self.aa_weights[aa] for aa in pep_seq[:i+1])
            y_ions[pep_len-1-i] = sum(self.aa_weights[aa] for aa in pep_seq[i:])

        self.b_ions = b_ions
        self.y_ions = y_ions

if __name__ == "__main__":
    for fpath in list(args.file_paths):
        spectrum = Spectrum(fpath)
        score = Score()
        annotated_spec = spectrum.parse_report()
        hscore = score.hyperscore(spectrum.pep_seq, annotated_spec)
        print('peptide: {}  hyperscore: {}'.format(spectrum.pep_seq, hscore))
