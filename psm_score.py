import argparse
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
                    keys = ls
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


if __name__ == "__main__":
    for fpath in list(args.file_paths):
        spectrum = Spectrum(fpath)
        score = Score()
        annotated_spec = spectrum.parse_report()
        hscore = score.hyperscore(spectrum.pep_seq, annotated_spec)
        print('peptide: {}  hyperscore: {}'.format(spectrum.pep_seq, hscore))
