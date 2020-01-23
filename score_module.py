import os
import argparse
import yaml
import random
from collections import Counter
import itertools as itt
import multiprocessing as mp
import numpy as np


class PSMReport:
    def __init__(self, report_path):
        self.report_path = os.path.abspath(report_path)
        psms = {}
        with open(self.report_path, 'r') as f:
            for n,l in enumerate(f):
                ls = l.split('\t')
                if n==0:
                    psm_keys = [f.strip() for f in ls[1:]]
                    continue
                psm_n = ls[0].split('.')[0]
                if psm_n == '':
                    ann_keys = [f.strip() for f in ls[1:]]
                    continue
                elif psm_n != '' and psm_n not in psms.keys():
                    psms[psm_n] = dict(zip(psm_keys, ls[1:]))
                    psms[psm_n]['Spectrum'] = np.array([
                        (float(x.split(',')[0]), float(x.split(',')[1])) \
                                for x in psms[psm_n]['Spectrum Array List'][2:-2].split('],[')
                    ], dtype = [('m/z', np.float64), ('Intensity', np.float64)])
                    psms[psm_n]['Annotation'] = []
                    continue
                if len(ls[0].split('.')) > 1:
                    ion              = dict(zip(ann_keys, ls[1:]))
                    ion['Subtype']   = bytes(ion['Subtype'], 'utf-8')
                    ion['m/z']       = float(ion['m/z'])
                    ion['Intensity'] = float(ion['Intensity'])
                    psms[psm_n]['Annotation'].append(ion)
        self.psms = psms

    def build_annotation_matrix(self, psm_n):
        ann_dicts = self.psms[psm_n]['Annotation']
        ann_mat = np.empty((len(ann_dicts), 2), dtype=object)
        for n,ann in enumerate(ann_dicts):
            ann_mat[n] = [ann[x] for x in ['Subtype', 'Intensity']]
        return ann_mat


class Score:
    '''
    Class to provide score calculators for parsed spectra
    '''
    def __init__(self, n_random=10000, tol=0.05, n_cpu=1):
        self.n_random = n_random
        self.tol      = tol
        self.n_cpu    = n_cpu

    def hyperscore(self, spec_ann, pep_seq=None, spectrum=None, compute_pval=False):
        if not spec_ann.any():
            return 0
        y    = spec_ann[:,0] == b'y'
        b    = spec_ann[:,0] == b'b'
        n_y  = max(1, sum(y))
        n_b  = max(1, sum(b))
        I_y  = max(1, sum(spec_ann[y,-1]))
        I_b  = max(1, sum(spec_ann[b,-1]))

        try:
            hscore = np.log(np.math.factorial(n_b)*np.math.factorial(n_y)*I_y*I_b)
        except:
            hscore = 1e8

        if compute_pval:
            if spectrum is None or pep_seq is None:
                return 'Spectrum and peptide sequence must be provided for P-value calulation.'
            pval = self.hyperscore_pval(pep_seq, spectrum, hscore)
            hscore = (hscore, pval)

        return hscore

    def hyperscore_pval(self, pep_seq, spectrum, hscore):
        rand_peps = self.get_random_peptides(pep_seq)
        T = TheoreticalSpectrum()

        if self.n_cpu > 1:
            pool = mp.Pool(self.n_cpu)
            args = list(zip(rand_peps, [spectrum,]*len(rand_peps), [T,]*len(rand_peps)))
            result = pool.map_async(self.align_and_score, args)
            hscores = result.get()
            pool.close()
            pool.terminate()
            pool.join()
        else:
            hscores = []
            for pep in rand_peps:
                T.compute_spectrum(pep)
                spec_ann  = self.get_peaks(T.ions, spectrum)
                hscores.append(self.hyperscore(spec_ann))
        hscores = np.array(hscores)
        pval = (1+sum(hscores > hscore))/len(hscores)
        return pval

    def align_and_score(self, pep_spec):
        pep, spec, T = pep_spec
        T.compute_spectrum(pep)
        spec_ann  = self.get_peaks(T.ions, spec)
        return self.hyperscore(spec_ann)

    def get_peaks(self, ions, spectrum):
        aligned_peaks = np.searchsorted(ions['m/z'], spectrum['m/z'])
        in_tol = np.abs(spectrum['m/z'] - ions['m/z'][aligned_peaks-1]) < self.tol
        if sum(in_tol)==0:
            return np.array([0])
        peaks = np.array([tuple(x) for x in ions[aligned_peaks-1][in_tol]], dtype=object)
        spec_int = spectrum['Intensity'][in_tol][:, np.newaxis]
        return np.concatenate((peaks, spec_int), axis=1)

    def get_random_peptides(self, pep_seq):
        pep_aas = list(pep_seq)
        rand_pep_seqs = set()
        n_permute = np.math.factorial(len(pep_seq)) / np.prod([np.math.factorial(n) for n in Counter(pep_seq).values()])
        if n_permute <= self.n_random:
            return list(permutations(pep_seq))
        while len(rand_pep_seqs)<self.n_random:
            random.shuffle(pep_aas)
            rand_pep = ''.join(pep_aas)
            if rand_pep != pep_seq:
                rand_pep_seqs.add(rand_pep)
        return list(rand_pep_seqs)


class TheoreticalSpectrum:
    def __init__(self, charges=[2,3,4]):
        dirname = os.path.dirname(__file__)
        aa_weights = os.path.join(dirname,'aa_weight.yml')
        self.ion_types  = ['y', 'b']
        self.frag_shift = {'y':19., 'b':1.}
        self.charges = charges
        with open(aa_weights, 'r') as f:
            self.aa_weights = yaml.full_load(f)

    def compute_spectrum(self, pep_seq):
        pep_len = len(pep_seq)
        frag_ions = []
        for i in range(pep_len):
            frags      = {'y':pep_seq[i:] ,'b':pep_seq[:i+1]}
            for ion_type in self.ion_types:
                frag = frags[ion_type]

                n_NH3 = sum(frag.count(x) for x in 'RKNQ')
                NH3_shifts = ['']
                if n_NH3:
                    NH3_shifts.extend([''.join(['-NH3']*(i+1)) for i in range(n_NH3)])

                n_H2O = sum(frag.count(x) for x in 'STED')
                H2O_shifts = ['']
                if n_H2O:
                    H2O_shifts.extend([''.join(['-H2O']*(i+1)) for i in range(n_H2O)])

                combined_shifts = [''.join(x) for x in itt.product(H2O_shifts, NH3_shifts)]
                frag_ions.extend([
                    {'frag':frag, 'ion_type':ion_type, 'shift':shift} for shift in combined_shifts
                    ])

        mono_ions = np.zeros(len(frag_ions), dtype=[
            ('Subtype', 'S1'),
            ('Ion Order', np.int32),
            ('Shift', 'S100'),
            ('m/z', np.float64),
            ('z', np.int32)
            ])

        for n,frag_ion in enumerate(frag_ions):
            subtype   = frag_ion['ion_type']
            ion_order = len(frag_ion['frag'])
            shift     = frag_ion['shift']
            mz        = sum(self.aa_weights[aa] for aa in frag_ion['frag']) \
                    + self.frag_shift[subtype] \
                    - (shift.count('NH3')*17.) \
                    - (shift.count('H2O')*18.)
            mono_ions[n] = (subtype, ion_order, shift, mz, 1)

        ions = mono_ions.copy()
        for c in self.charges:
            multi = mono_ions.copy()
            multi['z'] = c
            multi['m/z'] = (multi['m/z']+c-1)/c
            ions = np.concatenate((ions, multi))

        ions = np.sort(ions, order=['m/z'])
        self.ions = ions
        return ions


class ScoreReport:
    def __init__(self, out_report_path):
        self.out_report_path = os.path.abspath(out_report_path)

    def write(self, data):
        with open(self.out_report_path, 'a') as f:
            f.write('\t'.join([str(x) for x in data])+'\n')

def perm(unplaced, prefix):
    if unplaced:
        for element in unplaced:
            yield from perm(unplaced - Counter(element), prefix + element)
    else:
        yield prefix

def permutations(iterable):
    yield from perm(Counter(iterable), '')
