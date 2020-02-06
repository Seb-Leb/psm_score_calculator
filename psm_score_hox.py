import os
import sys
import argparse
import pickle
import multiprocessing as mp
from score_module import PSMReport, Score, ScoreReport, TheoreticalSpectrum

parser = argparse.ArgumentParser(description='Evaluate PSM score from PeptideShaker reports.')
parser.add_argument('pickle_path', metavar='path_to_pickle_object', type=str,
        help='path to pickle obj.')
parser.add_argument('--output_path', metavar='output_path', type=str,
        help='Output path to the psm score report.')
parser.add_argument('--partial_report', metavar='partial_report', type=str,
        help='Path to impcomplete output psm score report for completion.')
parser.add_argument('--n_cpu', metavar='n_cpu', type=int, default=1,
        help='Number of cpus available for parallel processing of P-values.')
parser.add_argument('--n_random', metavar='n_random', type=int, default=10000,
        help='Number of random peptide to generate for P-value calculation.')
parser.add_argument('--tol', metavar='tol', type=float, default=0.02,
        help='Tolerence for MS-2 peak assignment in DA.')
parser.add_argument('--compute_pval', metavar='compute_pval', type=int, default=1,
        help='Tolerence for MS-2 peak assignment in DA.')
parser.add_argument('--replicate', metavar='replicate', type=str,
        help='replicate identificartion i.e. "exp1".')

args = parser.parse_args()

def psm_score_worker(psm, q):
    rand_peps = pickle.load(open('/home/sleblanc/rand_peps.pkl', 'rb'))
    pep_seq, spectrum, scan_number = psm
    S = Score(tol=0.02, n_random=10000)
    T = TheoreticalSpectrum()
    T.compute_spectrum(pep_seq)
    spec_ann = S.get_peaks(T.ions, spectrum)
    hscore = S.hyperscore(spec_ann, pep_seq, spectrum, compute_pval=compute_pval, rand_peps=rand_peps)
    if hscore:
        hscore, pval = hscore
    else:
        hscore, pval = 'null', 'null'
    q.put([scan_number, pep_seq, hscore, pval])
    return [scan_number, pep_seq, hscore, pval]

def psm_report_writer(fname, mode, q):
    with open(fname, mode) as f:
        if mode == 'w':
            f.write('\t'.join(['scan_number', 'pep_seq', 'hscore', 'pval'])+'\n')
            f.flush()
        while 1:
            m = q.get()
            if m == 'Done.':
                break
            f.write('\t'.join([str(x) for x in m])+'\n')
            f.flush()

def filter_out_previous(fpath, psms, db_name, exp):
    scan_nums = set()
    with open(fpath, 'r') as f:
        for n,l in enumerate(f):
            db_name, exp, scan_number, pep_seq, hscore, pval = l.strip().split('\t')
            scan_nums.add('|'.join([db_name, exp, scan_number]))
    print('Continuing after {} spectra scored.'.format(len(scan_nums)), flush=True)
    psms = [x for x in psms if '|'.join([db_name, exp, x['Spectrum Scan Number']]) not in scan_nums]
    return psms

if __name__ == "__main__":
    args = vars(args)
    compute_pval = 1
    mode = 'w'
    if args['partial_report']:
        score_report_path = os.path.abspath(args['partial_report'])
        mode = 'a'
    elif args['output_path']:
        score_report_path = os.path.abspath(args['output_path'])
    else:
        score_report_path = os.path.abspath('./{}.psm_score_report.tsv'.format(ps_report.split('/')[-1]))
    print('Writing to file: {}\n'.format(score_report_path), flush=True)

    psm_rep = pickle.load(open(args['pickle_path'], 'rb'))

    manager = mp.Manager()
    q = manager.Queue()
    pool = mp.Pool(args['n_cpu'])

    watcher = pool.apply_async(psm_report_writer, (score_report_path, mode, q))

    jobs = []
    exp = args['replicate']
    psms = psm_rep[exp]
    print('{} | Scoring {} psms.'.format(exp, len(psms)), flush=True)
    if mode == 'a':
        psms = filter_out_previous(score_report_path, psms)
    for psm in psms:
        pep_seq     = psm['Sequence']
        spectrum    = psm['Spectrum']
        scan_number = psm['Spectrum Scan Number']
        args = [pep_seq, spectrum, scan_number]
        job = pool.apply_async(psm_score_worker, (args, q))
        jobs.append(job)

    for job in jobs:
        job.get()

    q.put('Done.')
    pool.close()
    pool.join()

    print('Done.')
