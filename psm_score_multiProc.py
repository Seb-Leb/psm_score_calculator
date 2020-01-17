import os
import sys
import argparse
import multiprocessing as mp
from score_module import PSMReport, Score, ScoreReport, TheoreticalSpectrum

parser = argparse.ArgumentParser(description='Evaluate PSM score from PeptideShaker reports.')
parser.add_argument('pepshaker_report_path', metavar='pepshaker_report_path', type=str,
        help='absolute filepath to the PeptideShaker report containing all spectra with annoations.')
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

args = parser.parse_args()

def psm_score_worker(psm, q):
    pep_seq, spectrum, scan_number = psm
    print('Scoring {} with spectrum {}.'.format(pep_seq, scan_number), flush=True)
    S = Score(tol=0.02, n_random=10000)
    T = TheoreticalSpectrum()
    T.compute_spectrum(pep_seq)
    spec_ann = S.get_peaks(T.ions, spectrum)
    hscore = S.hyperscore(spec_ann, pep_seq, spectrum, compute_pval=compute_pval)
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

def filter_out_previous(fpath, psms):
    scan_nums = set()
    with open(fpath, 'r') as f:
        for n,l in enumerate(f):
            scan_number, pep_seq, hscore, pval = l.strip().split('\t')
            scan_nums.add(scan_number)
    print('Continuing after {} spectra scored.'.format(len(scan_nums)), flush=True)
    psms = [x for x in psms if x['Spectrum Scan Number'] not in scan_nums]
    return psms

if __name__ == "__main__":
    args = vars(args)
    compute_pval = 1

    ps_report    = os.path.abspath(args['pepshaker_report_path'])
    if not os.path.exists(ps_report) or not os.path.isfile(ps_report):
        print('Provide full path to the peptide shaker report to be scored.')

    mode = 'w'
    if args['partial_report']:
        score_report_path = os.path.abspath(args['partial_report'])
        mode = 'a'
    elif args['output_path']:
        score_report_path = os.path.abspath(args['output_path'])
    else:
        score_report_path = os.path.abspath('./{}.psm_score_report.tsv'.format(ps_report.split('/')[-1]))
    print('Writing to file: {}\n'.format(score_report_path), flush=True)

    psm_rep = PSMReport(ps_report)
    psms = [v for k,v in psm_rep.psms.items()]
    if mode == 'a':
        psms = filter_out_previous(score_report_path, psms)
    print('Scoring {} psms.'.format(len(psms)), flush=True)

    manager = mp.Manager()
    q = manager.Queue()
    pool = mp.Pool(args['n_cpu'])

    watcher = pool.apply_async(psm_report_writer, (score_report_path, mode, q))

    jobs = []
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
