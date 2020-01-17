import os
import sys
import argparse
from score_module import PSMReport, Score, ScoreReport, TheoreticalSpectrum

parser = argparse.ArgumentParser(description='Evaluate PSM score from PeptideShaker reports.')
parser.add_argument('pepshaker_report_path', metavar='pepshaker_report_path', type=str,
        help='absolute filepath to the PeptideShaker report containing all spectra with annoations.')
parser.add_argument('output_path', metavar='output_path', type=str,
        help='Output path to the psm score report.')
parser.add_argument('--n_cpu', metavar='n_cpu', type=int, default=1,
        help='Number of cpus available for parallel processing of P-values.')
parser.add_argument('--n_random', metavar='n_random', type=int, default=10000,
        help='Number of random peptide to generate for P-value calculation.')
parser.add_argument('--tol', metavar='tol', type=float, default=0.02,
        help='Tolerence for MS-2 peak assignment in DA.')
parser.add_argument('--compute_pval', metavar='compute_pval', type=int, default=1,
        help='Tolerence for MS-2 peak assignment in DA.')

args = parser.parse_args()


if __name__ == "__main__":
    args = vars(args)
    compute_pval      = args['compute_pval']
    score_report_path = os.path.abspath(args['output_path'])
    if os.path.isfile(score_report_path):
        print('The psm scores report file {} already exists.\n'.format(score_report_path))
        while True:
            resp = input('Overwrite? (y/n) : ')
            if resp == 'n':
                sys.exit()
            if resp == 'y':
                os.remove(score_report_path)
                break
    out_report   = ScoreReport(score_report_path)
    ps_report    = os.path.abspath(args['pepshaker_report_path'])
    if not os.path.exists(ps_report) or not os.path.isfile(ps_report):
        print('Provide full path to the peptide shaker report to be scored.')

    psm_rep = PSMReport(ps_report)
    S = Score(
            n_cpu=args['n_cpu'],
            tol=args['tol'],
            n_random=args['n_random']
            )
    print('Scoring {} psms.'.format(len(psm_rep.psms)))
    T = TheoreticalSpectrum()
    for psm_n in psm_rep.psms:
        pep_seq     = psm_rep.psms[psm_n]['Sequence']
        spectrum    = psm_rep.psms[psm_n]['Spectrum']
        scan_number = psm_rep.psms[psm_n]['Spectrum Scan Number']
        T.compute_spectrum(pep_seq)
        spec_ann = S.get_peaks(T.ions, spectrum)
        hscore = S.hyperscore(spec_ann, pep_seq, spectrum, compute_pval=compute_pval)
        if compute_pval:
            hscore, pval = hscore
            out_report.write([scan_number, pep_seq, hscore, pval])
        else:
            out_report.write([scan_number, pep_seq, hscore])
    print('Done.')
