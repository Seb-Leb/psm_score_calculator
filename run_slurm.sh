#!/bin/bash

#SBATCH --account=def-xroucou
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=40
#SBATCH --mem-per-cpu=2G

python3 /home/sleblanc/psm_score_calculator/psm_score_multiProc_pickle.py /home/sleblanc/psms_to_score.pkl --partial_report /home/sleblanc/scored_psms.tsv --n_cpu 40
