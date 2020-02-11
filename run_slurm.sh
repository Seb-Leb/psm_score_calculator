#!/bin/bash

#SBATCH --account=def-xroucou
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=40
#SBATCH --mem-per-cpu=2G

python3 /home/sleblanc/psm_score_calculator/psm_score_multiProc_pickle.py /home/sleblanc/scores_Classic_Refs.pkl --output_path /home/sleblanc/Classic_Refs_scored_psms.tsv --n_cpu 40 --db Ref --method mvhscore
