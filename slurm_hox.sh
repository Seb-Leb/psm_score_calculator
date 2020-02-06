#!/bin/bash

#SBATCH --account=def-xroucou
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=40
#SBATCH --mem-per-cpu=2G

python3 /home/sleblanc/psm_score_calculator/psm_score_hox.py /home/sleblanc/merged_hox.pkl --output_path /home/sleblanc/hox_exp1_scores.tsv --n_cpu 40 --replicate exp1
