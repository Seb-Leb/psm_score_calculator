#!/bin/bash

#SBATCH --account=def-xroucou
#SBATCH --time=07:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=500M

python3 /nfs3_ib/ip32-ib/home/sleblanc/psm_score_calculator/psm_score_multiProc.py /nfs3_ib/ip32-ib/home/xroucou_group/echange_de_fichiers/Annotation/Annotation_OpenProt1peptide_3x_1_exp1.txt --partial_report /nfs3_ib/ip32-ib/home/xroucou_group/echange_de_fichiers/Annotation/Annotation_OpenProt1peptide_3x_1_exp1_scores.tsv --n_cpu 24 --n_random 10000
