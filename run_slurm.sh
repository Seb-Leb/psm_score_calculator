#!/bin/bash

#SBATCH --account=def-xroucou
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=500M

python3 /nfs3_ib/ip32-ib/home/sleblanc/psm_score_calculator/psm_score_multiProc.py /nfs3_ib/ip32-ib/home/xroucou_group/echange_de_fichiers/Annotation_PSM_complete_examplefile.txt --partial_report /nfs3_ib/ip32-ib/home/sleblanc/test_psm_scoring_slurm.tsv --n_cpu 24 --n_random 10000
