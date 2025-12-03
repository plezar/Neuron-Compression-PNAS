#!/bin/bash
#$ -M mzarodn2@nd.edu   # Email address for job notification
#$ -m abe            # Send mail when job begins, ends and aborts
#$ -pe smp 12     # Specify parallel environment and legal core size
#$ -q long           # Specify queue
#$ -N QC01       # Specify job name

conda activate bioinfo
Rscript 01_QC.r
