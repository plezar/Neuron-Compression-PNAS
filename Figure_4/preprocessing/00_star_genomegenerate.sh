#!/bin/bash

#$ -M jnajera2@nd.edu # Email address for job notification
#$ -m abe # Send mail when job begins, ends and aborts
#$ -pe mpi-24 48 # Specify parallel environment and legal core size
#$ -q debug # Specify queue
#$ -N genomeGenerate # Specify job name

module load bio/star/2.7.2

GENOMEDIR=&quot;/afs/crc/group/TIMELab/data/ILMN_1692_Datta_ND_RNAseq6_Apr2023/genome&quot;

STAR --runThreadN 12 \
    --runMode genomeGenerate \
    --genomeDir $GENOMEDIR/STAR \
    --genomeFastaFiles $GENOMEDIR/GRCh38.primary_assembly.genome.fa \
    --sjdbOverhang 149 \
    --sjdbGTFfile $GENOMEDIR/gencode.v43.primary_assembly.annotation.gtf