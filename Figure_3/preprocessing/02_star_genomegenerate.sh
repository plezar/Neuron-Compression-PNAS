#!/bin/bash

#$ -M mzarodn2@nd.edu   # Email address for job notification
#$ -m abe            # Send mail when job begins, ends and aborts
#$ -pe mpi-24 48     # Specify parallel environment and legal core size
#$ -q debug           # Specify queue
#$ -N genomeGenerate       # Specify job name

module load bio/star/2.7.2

GENOMEDIR=$1
genomeFastaFile=$2
sjdbGTFfile=$3

mkdir -p $GENOMEDIR

STAR --runThreadN 12 \
		--runMode genomeGenerate \
		--genomeDir $GENOMEDIR/STAR \
		--genomeFastaFiles $GENOMEDIR/$genomeFastaFile \
		--sjdbOverhang 149 \
		--sjdbGTFfile $GENOMEDIR/$sjdbGTFfile
