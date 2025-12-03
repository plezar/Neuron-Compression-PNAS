#!/bin/bash

#$ -M mzarodn2@nd.edu   # Email address for job notification
#$ -m abe            # Send mail when job begins, ends and aborts
#$ -pe mpi-24 24     # Specify parallel environment and legal core size
#$ -q debug           # Specify queue
#$ -N STAR       # Specify job name

# this remains the same for all jobs (specific to mapping)
OUTPREFIX=$1
NCPU=12
INPUTDIR=$2
OUTDIR=$3
GTFFILE=$4
GENOMEDIR=$5


module load bio/star/2.7.2
mkdir -p $OUTDIR
cd $OUTDIR

STAR \
	--outSAMattributes All \
	--outSAMtype BAM SortedByCoordinate \
	--quantMode GeneCounts \
	--runThreadN $NCPU \
	--sjdbGTFfile $GTFFILE \
	--outReadsUnmapped Fastx \
	--outMultimapperOrder Random \
	--outWigType wiggle \
	--genomeDir $GENOMEDIR \
	--readFilesIn ${INPUTDIR}/${OUTPREFIX}_R1_001.fastq ${INPUTDIR}/${OUTPREFIX}_R2_001.fastq \
	--outFileNamePrefix $OUTPREFIX
#--readFilesCommand zcat \
