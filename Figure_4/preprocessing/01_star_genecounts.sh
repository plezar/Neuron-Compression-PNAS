#!/bin/bash

#$ -M jnajera2@nd.edu # Email address for job notification
#$ -m abe # Send mail when job begins, ends and aborts
#$ -pe mpi-24 24 # Specify parallel environment and legal core size
#$ -q debug # Specify queue
#$ -N COUNT # Specify job name

# this remains the same for all jobs (specific to mapping)

OUTPREFIX=$1
NCPU=12
INPUTDIR=&quot;/afs/crc.nd.edu/group/TIMELab/ILMN_1857_Harker_ND_RNAseq6_Aug2023”
OUTDIR=&quot;/afs/crc.nd.edu/group/TIMELab/ILMN_1857_Harker_ND_RNAseq6_Aug2023/counts&quot;
GTFFILE=&quot;/afs/crc.nd.edu/group/TIMELab/ILMN_1857_Harker_ND_RNAseq6_Aug2023/gencode.v44
.basic.annotation.gtf&quot;
GENOMEDIR=&quot;/afs/crc.nd.edu/group/TIMELab/ILMN_1857_Harker_ND_RNAseq6_Aug2023/STAR”

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
    --readFilesIn ${INPUTDIR}/fastq/${OUTPREFIX}_R1_001.fastq
    ${INPUTDIR}/fastq/${OUTPREFIX}_R2_001.fastq \
    --outFileNamePrefix $OUTPREFIX

#--readFilesCommand zcat \