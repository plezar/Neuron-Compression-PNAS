#!/bin/bash

#$ -M mzarodn2@nd.edu   # Email address for job notification
#$ -m abe            # Send mail when job begins, ends and aborts
#$ -pe smp 24     # Specify parallel environment and legal core size
#$ -q long           # Specify queue
#$ -N BBmap       # Specify job name

dir="/scratch365/mzarodn2/Neurons_Sep2024/merged_fastq"
out_dir="/afs/crc/group/TIMELab/Neurons_Sep2024/bbmap"
s=$1

echo "SAMPLE: $s"
cd /afs/crc.nd.edu/user/m/mzarodn2/Private/rnaseq/Neurons_Sep2024

bash /afs/crc.nd.edu/user/m/mzarodn2/Private/soft/bbmap/bbsplit.sh \
                in1="$dir/${s}_R1_001.fastq.gz" \
                in2="$dir/${s}_R2_001.fastq.gz" \
                ref_hu=/afs/crc/group/TIMELab/genomes/hsapiens/refseq/ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna \
                ref_mmus=/afs/crc/group/TIMELab/genomes/mmus/refseq/GCF_000001635.27/ncbi_dataset/data/GCF_000001635.27/GCF_000001635.27_GRCm39_genomic.fna \
                ambiguous2=toss \
                out_hu="${out_dir}/${s}_hs_R1_001.fastq.gz" \
                out_mmus="${out_dir}/${s}_mmus_R1_001.fastq.gz"

# while read p; do qsub bbsplit.sh; done</afs/crc/group/TIMELab/Neurons_Sep2024/I27_sample_ids.txt
