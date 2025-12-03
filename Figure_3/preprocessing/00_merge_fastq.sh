#$ -M chef.maximiliano@gmail.com   # Email address for job notification
#$ -m abe            # Send mail when job begins, ends and aborts
#$ -pe smp 1     # Specify parallel environment and legal core size
#$ -q long           # Specify queue
#$ -N merge_fastq       # Specify job name

run1_dir=/afs/crc/group/TIMELab/Neurons_Sep2024/Datta-MZ-3042_240903/fastq
run2_dir=/afs/crc/group/TIMELab/Neurons_Sep2024/Datta-MZ-3042_240905/fastq
s=$1
cat "${run1_dir}/${s}_R1_001.fastq.gz" "${run2_dir}/${s}_R1_001.fastq.gz" > "/scratch365/mzarodn2/Neurons_Sep2024/merged_fastq/${s}_R1_001.fastq.gz"
cat "${run1_dir}/${s}_R2_001.fastq.gz" "${run2_dir}/${s}_R2_001.fastq.gz" > "/scratch365/mzarodn2/Neurons_Sep2024/merged_fastq/${s}_R2_001.fastq.gz"
