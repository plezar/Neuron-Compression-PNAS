1. The NCBI reference genomes were downloaded from [the NCBI Datasets portal](https://www.ncbi.nlm.nih.gov/datasets/genome) by searching for the appropriate species and retrieving the corresponding RefSeq dataset.

2. All raw FASTQ files were merged using the following script: `00_merge_fastq.sh`.

3. Reads were separated by species using BBSplit: `01_bbsplit.sh`. The output consists of paired, interleaved FASTQ files.

4. STAR genome indices were generated using: `02_star_genomegenerate.sh ${home_dir}/genomes/mmus/refseq/GCF_000001635.27/ncbi_dataset/data/GCF_000001635.27 GCF_000001635.27_GRCm39_genomic.fna genomic.gtf`

5. Gene counts were generated using STAR

```bash
while read p; do
 qsub 03_star_genecounts.sh \
    $p \
    ${home_dir}/Neurons_Sep2024/bbmap \
    ${home_dir}/Neurons_Sep2024/STAR_out \
    ${home_dir}/genomes/mmus/refseq/GCF_000001635.27/ncbi_dataset/data/GCF_000001635.27/genomic.gtf \
    ${home_dir}/genomes/mmus/refseq/GCF_000001635.27/ncbi_dataset/data/GCF_000001635.27/STAR/
done < ${home_dir}/Neurons_Sep2024/sample_names_1019_compress_mmus.txt
```

```bash
while read p; do
qsub 03_star_genecounts.sh \
    $p \
    ${home_dir}/Neurons_Sep2024/bbmap \
    ${home_dir}/Neurons_Sep2024/STAR_out \
    ${home_dir}/genomes/hsapiens/refseq/ncbi_dataset/data/GCF_000001405.40/genomic.gtf \
    ${home_dir}/genomes/hsapiens/refseq/ncbi_dataset/data/GCF_000001405.40/STAR/
done < ${home_dir}/Neurons_Sep2024/sample_names_hs.txt
```