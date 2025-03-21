#!/usr/bin/bash
# split_dedup.sh - a bash script that deduplicates .bam files, aligned with STARsolo
#                - deduplicates chr-by-chr, to reduce run time and memory requirements
# Usage:
# bash chr_split_dedup.sh /path/to/input.bam /path/to/whitelist num_cores /path/to/output/directory /path/to/tmp/directory file.LOG

# chr_split_dedup (v1.0) - input a .bam file, output the deduplicated .bam file in OUTDIR

#TODO: python script instead of bash to make parallelization a bit easier?

INBAM=$1 # path to .bam file (sorted & indexed already!)
BC=$2 # path to barcode whitelist
CORE=$3 # number of cores for parallelization
OUTBAM=$4 # output/deduped bam path
TMPDIR=$5
LOG=$6

mkdir -p ${TMPDIR}
cd ${TMPDIR}

PREFIX=$(basename ${INBAM} .bam)

echo ".bam file location:    ${INBAM}"
echo "Max cores:             ${CORE}"
echo "Output location:       ${OUTBAM}"
echo

# Check params...
if [ ! -f ${INBAM} ]; then
    echo "Can't find ${INBAM}"
    exit 1
fi

# Check for .bam index file
if [ ! -f ${INBAM}.bai ]; then
    echo "Indexing input .bam file because you forgot to..."
    samtools index -@ ${CORE} ${INBAM}
fi

# Remove reads that don't have a barcode (CB)
echo "Removing reads without 'CB' or 'UB' tags..."
date

samtools view -h -@ ${CORE} --tag-file CB:${BC} ${INBAM} \
| grep -v "UB:Z:-" \
| samtools view -bS \
> ${TMPDIR}/filter.bam

echo

# Split bam by chromosome/strand
echo "Splitting by chromosome..."
date

bamtools split -in ${TMPDIR}/filter.bam -reference

echo

BAMLIST=${TMPDIR}/BAMLIST.tmp
ls *REF_*.bam > ${BAMLIST}

# Index split .bam's
echo "Indexing split .bam files..."
date

while read BAM; do
  samtools index -@ ${CORE} ${BAM}
done < ${BAMLIST}

echo

# Dedup resulting bams one-by-one
echo "Deduplicating split .bam files..."
date

while read BAM; do
  samtools index -@ ${CORE} ${BAM}

  umi_tools dedup \
  -I ${BAM} \
  --extract-umi-method=tag \
  --umi-tag=UB \
  --cell-tag=CB \
  --method=unique \
  --per-cell \
  --unmapped-reads=discard \
  -S ${TMPDIR}/dedup_${BAM}
done < ${BAMLIST}

echo

# Merge, sort, and index dedup'ed .bams
echo "Merging, sorting, and indexing deduplicated .bam files..."
date

samtools merge -f -o ${TMPDIR}/dedup_merge.bam -@ ${CORE} --write-index ${TMPDIR}/dedup_*REF_*.bam

samtools sort ${TMPDIR}/dedup_merge.bam > ${OUTBAM}
samtools index -@ ${CORE} ${OUTBAM} && rm -r ${TMPDIR}
