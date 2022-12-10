#!/bin/bash
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -N sc_pip_1
#$ -pe smp 8
#$ -l h_vmem=4G
#$ -q processing.q

export LC_ALL=C
export MALLOC_ARENA_MAX=4

compression_level=5
fasta_version=hg38
bwa_version="0.7.17-r1188"
bwa_commandline="bwa mem -K 100000000 -p -v 3 -t 8 -Y"

ref_dir=~/references/homo_sapiens/GRCh38.broad

ref_fasta=${ref_dir}/references-hg38-v0-Homo_sapiens_assembly38.fasta
dbSNP_vcf=${ref_dir}/references-hg38-v0-Homo_sapiens_assembly38.dbsnp138.vcf
m1indel_vcf=${ref_dir}/references-hg38-v0-Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
knindel_vcf=${ref_dir}/references-hg38-v0-Homo_sapiens_assembly38.known_indels.vcf.gz

### references ###
## gatk4 data preprocessing:
# https://github.com/gatk-workflows/gatk4-data-processing
# https://gatk.broadinstitute.org/hc/en-us/articles/360035535912-Data-pre-processing-for-variant-discovery
## reference data:
# referene data: https://console.cloud.google.com/storage/browser/genomics-public-data/references/hg38/v0
## no read trimming any more:
# https://evodify.com/gatk-in-non-model-organism/
 
if [ "$#" -ne 2 ]; then
 echo "Usage: SAMPLEID bulk/cell" >&2
 exit 1
fi

sn=${1}

date
echo "==START==", ${sn}

mkdir -p ./${sn}/
mkdir -p ./${sn}/tmp
cd ./${sn}/

### sub-step 1. Quality control of raw sequencing read ###
fastqc -t 8 -o ./ ../fastq/${sn}_1.fq.gz ../fastq/${sn}_2.fq.gz

### sub-step 2. Convert fastq files to unmapped bam file ###
gatk --java-options "-Dsamjdk.compression_level=${compression_level} -Xms3000m" \
  FastqToSam \
  -F1 ../fastq/${sn}_1.fq.gz -F2 ../fastq/${sn}_2.fq.gz \
  --COMPRESSION_LEVEL 5 \
  --PLATFORM ILLUMINA \
  --OUTPUT ${sn}.unmapped.bam \
  --TMP_DIR ./tmp \
  --SAMPLE_NAME ${sn}

### sub-step 3. Mark sequencing adapters in unmapped bam file ###
# https://gatk.broadinstitute.org/hc/en-us/articles/360039568932--How-to-Map-and-clean-up-short-read-sequence-data-efficiently
gatk --java-options "-Dsamjdk.compression_level=${compression_level} -Xms3000m" \
  MarkIlluminaAdapters \
  -I ${sn}.unmapped.bam \
  -O ${sn}.markilluminaadapters.bam \
  -M ${sn}.markilluminaadapters_metrics.txt \
  --TMP_DIR `pwd`/tmp

### sub-step 4. Sequence alignment ###
gatk --java-options "-Dsamjdk.compression_level=${compression_level} -Xms3000m" \
  SamToFastq \
  --TMP_DIR `pwd`/tmp \
  --INPUT ${sn}.markilluminaadapters.bam \
  --FASTQ /dev/stdout \
  --CLIPPING_ATTRIBUTE XT \
  --CLIPPING_ACTION 2 \
  --INTERLEAVE true \
  --NON_PF true \
| \
bwa mem -K 100000000 -p -v 3 -t 8 -Y ${ref_fasta} /dev/stdin - 2> >(tee ${sn}.bwa.stderr.log >&2) \
| \
samtools view -1 - > ${sn}.${fasta_version}.bwa.bam

rm ${sn}.markilluminaadapters.bam

### sub-step 5. Merge mapped bam file with unmapped bam file (ubam) ###
gatk --java-options "-Dsamjdk.compression_level=${compression_level} -Xms3000m" \
  MergeBamAlignment \
  --TMP_DIR ./tmp \
  --VALIDATION_STRINGENCY SILENT \
  --EXPECTED_ORIENTATIONS FR \
  --ATTRIBUTES_TO_RETAIN X0 \
  --ALIGNED_BAM ${sn}.${fasta_version}.bwa.bam \
  --UNMAPPED_BAM ${sn}.unmapped.bam \
  --OUTPUT ${sn}.${fasta_version}.merged.bam \
  --REFERENCE_SEQUENCE ${ref_fasta} \
  --PAIRED_RUN true \
  --SORT_ORDER "unsorted" \
  --IS_BISULFITE_SEQUENCE false \
  --ALIGNED_READS_ONLY false \
  --CLIP_ADAPTERS false \
  --CLIP_OVERLAPPING_READS true \
  --INCLUDE_SECONDARY_ALIGNMENTS true \
  --MAX_RECORDS_IN_RAM 2000000 \
  --ADD_MATE_CIGAR true \
  --MAX_INSERTIONS_OR_DELETIONS -1 \
  --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
  --PROGRAM_RECORD_ID "bwamem" \
  --PROGRAM_GROUP_VERSION "${bwa_version}" \
  --PROGRAM_GROUP_COMMAND_LINE "${bwa_commandline}" \
  --PROGRAM_GROUP_NAME "bwamem" \
  --UNMAPPED_READ_STRATEGY COPY_TO_TAG \
  --ALIGNER_PROPER_PAIR_FLAGS true \
  --UNMAP_CONTAMINANT_READS true

rm ${sn}.${fasta_version}.bwa.* ${sn}.unmapped.*

### sub-step 6. Sort BAM file by coordinate order and fix tag values for NM and UQ ###
gatk --java-options "-Dsamjdk.compression_level=${compression_level} -Xms4000m" \
  SortSam \
  --TMP_DIR ./tmp \
  --INPUT ${sn}.${fasta_version}.merged.bam \
  --OUTPUT ${sn}.${fasta_version}.sortsam.bam \
  --SORT_ORDER "coordinate" \
  --CREATE_INDEX false \
  --CREATE_MD5_FILE false

gatk --java-options "-Dsamjdk.compression_level=${compression_level} -Xms4000m" \
  SetNmMdAndUqTags \
  --TMP_DIR ./tmp \
  --INPUT ${sn}.${fasta_version}.sortsam.bam \
  --OUTPUT ${sn}.${fasta_version}.setnmmdanduqtags.bam \
  --CREATE_INDEX true \
  --CREATE_MD5_FILE true \
  --REFERENCE_SEQUENCE ${ref_fasta}

rm ${sn}.${fasta_version}.merged.* ${sn}.${fasta_version}.sortsam.*

### sub-step 7. Mark duplicate reads to avoid counting non-independent observations ###
gatk --java-options "-Dsamjdk.compression_level=${compression_level} -Xms4000m" \
  MarkDuplicates \
  --TMP_DIR ./tmp \
  --INPUT ${sn}.${fasta_version}.setnmmdanduqtags.bam \
  --OUTPUT ${sn}.${fasta_version}.markduplicates.bam \
  --METRICS_FILE ${sn}.hg38.duplicate_metrics \
  --VALIDATION_STRINGENCY SILENT \
  --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
  --ASSUME_SORT_ORDER "queryname" \
  --CREATE_MD5_FILE true

rm ${sn}.${fasta_version}.setnmmdanduqtags.*

### sub-step 8. Generate Base Quality Score Recalibration (BQSR) model ###
gatk --java-options "-Xms4000m" \
  BaseRecalibrator \
  -R ${ref_fasta} \
  -I ${sn}.${fasta_version}.markduplicates.bam \
  --use-original-qualities \
  -O ${sn}.${fasta_version}.recal_data.csv \
  --known-sites ${dbSNP_vcf} \
  --known-sites ${m1indel_vcf} \
  --known-sites ${knindel_vcf}

### sub-step 9. Apply Base Quality Score Recalibration (BQSR) model ###
gatk --java-options "-Dsamjdk.compression_level=${compression_level} -Xms4000m" \
  ApplyBQSR \
  -R ${ref_fasta} \
  -I ${sn}.${fasta_version}.markduplicates.bam \
  -O ${sn}.${fasta_version}.bqsr.bam \
  -bqsr ${sn}.${fasta_version}.recal_data.csv \
  --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
  --add-output-sam-program-record \
  --create-output-bam-md5 \
  --use-original-qualities

rm ${sn}.${fasta_version}.markduplicates.*

### sub-step 10. Sort by coordinate, index, and calculate md5 ###
gatk --java-options "-Dsamjdk.compression_level=${compression_level}" \
  SortSam \
  --TMP_DIR ./tmp \
  --INPUT ${sn}.${fasta_version}.bqsr.bam \
  --OUTPUT ${sn}.${fasta_version}.gatk4.bam \
  --SORT_ORDER "coordinate" 

rm ${sn}.${fasta_version}.bqsr.*

gatk --java-options "-Xms4000m" \
  BuildBamIndex \
  -I ${sn}.${fasta_version}.gatk4.bam \
  --TMP_DIR ./tmp \

md5sum ${sn}.${fasta_version}.gatk4.bam > ${sn}.${fasta_version}.gatk4.bam.md5

rm -r ./tmp

### sub-step 11. Calculate depth and coverage ###
mapq=40
mkdir -p ../coverage/
cd ../coverage/
samtools depth -Q ${mapq} ../${sn}/${sn}.hg38.gatk4.bam | bzip2 -z9 > ${sn}.depth.bz2

a20=$(bzcat ${sn}.depth.bz2 | awk '$3>=20' | wc -l)
a15=$(bzcat ${sn}.depth.bz2 | awk '$3>=15' | wc -l)
a10=$(bzcat ${sn}.depth.bz2 | awk '$3>=10' | wc -l)
a5=$(bzcat ${sn}.depth.bz2 | awk '$3>=5' | wc -l)
a1=$(bzcat ${sn}.depth.bz2 | awk '$3>=1' | wc -l)
d=$(bzcat ${sn}.depth.bz2 | awk '{sum += $3};END {print sum}')
rm ${sn}.depth.bz2

echo ${sn},$d,$a20,$a15,$a10,$a5,$a1 >> summary_coverage.csv

### sub-step 12. Call germline heterozygous SNVs and INDELs (only for bulk DNA) ###
if [ "$2" = "bulk" ]
then
  cd ../
  mkdir -p ./ht
  cd ./ht

  gatk --java-options "-Xms4000m" \
    HaplotypeCaller \
    -R ${ref_fasta} \
    -O ./${sn}.gvcf.gz \
    --dbsnp ${dbSNP_vcf} \
    -I ../${sn}/${sn}.${fasta_version}.gatk4.bam \
    -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation \
    -ERC GVCF

  zcat ${sn}.gvcf.gz | grep -v '##' \
    | awk 'length($3)>1 && length($4)==1 && $5 ~ /^[[:alpha:]],<NON_REF>/ && $7=="."' \
    | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10}' \
    | grep '0/1' \
    | sort -k1,1d -k2,2n > ${sn}.${fasta_version}.hsnp.biallelic.dbsnp.vcf

  zcat ${sn}.gvcf.gz | grep -v '##' \
    | awk 'length($3)>1 && $7=="." && ((length($4)>1 && $5 ~ /^[[:alpha:]]+,<NON_REF>/) || $5 ~ /^[[:alpha:]][[:alpha:]]+,<NON_REF>/)' \
    | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10}' \
    | grep '0/1' \
    | sort -k1,1d -k2,2n > ${sn}.${fasta_version}.hindel.biallelic.dbsnp.vcf

  cd ../
  mkdir -p ./tmp
  mkdir -p ./sensitivity

  # for sensitivity calculation only limit to 20x in bulk to make sure the hSNPs and hINDELs are real
  awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $10}' ./ht/${sn}.${fasta_version}.hsnp.biallelic.dbsnp.vcf \
    | awk -F "[:,]" '$3 + $4 >= 20' | sort -T ./tmp -k1,1d -k2,2n | grep -v "#" | awk '{print $1 "\t" $2 "\t" $2 "\t" $4 "\t" $5 "\t" $6}' \
    > ./sensitivity/${sn}.hsnp.biallelic.dbsnp.20x.sorted.bed

  awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $10}' ./ht/${sn}.${fasta_version}.hindel.biallelic.dbsnp.vcf \
    | awk -F "[:,]" '$3 + $4 >= 20' | sort -T ./tmp -k1,1d -k2,2n | grep -v "#" | awk '{print $1 "\t" $2 "\t" $2 "\t" $4 "\t" $5 "\t" $6}' \
    > ./sensitivity/${sn}.hindel.biallelic.dbsnp.20x.sorted.bed

else
  echo "skip haplotype caller for non bulk samples"
fi

echo "===END===",${sn}
date

