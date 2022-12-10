#!/bin/bash
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -N sc_pip_2
#$ -pe smp 8
#$ -l h_vmem=8G
#$ -q processing.q

export LC_ALL=C
export MALLOC_ARENA_MAX=4

if [ "$#" -ne 2 ]; then
  echo "Usage: $0 BULKID CELLID" >&2
  exit 1
fi

fasta_version=hg38

bulkid=${1}
cellid=${2}

mkdir -p ./sccaller

# reference genome (fasta)
ref_dir=~/references/homo_sapiens/GRCh38.broad
ref_fasta=${ref_dir}/references-hg38-v0-Homo_sapiens_assembly38.fasta

# dbSNP
dbSNP_bed=${ref_dir}/references-hg38-v0-Homo_sapiens_assembly38.dbsnp138.sorted.bed

# cell (bam)
cell=`pwd`/${cellid}/${cellid}.${fasta_version}.gatk4.bam
# cell (output vcf)
cell_vcf=`pwd`/sccaller/${cellid}.sccaller.vcf
# bulk (bam)
bulk=`pwd`/${bulkid}/${bulkid}.${fasta_version}.gatk4.bam
# hSNP detected from bulk (vcf / bed)
hsnp=`pwd`/ht/${bulkid}.${fasta_version}.hsnp.biallelic.dbsnp.vcf
# number of CPUs
ncpu=8

date
echo "==START==",${bulkid},${cellid}

### sub-step 1: Run SCcaller on all candidate mutations ###
python sccaller_v2.0.0.py \
  --bam ${cell} \
  --fasta ${ref_fasta} \
  --bulk ${bulk} \
  --output ${cell_vcf} \
  --snp_type hsnp \
  --snp_in ${hsnp} \
  --cpu_num ${ncpu} \
  --bias 0.75 \
  --wkdir ./sccaller \
  --engine samtools


### sub-step 2: Filter out low quality and germline calls and keep only somatic mutations ###
input=${cellid}.sccaller.vcf
o1=${cellid}.somatic.snv.vcf
o2=${cellid}.somatic.indel.vcf

cd ./sccaller/
# for snv calling
grep '0/1' ${input} | grep 'True' | awk '$7=="." && length($5)==1' | awk -F "[:,]" '$8+$9>=20' > ${o1}
# for indel requiring variant calling Q>=30
grep '0/1' ${input} | grep 'True' | awk '$7=="." && length($5)>1 && $6>=30' | awk -F "[:,]" '$8+$9>=20' > ${o2}


### sub-step 3: Calculate sensitivity ###
cd ../
sort -T ./tmp -k1,1d -k2,2n ./sccaller/${cellid}.sccaller.vcf \
  | grep -v "#" | awk '{print $1 "\t" $2 "\t" $2 "\t" $4 "\t" $5 "\t" $7 "\t" $10}' \
  | awk '$6=="."' > ./sensitivity/${cellid}.sccaller.sorted.bed

bedtools intersect -sorted -f 1 -wo \
  -a ./sensitivity/${cellid}.sccaller.sorted.bed \
  -b ./sensitivity/${bulkid}.hsnp.biallelic.dbsnp.20x.sorted.bed \
  | awk -F "[:,]" '$3+$4>=20' \
  | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7}' \
  > ./sensitivity/${cellid}.hsnp.20x.bed

x1=$(grep "0/1" ./sensitivity/${cellid}.hsnp.20x.bed | wc -l)
x2=$(wc -l ./sensitivity/${cellid}.hsnp.20x.bed | awk '{print $1}')

bedtools intersect -sorted -f 1 -wo \
  -a ./sensitivity/${cellid}.sccaller.sorted.bed \
  -b ./sensitivity/${bulkid}.hindel.biallelic.dbsnp.20x.sorted.bed \
  | awk -F "[:,]" '$3+$4>=20' \
  | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7}' \
  > ./sensitivity/${cellid}.hindel.20x.bed

x3=$(grep "0/1" ./sensitivity/${cellid}.hindel.20x.bed | wc -l)
x4=$(wc -l ./sensitivity/${cellid}.hindel.20x.bed | awk '{print $1}')

echo ${cellid},${x1},${x2},${x3},${x4} >> ./sensitivity/sensitivity.txt


### sub-step 4: Calculate coverage of both bulk and single cell ###
mapq=40
cd ./coverage/

a=$(samtools mpileup -q ${mapq} ../${bulkid}/${bulkid}.hg38.gatk4.bam ../${cellid}/${cellid}.hg38.gatk4.bam | awk '{print $1 "\t" $2 "\t" $4 "\t" $7}' | awk '$3>=20 && $4>=20' | wc -l)

echo ${cellid},$a >> summary_coverage_both.csv


### sub-step 5: Estimate mutation burden per single cell ###
cd ../
mkdir -p mutationburden/
cd mutationburden/

ref_size=$(awk '{sum += $2};END {print sum}'  ${ref_fasta}.fai)
co_size=${a}

no_snv=$(wc -l ../sccaller/${cellid}.somatic.snv.vcf | awk '{print $1}')
se_snv=$(echo "${x1} ${x2}" | awk '{printf ("%.2f\n",$1/$2)}')
burden_snv=$(echo "${no_snv} ${x1} ${x2} ${ref_size} ${co_size}" | awk '{printf ("%.1f\n",$1/$2*$3*$4/$5)}')

no_indel=$(wc -l ../sccaller/${cellid}.somatic.indel.vcf | awk '{print $1}')
se_indel=$(echo "${x3} ${x4}" | awk '{printf ("%.2f\n",$1/$2)}')
burden_indel=$(echo "${no_indel} ${x3} ${x4} ${ref_size} ${co_size}" | awk '{printf ("%.1f\n",$1/$2*$3*$4/$5)}' )

if [ -e mutation_burden.csv ]
then
  echo "ok"
else
  echo "bulk_id,cell_id,ref_size,bp_covered(20x),snv_sensitivity,indel_sensitivity,no_snv_observed,no_indel_observed,snv_burden,indel_burden" >> mutation_burden.csv
fi

echo ${bulkid},${cellid},${ref_size},${co_size},${se_snv},${se_indel},${no_snv},${no_indel},${burden_snv},${burden_indel} >> mutation_burden.csv

echo "==END==",${bulkid},${cellid}
date

