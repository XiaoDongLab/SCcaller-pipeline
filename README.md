# SCcaller-pipeline
An pipeline for analyzing single-cell whole-genome sequencing data: from fastq files to mutation calls and mutation burdens.

Version: 1.0.0

Updated date: 2022.12.10

Citation:

Lei Zhang, Moonsook Lee, Alexander Y. Maslov, Jan Vijg, and Xiao Dong. Single-cell whole-genome sequencing for discovering somatic mutations. In submission.

#####
## Author and License

Author: Xiao Dong

Email: biosinodx@gmail.com (X.D.), dong0265@umn.edu (X.D.), spsc83@gmail.com (Y.W.)

Licensed under the GNU Affero General Public License version 3 or later

#####
## Dependencies

• bwa v0.7.17 (https://bio-bwa.sourceforge.net/)

• bedtools v2.30.0 (https://bedtools.readthedocs.io/en/latest/)

• FastQC v0.11.9 (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

• GATK v4.3.0.0 (https://gatk.broadinstitute.org/hc/en-us)

• Java v1.8 (https://www.java.com/en/)

• python v2.7.18 (https://www.python.org/)

• samtools v1.9 (http://www.htslib.org/)

• SCcaller v2.0.0 (https://github.com/biosinodx/SCcaller)

#####
## Reference files

All the reference files can be downloaded from the cloud storage of the Broad Institute at:
https://console.cloud.google.com/storage/browser/genomics-public-data/references/hg38/v0

• Reference genome: Homo_sapiens_assembly38.fasta

• dbSNP dataset: Homo_sapiens_assembly38.dbsnp138.vcf

• INDEL datasets: Homo_sapiens_assembly38.known_indels.vcf, Mills_and_1000G_gold_standard.indels.hg38.vcf

#####
## Pipeline

Current SCcaller-pipeline includes two major steps, and 17 sub-steps. See the flowchart, and Usage below.

![alt text](https://github.com/XiaoDongLab/SCcaller-pipeline/blob/main/sccaller_pipeline_flowchart_v1.0.0.png)

#####
## Usage

### Before usage

0.1 Download and setup all required software tools.

0.2 Download the reference files, and edit the pipeline files (“sccaller_pipeline_1.sh” and “sccaller_pipeline_2.sh”) to reflect the  directories and file names of the reference files in your system.

0.3 Edit the pipeline files (“sccaller_pipeline_1.sh” and “sccaller_pipeline_2.sh”) to use your job scheduler (As shown, it is an example to usge an SGE scheduler).

### Step 1. Quality control and alignment

This step uses 8 CPU cores and 32 GB RAM per sample, and takes 2-5 days (may flacturate more) depending on the performance of your computer cluster.

1.1 Deposit fastq files under folder for the pair-end reads for a specific sample:
```shell
./fastq/${sample_id}_1.fq # read 1
./fastq/${sample_id}_2.fq # read 2
```
1.2 Submit this job to a computer cluster as the following (shown for SGE):

For a single cell sample:
```shell
qsub sccaller_pipeline_1.sh ${sample_id} cell
```
For a bulk DNA sample: 
```shell
qsub sccaller_pipeline_1.sh ${sample_id} bulk
```

### Step 2. Variant calling and mutation burden estimation

2.1 Submit this job to a computer cluster as the following (shown for SGE):
```shell
qsub sccaller_pipeline_2.sh ${bulk_id} ${cell_id}
```
"bulk_id" and "cell_id" are the "sample_id" of the bulk DNA and the single cell, respectively.



