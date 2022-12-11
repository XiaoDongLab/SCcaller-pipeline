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

![alt text](https://github.com/XiaoDongLab/SCcaller-pipeline/blob/main/sccaller_pipeline_flowchart_v1.0.0.png)

