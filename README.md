# SCcaller-pipeline
A computational pipeline for analyzing single-cell whole-genome sequencing data: from fastq files to mutation calls and mutation burdens.

Version: 1.0.0

Updated date: 2022.12.10

Citation:

Lei Zhang, Moonsook Lee, Alexander Y. Maslov, Cristina Montagna, Jan Vijg, and Xiao Dong. Single-cell whole-genome sequencing for discovering somatic mutations. In submission.

#####
## Author and License

Author: Xiao Dong

Email: dong0265@umn.edu (X.D.)

Licensed under the GNU Affero General Public License version 3 or later.

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

Current SCcaller-pipeline includes two major steps, and 17 sub-steps. See the flowchart, and usage below.

![alt text](https://github.com/XiaoDongLab/SCcaller-pipeline/blob/main/sccaller_pipeline_flowchart_v1.0.0.png)

#####
## Usage and Results

### Before usage

• 0.1 Download and setup all required software tools.

• 0.2 Download the reference files, and edit the pipeline files (“sccaller_pipeline_1.sh” and “sccaller_pipeline_2.sh”) to reflect the  directories and file names of the reference files in your system.

• 0.3 Edit the pipeline files (“sccaller_pipeline_1.sh” and “sccaller_pipeline_2.sh”) to use your job scheduler (As shown, it is an example to usge an SGE scheduler).

### Step 1. Quality control and alignment

#### Usage of step 1

This step uses 8 CPU cores and 32 GB RAM per sample (either cell or bulk), and takes 2-5 days (may flacturate more) depending on the performance of your computer cluster.

• 1.1 Deposit fastq files under folder for the pair-end reads for a specific sample:
```shell
./fastq/${sample_id}_1.fq.gz # read 1
./fastq/${sample_id}_2.fq.gz # read 2
```
• 1.2 Submit this job to a computer cluster as the following (shown for SGE):

For a single cell sample:
```shell
qsub sccaller_pipeline_1.sh ${sample_id} cell
```
For a bulk DNA sample: 
```shell
qsub sccaller_pipeline_1.sh ${sample_id} bulk
```
#### Understanding results of step 1

Result files of step 1 is explained in the table below.

|Directory|File(s)|Note and usage|
| --- | --- | --- |
|./${sample_id}|${sample_id}_1 (or 2) _fastqc.html (or .zip)|FastQC reports|
|./${sample_id}|${sample_id}.markilluminaadapters_metrics.txt|Markadaptor report|
|./${sample_id}|${sample_id}.bwa.stderr.log|Alignment report|
|./${sample_id}|${sample_id}.hg38.duplicate_metrics|Markduplicate report|
|./${sample_id}|${sample_id}.hg38.recal_data.csv|RBQS report|
|./${sample_id}|${sample_id}.hg38.gatk4.bam|Alignment bam file; input of the next step|
|./${sample_id}|${sample_id}.hg38.gatk4.bai|bam index|
|./${sample_id}|${sample_id}.hg38.gatk4.bam.md5|bam md5|
|./coverage|summary_coverage.csv|Summary of sequencing coverage|
|./ht|${sample_id}.gvcf.gz|Germline variants (vcf; only for bulk DNA)|
|./ht|${sample_id}.gvcf.gz.tbi|Index of germline variant vcf file (only for bulk DNA)|
|./ht|${sample_id}.hg38.hsnp.biallelic.dbsnp.vcf|Heterozygous germline SNVs reported in dbSNP (vcf; only for bulk DNA); input of the next step|
|./ht|${sample_id}.hg38.hindel.biallelic.dbsnp.vcf|Heterozygous germline INDELs reported in dbSNP (vcf; only for bulk DNA); input of the next step|
|./sensitivity|${sample_id}.hsnp.biallelic.dbsnp.20x.sorted.bed|Heterozygous germline SNVs reported in dbSNP (bed; only for bulk DNA); input of the next step|
|./sensitivity|${sample_id}.hindel.biallelic.dbsnp.20x.sorted.bed|Heterozygous germline INDELs reported in dbSNP (bed; only for bulk DNA); input of the next step|

#### Check points of step 1

• For cells of substantially less sequencing coverage (e.g., 30% genome covered with 20x depth), first check if enough sequencing data was obtained for both single cells and their bulk DNA, and if so, this suggests severe allelic dropouts of the cells, which are rare for amplicons of enough yield and passed the LDO test.

• For most bulk DNA, typically >1 million heterozygous germline SNVs and >100 thousand heterozygous germline INDELs are observed from its sequencing data and reported by the dbSNP database before.

### Step 2. Variant calling and mutation burden estimation

#### Usage of step 2

This step uses 8 CPU cores and 64 GB RAM per cell, and takes 1-2 days (may flacturate more) depending on the performance of your computer cluster.

• 2.1 Submit this job to a computer cluster as the following (shown for SGE):
```shell
qsub sccaller_pipeline_2.sh ${bulk_id} ${cell_id}
```
"bulk_id" and "cell_id" are the "sample_id" of the bulk DNA and the single cell, respectively.

#### Understanding results of step 2

Result files of step 2 is explained in the table below.

|Directory|File(s)|Note and usage|
| --- | --- | --- |
|./sccaller|${cell_id}.sccaller.vcf|SCcaller results for all candidate variant positions|
|./sccaller|${cell_id}.somatic.snv.vcf|SCcaller results for somatic SNVs|
|./sccaller|${cell_id}.somatic.indel.vcf|SCcaller results for somatic INDELs|
|./sccaller|sc_${cell_id}.sccaller_01to-1.log|SCcaller log file|
|./sensitivity|sensitivity.txt|Estimated variant calling sensitivity based on germline heterozygous SNVs and INDELs|
|./coverage|summary_coverage_both.csv|No. base pairs covered with at least 20x in both single cell and its bulk|
|./mutationburden|mutation_burden.csv|A summary of coverage, sensitivity, and SNV and INDEL burdens|

#### Check pionts of step 2

• Regarding variant calling quality, although we did not observe any cells of significantly low variant calling quality in all our previous projects, it is still important to monitor its quality. The η value estimated by SCcaller indicates bias in genome amplification and can be used as an important quality measurement: a lower η value corresponds to higher bias and results a more conservative cutoff of variant calling in SCcaller. In the figure below, we provide η values of all cells in three major datasets above including human lymphocytes (Zhang et al), hepatocytes (Brazhnik et al), and lung epithelium cells (Huang et al) for users to compare with their own data.

![alt text](https://github.com/XiaoDongLab/SCcaller-pipeline/blob/main/sccaller_pipeline_eta_v1.0.0.png)

In the figure above, average η value of autosomes per cell was calculated for each single cell in the three studies. Boxplot elements are defined as: center line indicates median, box limits indicate upper and lower quartiles, whiskers indicate 1.5× interquartile range, and points indicate the outliers of the boxplots.

#####
## Release Notes
• v1.0.0, 2022.12.10, 1st release version.
