# RNASeq Differential Expression Analysis.
## Project description
Comparative transcriptome profile of genes differentially expressed in longissimus dorsi muscles between Japanese black (Wagyu) and Chinese Red Steppes cattle by RNA-seq
## Data Description.
RNAseq samples from two cattle breeds i.e, Japanese black (Wagyu) cattle and Chinese Red Steppes cattle.

Sample      Breed
SRR13107018	Japanese black (Wagyu) cattle
SRR13107019	Japanese black (Wagyu) cattle
SRR13107020	Japanese black (Wagyu) cattle
SRR13107021	Chinese Red Steppes cattle
SRR13107022	Chinese Red Steppes cattle
SRR13107023	Chinese Red Steppes cattle

*** For practise, you can obtain the fastq files from NCBI SRA website using the `fastq-dump`command from the SRA Toolkit to download the data in FASTQ format and store it in the Data folder.

## RNASeq anaylsis workflow.

### The RNAexpression.sh script.
This script is a pipeline for RNA-Seq data analysis, and it includes quality control, alignment to a reference genome, transcript assembly, and abundance estimation. 

#### The RNAexpression.sh Script Input:
Fasta files
#### The RNAexpression.sh Script output: 
The abundance estimation files for each sample will be stored in the ballgown folder within the Data folder. Theses will serve as an input for the R script

### The RNASeq.R Script.
This script is used to analyze RNA-Seq data for differential expression between different breeds of cattle. It filters the data, performs statistical tests, and then creates visualizations (plots and heatmaps) for the results. 

#### The RNASeq.R Script Input: 
ballgown folder abundance estimation files.
#### The RNASeq.R Script output:
Results plots folder containing Visualization plots heatmaps.