# PORE-cupine
## First version of PORE-cupine. Detecting SHAPE modification using direct RNA sequencing

The codes in the two folders are similar.  
For single gene is suitable for running one gene at a time.  
For transcriptome is recommeded for running multiple genes.  
Both will yield the same results.  

### Programs needed to run the analysis:
Albacore (Oxford nanopore)  
Nanopolish (https://github.com/jts/nanopolish) a modified copy that removes the outliers from fast5 is included here  
Graphmap (https://github.com/isovic/graphmap)  
R (https://www.r-project.org/)  

### R packages required:
dplyr  
e1071  
data.table  
ggplot2  
ggpubr  
optparse  
Rcpp  


## Steps:

### To basecall raw fast5, output for both fast5 and fastq is required 
read_fast5_basecaller.py -i "location of fast5" -s "output_location" -r -k SQK-RNA001 -f FLO-MIN106 -o fast5,fastq --disable_filtering

### To map
cat fastq* | sed 's/U/T/g' > coverted.fastq
graphmap align -r "reference.fa" -d coverted.fastq -o gene.sam  --double-index
samtools view -bT "reference.fa" -F 16 gene.sam > gene.bam
samtools sort gene.bam > gene.s.bam
samtools index gene.s.bam

### aligning of raw signal with nanopolish
nanopolish index -d "location of basecalled fast5" converted.fastq
### scaling of events current to the model current is required
nanopolish eventalign  --reads converted.fastq --bam gene.s.bam --genome "reference.fa" --print-read-names --scale-events > gene.event

## For single genes
### To combine mulitple events from same position and strands
~/Read_events.R -f gene.event -o combined.RData

### To generate reacitvity
~/SVM.R -m "modified_gene.RData" -u "unmodified_gene.RData" -o "output file names.csv"
