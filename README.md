# PORE-cupine
## Chemical utilized probing interrogated using nanopores (PORE-cupine)
 Detecting SHAPE modification using direct RNA sequencing

For single gene is suitable for running one gene at a time.  
For transcriptome is recommeded for running multiple genes.  
Both will yield the same results.  

### Programs needed to run the analysis:
Albacore (Oxford nanopore)  
Nanopolish (https://github.com/jts/nanopolish) a modified copy that removes the outliers from fast5 is included (nanopolish-edited.zip)   
Graphmap (https://github.com/isovic/graphmap)  
R (https://www.r-project.org/)  

### R packages required:
dplyr  
e1071  
data.table  
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
./Read_events.R -f gene.event -o combined.RData

### To generate reacitvity
./SVM.R -m "modified_gene.RData" -u "unmodified_gene.RData" -o "output file names.csv"

## For transcriptome

### split events to individual transcript
./split_events.sh "folder to store tmp files" gene.event

### Optional step run if needed to combine tmp files from multiple flowcells
#### combined tmp files will be found in folder named combined
./combine.sh 

### To combine mulitple events from same position and strands
./loop_for_Read_files.sh "number of parts" "input folder" "output folder"

### To generate reactivity profile for mulitple transcript
./SVM_multi.R -s "number of parts" "RData folder containing modified samples" "RData folder containing unmodified samples" "Output folder"

# Acknowledgments
Li Chenhao for his help in getting me started and the calculation of error per strands  (https://github.com/lch14forever)  
Shen Yang for his code for aligning transcipt positions to genomic position and for the TRipseq anaylsis (https://github.com/shenyang1981)
Zhang Yu for the calucation of error rates.

For combining of standard devations with mean, standard devations and number of samples.
Headrick, T. C. (2010). Statistical Simulation: Power Method Polynomials and other Transformations. Boca Raton, FL: Chapman & Hall/CRC.
