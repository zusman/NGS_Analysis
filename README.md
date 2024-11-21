# NGS_Analysis
A comprehensive pipeline performing Next-Generation Sequencing (NGS) data analysis.

**Input**
A folder containing sequence fastq files of Patients sufferening from different Hepatitis B patients undering REP2139 Therapy. 
The files are paired read fastq files (_R1.fastq and _R2.fastq) from patients at varying timepoints  (i.e 3-4 timepoints). A reference sequence file containing reference sequence of HBV downloaded from REFSEQ NCBI (https://www.ncbi.nlm.nih.gov/nuccore/NC_003977.2).
A Rscript to Plot Coverage distribution Plot.

**Method**
NGS Sequence analysis includes 4 key steps.
1. View Quality control using FASTQC 
2. Filtering reads of low quality bases(Phred=33, Quality=20) using Trim galore, trimmomatic etc
3. Create an alignment using BWA , TopHat, etc. Run Picard for sorting, indexing Sam and Bam files.
4. Create a consensus sequence using different time points sequences and generating one consensus sequence to visualize the mutations trajectory over the course of infection.

**Results**
1. Filtered Paired Fastq files
2. Sam and Bam Mapped read File.
3. Nucleotide and Aminoacid consensus sequence for each patient across different time points.
4. Coverage plot of mapped read file for every sample across different time points.
   
