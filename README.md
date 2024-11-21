# NGS_Analysis
A comprehensive pipeline performing Next-Generation Sequencing (NGS) data analysis.

Input
A folder containing sequence fastq files of Patients sufferening from different Hepatitis B patients undering REP2139 Therapy. 
The files are paired read fastq files (_R1.fastq and _R2.fastq) from patients at varying timepoints  (i.e 3-4 timepoints). A reference sequence file containing reference sequence of HBV downloaded from REFSEQ NCBI (https://www.ncbi.nlm.nih.gov/nuccore/NC_003977.2).
A Rscript to Plot Coverage distribution Plot.

Method
NGS Sequence analysis includes 5 key steps.
1- View Quality control using FASTQC 
2- Filtering reads of low quality bases using Trim galore, trimmomatic etc
3- Create an alignment using BWA , TopHat, etc
4- Create a consensus sequence using different time points sequences and generating one consensus sequence to visualize the mutations trajectory over the course of infection.



