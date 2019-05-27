# Multitaxon inventory reveals highly consistent biodiversity responses to ecospace variation - bioinformatic analyses  
___

This repository (R project) contains all data, scripts and links necessary to run the bioinformatic eDNA analyses from the study  **Multitaxon inventory reveals highly consistent biodiversity responses to ecospace variation - bioinformatic analyses**  (Brunbjer et al. (in prep).  
All steps/processes for this study can be carried out on the same computer/platform. But, in practise all analyses were carried out on a linux server setup with 64 processors (AMD Opteron(tm) 6380), except R-scripts, which were run on a MacBook Pro (2.9 GHz Intel Core i7, 16 GB 2133 MHz LPDDR3).
All analyses were carried out in one directory and sub-directories of this.

## Sequencing data
Links to the Illumina sequence data is referenced in the relevant scripts  

## Bioinformatic tools
### most important CLI tools were used for this study  

 * VSEARCH v.2.9 (or later) (https://github.com/torognes/vsearch) 
 * Cutadapt v 1.17 (https://cutadapt.readthedocs.io/en/stable/)  
 * blastn v2.4.0+ (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) 
 
Various R-packages were used for this study (see in the relevant markdown files).  

## Description of the sub directories  

 * data : contains all the data produced as part of the analyses here  
 * seq_processing : contains all the scripts and files necessary to perform the initial sequence processing  
 * R : contains a few functions used in the analyses  
 * plots : output directory for the plots/figures produced in the analyses
