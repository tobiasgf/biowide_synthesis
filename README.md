# Multitaxon inventory reveals highly consistent biodiversity responses to ecospace variation - bioinformatic analyses  
___

This repository (R project) contains all data, scripts and links necessary to run the ecological/statistical analyses as well as the bioinformatic eDNA analyses from the study  **Multi-taxon inventory reveals highly consistent biodiversity responses to ecospace variation**  (Brunbjerg et al.)  

## Preparation of OTU tables / sequencing data
Links to the Illumina sequence data is referenced in the relevant scripts  

### most important CLI tools were used for this study  

 * VSEARCH v.2.9 (or later) (https://github.com/torognes/vsearch) 
 * Cutadapt v 1.17 (https://cutadapt.readthedocs.io/en/stable/)  
 * blastn v2.4.0+ (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) 
 
Various R-packages were used for this study (see in the relevant markdown files).  

### Description of the sub directories for the sequence data processing  

 * seq_data : contains all the data produced as part of the analyses here  
 * seq_processing : contains all the scripts and files necessary to perform the initial sequence processing  
 * R : contains a few functions used in the sequence data analyses  
 * seq_plots : output directory for the plots/figures produced in the sequence analyses

### Sequence of analyses
The preparation of the three different OTU tables used in this study can be carried out in any order as they are not dependent on each other. But they need to be carried out before doind the statistical analyses and modelling.

## Statistical analyses and modelling
Here AKB will insert an overall description of the scripts necessary to reprosuce the study. Something like:  

The ecological/statistical analyses and modelling employied previously published data and new data (including the sequencing data described above). All data and scripts necessary to reproduce the analyses have been included in this repository and a structured in the following way.  

  * eco_data: contains all the data produced as part of the analyses here, as well as copies/version of previously published data necessary for the analyses  
  * eco_processing: contains all  the scripts necessary to perform the analyses/modelling  
  * eco_plots: contains all figures produced in the ecological analyses  
  