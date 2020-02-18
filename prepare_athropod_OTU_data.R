# Analysis Biowide: preparing mt16S insect OTU data
# Manuscript: Multi-taxon inventory reveals highly consistent biodiversity responses to ecospace variation
# Author: Tobias Guldberg Frøslev
# Date: 30-01-2020


#### get data and make initial processing ####
#intermediate data from these steps not included in this GitHub repository

#download sequence data from: To be added upon acceptance for publication

#demultiplex the data using the tag-info files available in the seq_processing sub-directory "mt16S_insect" and the script: https://github.com/tobiasgf/man_vs_machine/blob/master/seq_processing/DADA2_demultiplex.sh

#analyse using DADA2 implemented in the script (seq_processing) : dada2_v3.1.R
#   R<dada2_v3.1.R --no-save &>log.txt
#(needs info from file variables.txt in dir seq_processing/arthropod_data)

#Rename and proceed with the dada2 table with chimeras removed:
#   mv DADA2_nochim.table art_DADA2_nochim.table

#Get taxonomic affiliation for all OTUs (doint it for the data where chimeras have not been removed.
#   bash spilt_and_blast.sh DADA2_raw.otus # You could just blast the NOCHIM version if you are sure that is the one you will be using.
#(resulting file has been zipped due to size)

#lulu matchlist
#format sequences to only contain sequence ids  
#   makeblastdb -in DADA2_nochim.centroids -dbtype nucl
#   blastn -db DADA2_nochim.centroids -num_threads 50 -outfmt '6 qseqid sseqid pident' -out art_nochim.matchlist -qcov_hsp_perc 80 -perc_identity 84 -query DADA2_nochim.centroids
#(matchlist has been zipped due to size)

#### environment ####
library(here)
library(lulu)
library(dplyr)
library(readr)
library(stringr)

#### lulu ####
#read otu tab and matchlist and perform lulu on the data
otutab <- read.csv(here::here("seq_processing/arthropod_data","DADA2_nochim.otutable"),sep='\t',header=T,row.names = 1, as.is=TRUE)
matchlist <-  read.csv(here::here("seq_processing/arthropod_data","art_nochim.matchlist"),sep='\t',header=F,as.is=TRUE)
lulified_tab <- lulu(otutab,matchlist, minimum_match = 84, minimum_relative_cooccurence = 1)
saveRDS(lulified_tab,here::here("data","lulified_art_nochim.RDS"))

#### dereplicate ####
source(here::here("R","dereplicate.r"))
#collapse the two replicates for each sample
otutab <- read.csv(here::here("seq_processing/arthropod_data","DADA2_nochim.otutable"),sep='\t',header=T,row.names = 1, as.is=TRUE)
samplelist <- read.table(here::here("data","sampleid_arthropod.txt"), header=TRUE, as.is = TRUE)
derep_tab <- dereplicate(otutab, samplelist)

saveRDS(derep_tab,here::here("data","derep_nochim_arthropod.RDS"))


#### assign taxonomy ####
source(here::here("R","unitax_lineage.R"))

options(ENTREZ_KEY = "PUT_YOUR_KEY_HERE")

IDtable=read.csv(here::here("data","DADA2_raw.centroids.blasthits"),sep='\t',header=F,as.is=TRUE)
names(IDtable) <- c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","qlen","qcovs","sgi","sseq","ssciname","staxid")

pf <- prefilter(IDtable)
gc <- get_classification(IDtable2=pf$preprocessed_table)
cf <- evaluate_classification(classified=gc$taxonomic_info)
ac <- adjust_classification(class_table=cf$taxonon_table)
result <- list(adjusted_classified_table = ac$adjusted_classification, classified_table=cf$taxonon_table, all_classifications=cf$all_taxa_table, all_classifications_summed=cf$all_taxa_table_summed)

saveRDS(result, here::here("data","myclassifiedresults_arthropods.RDS"))
#write.table(result$adjusted_classified_table, here::here("data","my_classified_otus_adjusted.txt"), sep = "\t", quote = F, row.names = F)

#### focus tables on arthropod otus ####
lulified_tab <- readRDS(here::here("data","lulified_art_nochim.RDS"))
derep_tab <- readRDS(here::here("data","derep_nochim_arthropod.RDS"))
taxonomy <- readRDS(here::here("data","myclassifiedresults_arthropods.RDS"))
tax_tab <- taxonomy$adjusted_classified_table
tax_tab$OTU_ID <- tax_tab$qseqid
  
tax_tab2 <- tax_tab %>% group_by(species) %>% mutate(redundancy = n()) # add count of each species name

lulified_to_keep <- tax_tab2$OTU_ID[tax_tab2$OTU_ID %in% lulified_tab$discarded_otus & tax_tab2$redundancy == 1 & tax_tab2$pident == 100] # which lulu-flagged otus should be kept (high match, non-redundant)
keep_curated <- c(lulified_tab$curated_otus,lulified_to_keep) # which flagged errors should be removed
tax_tab3 <- tax_tab2[tax_tab2$OTU_ID %in% keep_curated,c(6,9,10,11,12,13,14,15,23,24,25)] # focus taxonomy on otus to keep

arthropod_tax <- tax_tab3 %>% filter(phylum == "Arthropoda")

otu_tab <- derep_tab[arthropod_tax$OTU_ID,] # focus otu_tab on otus to keep
otu_tab$OTU_ID <- row.names(otu_tab)

edna_tab <- left_join(otu_tab, arthropod_tax, by = "OTU_ID") # add taxonomy

#table with sampling rounds separate
write.table(edna_tab, here::here("data","arthropod_table_final.txt"), sep = "\t", quote = F, row.names = F)
