# The scripts needs the following packages
library(taxize)
library(dplyr)
library(tidyr)

#  IDtable=read.csv("BLAST_HIT_OUTPUT",sep='\t',header=F,as.is=TRUE)
# names(IDtable) <- c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","qlen","qcovs","sgi","sseq","ssciname","staxid")

# take a look at the data and save it
#     my_clasified_result$classified_table
#     my_clasified_result$all_classifications # Here you can review all scores for all matches per OTU
#     my_clasified_result$all_classifications_summed # Here you can review all summed scores for all taxa per OTU
#     write.table(my_clasified_result$classified_table, "my_classified_otus.txt", sep = "\t", quote = F, row.names = F)
#     write.table(my_clasified_result$adjusted_classified_table, "my_classified_otus_adjusted.txt", sep = "\t", quote = F, row.names = F)
#     saveRDS(my_clasified_result,"myclassifiedresultsRDS")

#pf <- prefilter(IDtable)
#gc <- get_classification(IDtable2=pf$preprocessed_table)
#cf <- evaluate_classification(classified=gc$taxonomic_info)
#ac <- adjust_classification(class_table=cf$taxonon_table)
#result <- list(adjusted_classified_table = ac$adjusted_classification, classified_table=cf$taxonon_table, all_classifications=cf$all_taxa_table, all_classifications_summed=cf$all_taxa_table_summed, upper_margin=upper_margin, lower_margin=lower_margin, remove=remove, useDB=useDB, appendDB=appendDB, db_path=db_path, tax_db=tax_db, id_cut=id_cut, ambiguity_cutoff=ambiguity_cutoff, used_from_database=gc$used_from_database, added_to_database=gc$added_to_database)

#write.table(result$adjusted_classified_table, "my_classified_otus_adjusted.txt", sep = "\t", quote = F, row.names = F)
#saveRDS(result,"myclassifiedresultsRDS")

#table(result$adjusted_classified_table$order)

#Function2
#Filter data OTU wise according to upper and lower margin set, and taxa to exclude
prefilter <- function(IDtable, upper_margin=5, lower_margin=10, remove = c("uncultured","environmental","N/A"), min_cov = 50){
 print(paste0("Collapsing taxids and preprocessing")) # make a progressline (indicating the index the loops needs to be
 new_IDtable <- IDtable[0,] # prepare filtered matchlist
 IDtable <- IDtable[!IDtable$staxid == "N/A",]
 IDtable <- IDtable[which(IDtable$qcovs >= min_cov),]
 ids <- names(table(IDtable$qseqid))
 i=1
 o=length(ids)
 pb <- txtProgressBar(min = 0, max = o, style = 3)
 for (name in ids){
  test <- IDtable[which(IDtable$qseqid == name),] # select all lines for a query
  if (nchar(remove[1])>0){
   test2 <- test
   for (rm in 1:length(remove)){
    test2 <- test2[!grepl(remove[rm], test2$ssciname,ignore.case = TRUE),]
   }
   if (nrow(test2) > 1) {test <- test2}
  }
  max <- max(test$bitscore)
  upper <- max-upper_margin
  lower <- max-lower_margin
  test <- test[which(test$bitscore >= lower),] # select all lines for a query
  test$margin <- "lower"
  test[test$bitscore >= upper,"margin"] <- "upper"
  
  new_IDtable = rbind(new_IDtable,test) # add this row to the filtered IDtable
  i=i+1
  setTxtProgressBar(pb, i)
 }
 close(pb)
 result <- list(preprocessed_table = new_IDtable, upper_margin=upper_margin, lower_margin=lower_margin, remove=remove)
 return(result)
}

#Function3
# Get full taxonomic path for all hits within the upper limit of each OTU. Identical species are only queried once....
get_classification <- function(IDtable2, lineages = "~/tax_db/lineages-2019-02-20.csv", new_lineages="~/tax_db/extraTaxidsRDS", tax_levels=c("kingdom","phylum","class","order","family","genus","species")){
 
 require(taxize)
 reused_taxids <- vector()
 classified_taxids <- vector()
 all_staxids <- names(table(IDtable2$staxid[IDtable2$margin=="upper"])) # get all taxids for table
 all_classifications <- list() # prepare list for taxize output
 
 if(!file_test("-f", lineages)){
  stop(paste0("No lineage file (", lineages,") available, download latest lineage file from https://github.com/zyxue/ncbitax2lin"))
 }
 require(data.table)
 taxdbncbi <- as.data.frame(fread(lineages, header = T, sep = ','))
 lineage_taxids <- taxdbncbi$tax_id[taxdbncbi$tax_id %in% all_staxids]
 unknown_taxids <- all_staxids[!all_staxids %in% taxdbncbi$tax_id]
 print(paste0("Found ", length(lineage_taxids), " of ", length(all_staxids)," taxids in ",lineages))
 grapped_lineages <- taxdbncbi[taxdbncbi$tax_id %in% all_staxids, tax_levels]
 grapped_lineages$staxid <- taxdbncbi$tax_id[taxdbncbi$tax_id %in% all_staxids]
 if(!file_test("-f", new_lineages)){
  print(paste0("No extra database file (", new_lineages,") available, classifying all ", length(all_staxids)-length(lineage_taxids), " taxids not in ",lineages," from scratch"))
 } else {
  storedTaxids <- readRDS(new_lineages)
  stored_vector <- storedTaxids$taxids %in% unknown_taxids
  classify_vector <- !unknown_taxids %in% storedTaxids$taxids
  print(paste0("Using ",sum(stored_vector)," stored extra taxids from database"))
  reused_taxids <- storedTaxids$taxids[stored_vector]
  reused_classifications <- storedTaxids$classifications[stored_vector]
  classify_taxids <- unknown_taxids[classify_vector]
 }
 
 o=length(classify_taxids) # number of taxids
 
 Start_from <- 1 # change if loop needs to be restarted due to time-out
 new_classifications <- list()
 #Get ncbi classification of each entry
 if(o>0){
  print(paste0("Fetching classifications for ",o," taxids from scratch"))
  pb <- txtProgressBar(min = 0, max = o, style = 3)
  for (cl in Start_from:o){ # the taxize command "classification" can be run on the all_staxids vector in one line, but often there is
   #a timeout command, therefor this loop workaround.
   #restarted from if it quits)
   new_classifications[cl] <- classification(classify_taxids[cl], db = "ncbi")
   if(round(cl/10) == cl/10){Sys.sleep(2)}
   setTxtProgressBar(pb, cl)
   
   if(!file_test("-f", new_lineages)){
    print(paste0("No database file found. Saving all ",length(classify_taxids)," classified taxids in a new database ",new_lineages,")"))
    append_tax <- list(taxids=all_staxids[cl],classifications=all_classifications[cl])
    saveRDS(append_tax,new_lineages)
   } else {
    storedTaxids <- readRDS(new_lineages)
    append_tax <- list(taxids=c(storedTaxids$taxids,classify_taxids[cl]),classifications=c(storedTaxids$classifications,new_classifications[cl]))
    saveRDS(append_tax,new_lineages)
   }
  }
  close(pb)
  print(paste0("Appended classification for ",length(classify_taxids)," taxids to extra database"))
 }
 
 unknown_staxids <- c(classify_taxids,reused_taxids)
 unknown_classifications <- c(new_classifications,reused_classifications)
 if(length(classify_taxids)>0){
 #Construct a taxonomic path from each classification
 output <- data.frame(matrix(ncol = 1+length(tax_levels), nrow = length(unknown_staxids)))
 x <- c("staxid", tax_levels)
 colnames(output) <- x
 
 totalnames <- length(unknown_staxids)
 pb <- txtProgressBar(min = 0, max = totalnames, style = 3)
 for (curpart in seq(1:totalnames)){
  currenttaxon <- unknown_classifications[curpart][[1]]
  if (nchar(currenttaxon[1]) > 0 & !is.na(currenttaxon[[1]][1])) {
   spec <- unknown_staxids[curpart]
   for (lev in seq(1:length(tax_levels))){
    output[curpart,tax_levels[lev]] <- currenttaxon[which(currenttaxon$rank == tax_levels[lev]),"name"][1]
    #output[curpart,"phylum"] <- currenttaxon[which(currenttaxon$rank == "phylum"),"name"][1]
    #output[curpart,"class"] <- currenttaxon[which(currenttaxon$rank == "class"),"name"][1]
    #output[curpart,"order"] <- currenttaxon[which(currenttaxon$rank == "order"),"name"][1]
    #output[curpart,"family"] <- currenttaxon[which(currenttaxon$rank == "family"),"name"][1]
    #output[curpart,"genus"] <- currenttaxon[which(currenttaxon$rank == "genus"),"name"][1]
    #output[curpart,"species"] <- currenttaxon[which(currenttaxon$rank == "species"),"name"][1]
   }
   output[curpart,"staxid"] <-  spec # add that row to the filtered IDtable
  }
  setTxtProgressBar(pb, curpart)
 }
 close(pb)
 }
 
 if(length(classify_taxids)>0){
 output2 <- rbind(output,grapped_lineages)
 } else {
  output2 <- grapped_lineages
 }
 
 taxonomic_info <- merge(IDtable2,output2,by = "staxid", all=TRUE)
 taxonomic_info$species[is.na(taxonomic_info$species)] <- taxonomic_info$ssciname[is.na(taxonomic_info$species)]
 result <- list(taxonomic_info=taxonomic_info)
 return(result)
}

#Function4
#Function for evaluating the taxonomic assignment of each OTU. All hits within the upper margin are used in the evaluation weithted by thei evalue, so that suboptimal matches has a lower weight. All hits within the lower margin are put into the output (but not used for evaluating classification)
evaluate_classification <- function(classified, remove=c("uncultured","environmental","N/A")){
 require(tidyr)
 require(dplyr)
 ids <- names(table(classified$qseqid))
 i <- 1
 print("Evaluating the classifications")
 pb <- txtProgressBar(min = 0, max = length(ids), style = 3)
 for (name in ids){
  setTxtProgressBar(pb, i)
  test <- classified[which(classified$qseqid == name),]
  test2 <- test %>% filter(margin == "upper")
  if (nchar(remove[1])>0){
   test2x <- test2
   for (rm in 1:length(remove)){
    test2x <- test2x[!grepl(remove[rm], test2x$species,ignore.case = TRUE),]
   }
   if (nrow(test2x) > 1) {test2 <- test2x}
  }
  
  test2$score <- 100*(1/(test2$evalue+2e-300))/sum(1/(test2$evalue+2e-300))  # HER BEREGSES SCOREN FOR ALLE MATCHES PER OTU
  #test2 %>% group_by(species) %>% mutate(n_obs=n())
  
  test4 <- test2 %>% filter(margin == "upper") %>%
   dplyr::select(margin,qseqid,sgi,sseq,staxid,pident,score,qcovs,kingdom,phylum,class,order,family,genus,species) %>% 
   group_by(qseqid,kingdom, phylum,class,order,family,genus,species) %>% 
   mutate(species_score=sum(score)) %>% 
   group_by(qseqid,kingdom, phylum,class,order,family,genus) %>% 
   mutate(genus_score=sum(score)) %>%
   group_by(qseqid,kingdom, phylum,class,order,family) %>% 
   mutate(family_score=sum(score))%>%
   group_by(qseqid,kingdom, phylum,class,order) %>% 
   mutate(order_score=sum(score)) %>%
   group_by(qseqid,kingdom, phylum,class) %>% 
   mutate(class_score=sum(score)) %>%
   group_by(qseqid,kingdom, phylum) %>% 
   mutate(phylum_score=sum(score)) %>%
   group_by(qseqid,kingdom) %>% 
   mutate(kingdom_score=sum(score)) %>% ungroup() %>%
   arrange(-kingdom_score,-phylum_score,-class_score,-order_score,-family_score,-genus_score,-species_score)
  test3 <- test4 %>% slice(1)
  test5 <- test4 %>% distinct(qseqid,sgi,sseq,pident,qcovs,kingdom, phylum,class,order,family,genus,species,kingdom_score,phylum_score,class_score,order_score,family_score,genus_score,species_score) 
  string1 <- test %>% dplyr::group_by(species,pident) %>%  summarize(count=n()) %>% dplyr::select(species,count,pident) %>% arrange(-pident,-count) %>% t()
  string2 <- toString(unlist(string1))
  test3$alternatives <- string2
  if (i == 1){result <- test3} else{
   result <- rbind(result,test3)
  }
  if (i == 1){result2 <- test2} else{
   result2 <- rbind(result2,test2)
  }
  if (i == 1){result3 <- test5} else{
   result3 <- rbind(result3,test5)
  }
  i=i+1
 }
 close(pb)
 total_result <- list(taxonon_table = result, all_taxa_table=result2, all_taxa_table_summed=result3, remove=remove)
 return(total_result)
}

# Function for adjusting taxonomic annotation. Annotation is adjusted based on the level of match (id_cut). Default settings, applicable for ITS2 data, 
# is 98, 90, 85, 80, 75, 70, 50 for species, genus, family, order, class, phylum, kingdom assignment.
# Assignment is also adjusted for taxonomic agreement among reference database matches. Default threshod scores are 90 for accepting a match at any taxoomic level.

adjust_classification <- function(class_table, id_cut=c(98, 90, 85, 80, 75, 70, 50), ambiguity_cutoff=c(50,50,50,50,50,50,50)){
 
 my_classifications <- class_table
 
 cutoff_index <- which(my_classifications$pident < id_cut[7])
 reclassify_index <- c("kingdom","phylum","class","order","family","genus")
 my_classifications[cutoff_index,reclassify_index] <- "unmatched"
 my_classifications$species[cutoff_index] <- paste0("unmatched_",my_classifications$phylum[cutoff_index], "_sp")
 
 cutoff_index <- which(my_classifications$pident < id_cut[6]  & my_classifications$pident >= id_cut[7] & my_classifications$kingdom != "unmatched")
 reclassify_index <- c("phylum","class","order","family","genus")
 my_classifications[cutoff_index,reclassify_index] <- "unmatched"
 my_classifications$species[cutoff_index] <- paste0("unmatched_",my_classifications$phylum[cutoff_index], "_sp")
 
 cutoff_index <- which(my_classifications$pident < id_cut[5]  & my_classifications$pident >= id_cut[6] & my_classifications$phylum != "unmatched")
 reclassify_index <- c("class","order","family","genus")
 my_classifications[cutoff_index,reclassify_index] <- "unmatched"
 my_classifications$species[cutoff_index] <- paste0("unmatched_",my_classifications$phylum[cutoff_index], "_sp")
 
 cutoff_index <- which(my_classifications$pident < id_cut[4] & my_classifications$pident >= id_cut[5] & my_classifications$class != "unmatched")
 reclassify_index <- c("order","family","genus")
 my_classifications[cutoff_index,reclassify_index] <- "unmatched"
 my_classifications$species[cutoff_index] <- paste0("unmatched_",my_classifications$class[cutoff_index], "_sp")
 
 cutoff_index <- which(my_classifications$pident < id_cut[3] & my_classifications$pident >= id_cut[4] & my_classifications$order != "unmatched")
 reclassify_index <- c("family","genus")
 my_classifications[cutoff_index,reclassify_index] <- "unmatched"
 my_classifications$species[cutoff_index] <- paste0("unmatched_",my_classifications$order[cutoff_index],"_sp")
 
 cutoff_index <- which(my_classifications$pident < id_cut[2] & my_classifications$pident >= id_cut[3] & my_classifications$family != "unmatched")
 reclassify_index <- c("genus")
 my_classifications$genus[cutoff_index] <- "unmatched"
 my_classifications$species[cutoff_index] <- paste0("unmatched_",my_classifications$family[cutoff_index], "_sp")
 
 cutoff_index <- which(my_classifications$pident < id_cut[1] & my_classifications$pident >= id_cut[2] & my_classifications$genus != "unmatched")
 my_classifications$species[cutoff_index] <- paste0("unmatched_",my_classifications$genus[cutoff_index], "_sp")
 
 my_classifications$kingdom[my_classifications$kingdom_score<ambiguity_cutoff[7]] <- "ambiguous"
 my_classifications$phylum[my_classifications$phylum_score<ambiguity_cutoff[6]] <- "ambiguous"
 my_classifications$class[my_classifications$class_score<ambiguity_cutoff[5]] <- "ambiguous"
 my_classifications$order[my_classifications$order_score<ambiguity_cutoff[4]] <- "ambiguous"
 my_classifications$family[my_classifications$family_score<ambiguity_cutoff[3]] <- "ambiguous"
 my_classifications$genus[my_classifications$genus_score<ambiguity_cutoff[2]] <- "ambiguous"
 
 print("Adjusting the classifications")
 pb <- txtProgressBar(min = 0, max = nrow(my_classifications), style = 3)
 
 for(i in 1:nrow(my_classifications)){
  setTxtProgressBar(pb, i)
  if (my_classifications[i,"pident"] > id_cut[1] & my_classifications[i,"species_score"] < ambiguity_cutoff[1]){
   if (!my_classifications[i,"genus"] %in% c("ambiguous", "unmatched", NA)){
    my_classifications[i,"species"] <- paste0("ambiguous_",my_classifications[i,"genus"],"_sp")
   } else if (!my_classifications[i,"family"] %in% c("ambiguous", "unmatched", NA)){
    my_classifications[i,"species"] <- paste0("ambiguous_",my_classifications[i,"family"],"_sp")
   } else if (!my_classifications[i,"order"] %in% c("ambiguous", "unmatched", NA)){
    my_classifications[i,"species"] <- paste0("ambiguous_",my_classifications[i,"order"],"_sp")
   } else if (!my_classifications[i,"class"] %in% c("ambiguous", "unmatched", NA)){
    my_classifications[i,"species"] <- paste0("ambiguous_",my_classifications[i,"class"],"_sp")
   } else if (!my_classifications[i,"phylum"] %in% c("ambiguous", "unmatched", NA)){
    my_classifications[i,"species"] <- paste0("ambiguous_",my_classifications[i,"phylum"],"_sp")
   } else if (!my_classifications[i,"kingdom"] %in% c("ambiguous", "unmatched", NA)){
    my_classifications[i,"species"] <- paste0("ambiguous_",my_classifications[i,"kingdom"],"_sp")
   }
  }
 }
 
 close(pb)
 result <-  list(adjusted_classification = my_classifications, id_cut=id_cut, ambiguity_cutoff=ambiguity_cutoff)
 return(result)
}
