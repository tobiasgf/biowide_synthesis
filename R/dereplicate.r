dereplicate <- function(otutable,samplelist,exclude_samples = NULL,exclude_sites = NULL){
 #remove everything but Sxxx from headers
 col_names <- names(otutable)
 tax_cols <-  which(substring(col_names,1,1) != "S")
 sam_colno <-  which(substring(col_names,1,1)=="S")
 names(otutable)[sam_colno] <- substring(names(otutable)[sam_colno], 1, 4)
 sam_cols = col_names[which(substring(col_names,1,1)=="S")]
 
 #exdluding sample numbers that does not appear in the sequencing results (OTU table) from the samplelist
 samplelist_samples <- samplelist[,1]
 row.names(samplelist)=samplelist[,1]
 samplelist=samplelist[sam_cols,]
 samplelist_only <- setdiff(samplelist_samples,sam_cols)
 otutable_only <- setdiff(sam_cols,samplelist_samples)
 
 #excluding sample numbers that belongs to blanks/negatives from the samplelist
 for (ex in exclude_sites){
  ex_len <- nchar(ex)
  samplelist=samplelist[which(substring(samplelist$Site,1,ex_len) != ex),]
 }
 
 #excluding sample numbers that belongs to blanks/negatives from the samplelist
 for (ex in exclude_samples){
  ex_len <- nchar(ex)
  samplelist=samplelist[which(substring(samplelist$Sample,1,ex_len) != ex),]
 }
 
 #get a non-redundant list of site (biological sample) names
 site_names <- names(table(samplelist[,2]))
 sort_order <- order(substring(site_names, 3))
 site_names <- site_names[sort_order]
 
 otubuild = otutable[,0] # a rebuild OTU-table (with sitenames as headers insted of samplenames) but keeping replicates (samples) separate
 otuPA=otutable[,0] # A total presence/absence table
 OTU_sample_rep_count=otutable[,0] # OTU table (sites) with information about number of positive replicates (readcount > x) pr site.
 OTU_totalcount=otutable[,0] # a dereplicated OTU table (sites) with total readcounts from all replicates (samples)
 for (site in site_names){ # taking the sites one by one
  cursite <- samplelist[which(samplelist$Site == site),"Sample"]    # get the sampleIDs that belong to the current site  
  cursamples <- as.data.frame(otutable[,cursite]) # get the OTU colums matching sampleIDs
  colnames(cursamples) <- rep(site,dim(cursamples)[2])  # set the column names to siteID
  cur_pres_abs <- cursamples  # make a pres/abs version of the current select
  cur_pres_abs[cur_pres_abs>0]=1 # ----
  PA_count = rowSums(cur_pres_abs)  # count the number of presences pr OTU for the current site
  cur_total = rowSums(cursamples)  # find the total read count for the OTUs of the current site 
  otuPA=cbind(otuPA,cur_pres_abs)  # build current info onto otuPA table 
  otubuild=cbind(otubuild,cursamples)  # build current info onto otubuild table
  OTU_sample_rep_count=cbind(OTU_sample_rep_count,PA_count) # # build current info onto OTU_sample_rep_count
  OTU_totalcount=cbind(OTU_totalcount,cur_total) # build current info onto OTU_totalcount
 }
 names(OTU_totalcount) = site_names
 return(OTU_totalcount)
}
