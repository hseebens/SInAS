
####### SInAS workflow: Integration and standardisation of alien species data ###########
##
## Step 3: Merging of standardised databases of alien species data
##
## Manuela Gómez-Suárez, Hanno Seebens, Giessen, 25.09.2025
#########################################################################################


MergeDatabases <- function(FileInfo=NULL, 
                           version=NULL,
                           outputfilename=NULL,
                           output=NULL){
  
  ## identify input datasets based on file name "StandardSpec_....csv"
  allfiles <- list.files(file.path("Output","Intermediate"))
  inputfiles_all <- allfiles[grep("Step5_StandardIntroYear_",allfiles)]
  inputfiles <- vector()
  for (i in 1:nrow(FileInfo)){
    # inputfiles <- c(inputfiles,grep(FileInfo[i,"Dataset_brief_name"],inputfiles_all,value=T))
    inputfiles <- c(inputfiles,paste("Step5_StandardIntroYear_",FileInfo[i,"Dataset_brief_name"],".csv",sep=""))
  }
  inputfiles <- inputfiles[!is.na(inputfiles)]
  
  ## merge columns #######################################
  
  for (i in 1:length(inputfiles)){#
    
    dat <- read.table(file.path("Output","Intermediate",paste0(inputfiles[i])),header=T,stringsAsFactors = F)

    cnames <- colnames(dat)
    cnames <- cnames[!cnames%in%c("verbatimTaxonRank","location_orig","Kingdom_user","Country_ISO","ISO3","eventDate_orig","Taxon_group")]
    dat <- dat[,colnames(dat)%in%cnames]

    dat$datasetName <- FileInfo[i,1]
    
    # remove records whose location could not be standardized 
    dat <- dat[!is.na(dat$locationID), ]
 
    # check each source before merging and remove records which are flagged as absent, 
    # if a present record for that taxon and location is available
    
    if (any(grepl("occurrenceStatus", cnames))) { # check availability of occurrenceStatus information
      
      # create abs_records twice (for alien and for native) to keep only "absent" records (for later removal) 
      #if a taxon is reported as both present and absent in the same location (no conflict with other DBs)
      
      # for introduced records
      abs_int_records <- dat %>%
        filter(establishmentMeans %in% c("uncertain", "introduced")) %>% 
        group_by(location, locationID, taxonID, taxon) %>%
        filter(all(c("absent", "present") %in% occurrenceStatus)) %>% # keep records where a taxon in a location is flagged as both absent and present
        filter(occurrenceStatus == "absent") %>% # keep absent record to remove it in a later step
        ungroup()
      
      # for native records 
      abs_nat_records <- dat %>%
        filter(establishmentMeans %in% "native") %>% 
        group_by(location, locationID, taxonID, taxon) %>%
        filter(all(c("absent", "present") %in% occurrenceStatus)) %>% # keep records where a taxon in a location is flagged as both absent and present
        filter(occurrenceStatus == "absent") %>% # keep absent record to remove it in a later step
        ungroup()
      
      # remove the absent record if there is already a present record for a taxon in a location
      dat <- dat %>%
        anti_join(abs_int_records, by = c("location", "locationID", "taxonID", "taxon", "establishmentMeans", "occurrenceStatus"))
      dat <- dat %>%
        anti_join(abs_nat_records, by = c("location", "locationID", "taxonID", "taxon", "establishmentMeans", "occurrenceStatus"))
    }
    
    if (i == 1) {
      
      alldat <- dat
      
    } else {
      
      ## Merge databases #####################
      
      alldat <- merge(alldat, dat, by = c("location", "locationID", "taxon", "scientificName", "taxonID","establishmentMeans"), all = TRUE)
      
      ## Check for and resolve duplicate columns
      if ("taxonID.x" %in% colnames(dat) && "taxonID.y" %in% colnames(dat)) {
        
        if (all(dat$taxonID.x == dat$taxonID.y, na.rm = TRUE)) { # Confirm the columns are identical
          dat$taxonID <- dat$taxonID.x  # consolidate columns
          dat <- dat[, !colnames(dat) %in% c("taxonID.x", "taxonID.y")]  # Drop now redundant columns
        } else {
          stop("The columns 'taxonID.x' and 'taxonID.y' are not identical.")
        }
      }
      
      alldat <- alldat %>%
        mutate(across(starts_with("eventDate"), as.character))
      
      ## treat duplicated columns (named by R default ".x" and ".y")
      if (any(grepl("\\.y",colnames(alldat)))){
        
        ## solve multiple first records
        if (any(colnames(alldat)=="eventDate.x")){ 
          
          ## select the minimum of multiple first records
          alldat$eventDate <- apply(alldat[,c("eventDate.x","eventDate.y")],1,function(s) ifelse(all(is.na(s)),NA,min(s,na.rm=T)))
          alldat <- alldat[,!colnames(alldat)%in%c("eventDate.x","eventDate.y")]
        }
      
        while(any(grepl("\\.y",colnames(alldat)))){ # check if discrepancy still exists
          colname_dupl <- colnames(alldat)[grep("\\.y",colnames(alldat))][1] # identify duplicated column
          colname_dupl <- gsub("\\..+$","",colname_dupl) # create new column new
          colnames_dupl <- colnames(alldat)[grep(colname_dupl,colnames(alldat))] # identify .x and .y component
          
          ## merge both columns into a new one by combining the content separated by ';'
          eval(parse(text=paste0("alldat$",substitute(a,list(a=colname_dupl)),"<-paste(alldat$",substitute(a,
                        list(a=colname_dupl)),".x,alldat$",substitute(a,list(a=colname_dupl)),".y,sep=\"; \")"))) # add column with database information
          alldat <- alldat[,!colnames(alldat)%in%colnames_dupl] # remove .x and .y columns

        }
      }
    }
    
   ## merge records flagged as both: introduced and uncertain
    alldat <- alldat %>%
      mutate(dist_group = case_when( # create establishmentMeans groups for merging
        establishmentMeans == "introduced" ~ "a",
        establishmentMeans == "uncertain" ~ "a",
        establishmentMeans == "native" ~ "b",
        establishmentMeans == "vagrant" ~ "c",
        TRUE ~ NA_character_  # Assign NA for any other values
      ))
    
    # ## reframe (i.e. mutate) based on the created flag
    alldat <- alldat %>%
      reframe(.by = c(location, locationID, taxon, scientificName, taxonID, dist_group), # group by taxon, location, and dist_group 
              across(everything(),  # concatenate values across all other columns after the defined grouping (taxon, location, and dist_group)
                     ~ paste(unique(.x), collapse = "; ")  # concatenate unique values  
      )) %>%
      select(-dist_group)  # Remove custom column to flag records
  }

  #### merge rows ############################
  
  # clean occurrenceStatus for taxa reported as both absent and present in a location (for conflict between DBs)
  alldat$occurrenceStatus <- ifelse(alldat$occurrenceStatus == "absent; present", "present", alldat$occurrenceStatus)
  
  ## merge remaining content of rows defined by the "by" statement
  dt <- as.data.table(alldat)
  dt$eventDate <- as.character(dt$eventDate)
  alldat <- dt[,lapply(.SD,function(s) paste(s,collapse="; ")),by=list(taxon,location,locationID,taxonID,scientificName,establishmentMeans)]

  ## standardise column entries
  alldat <- as.data.frame(alldat,stringsAsFactors=F)  
  alldat$eventDate <- as.character(alldat$eventDate)

  # clean column entries
  ind_col <- which(!colnames(alldat)%in%c("taxon","taxonID","location","locationID","scientificName","establishmentMeans"))

  for (i in ind_col){
    alldat[,i] <- unlist(lapply(strsplit(alldat[,i],"; "),function(s) paste(sort(unique(s)),collapse="; ")))
    alldat[,i] <- gsub("; NA","",alldat[,i]) # clean new column
    alldat[,i] <- gsub("NA; ","",alldat[,i]) # clean new column
    alldat[,i] <- gsub("^; ","",alldat[,i]) # clean new column
    alldat[,i] <- gsub("; $","",alldat[,i]) # clean new column
    alldat[,i] <- gsub(";$","",alldat[,i]) # clean new column
    alldat[,i] <- gsub(" ; "," ",alldat[,i]) # clean new column
    alldat[,i] <- gsub("^NA$","",alldat[,i]) # clean new column
  }
  
  ## select minimum first record
  ind_dupl_fr <- grep("; ",alldat$eventDate)
  single_fr <- lapply(strsplit(alldat$eventDate[ind_dupl_fr],"; "),function(s) min(as.integer(s)))
  alldat$eventDate[ind_dupl_fr] <- unlist(single_fr)
  
  ## remove duplicated entries
  ind_dupl <- duplicated(alldat) # remove identical lines
  alldat <- alldat[!ind_dupl,]
  
  ## specify columns to check for duplicates 
  ind_dupl <- duplicated(alldat[,c("taxon","location","establishmentMeans")]) # remove duplicated species-region combinations
  all_dat_dupl <- unique(alldat[ind_dupl,c("taxon","location","establishmentMeans")]) # remove duplicated species-region combinations

  ind_rm <- col_dupl <- vector()
  for (j in 1:nrow(all_dat_dupl)){ # loop over duplicated entries
    
    ind_each <- which(alldat$taxon == all_dat_dupl$taxon[j] & alldat$location == all_dat_dupl$location[j] & alldat$establishmentMeans == all_dat_dupl$establishmentMeans[j])

    for (k in 1:dim(alldat)[2]){ # loop over columns
      if (colnames(alldat)[k]%in%c("location","locationID","taxon","scientificName","taxonID","establishmentMeans")) next # ignore these columns
      if (all(is.na(alldat[ind_each,k]))) next # skip if all NA
      if (all(duplicated(alldat[ind_each,k])[-1])){ # skip if all equal (non-NA)
        next 
      } else {
        ind_NA <- is.na(alldat[ind_each,k]) # avoid NAs
        if (length(unique(alldat[ind_each,k][!ind_NA]))==1){ # if only a single value is non-NA, add non-NA information to first row
          alldat[ind_each,k][1] <- alldat[ind_each,k][!ind_NA][1]
        } else { # if deviating entries exit, merge content into first row
          if (colnames(alldat)[k]=="eventDate"){ # treat first records differently (no merging)
            alldat[ind_each,k][1] <- min(alldat[ind_each,k],na.rm = T) # select earliest first record
          } else {
            entries <- unique(alldat[ind_each,k][!ind_NA])
            entries <- entries[entries!=""]
            alldat[ind_each,k][1] <- paste(entries,collapse="; ") # concatenate unequal row entries
            col_dupl <- c(col_dupl,colnames(alldat)[k]) # store column with deviating information for report
          }
        }
        ind_rm <- c(ind_rm,ind_each[-1]) # collect rows to remove (all except the first one)
      }
    }
  }
  if (length(ind_rm)>0) alldat <- alldat[-unique(ind_rm),] # remove all duplicates
  if (length(col_dupl)>0) cat(paste0("\n    Note: Multiple entries for the same record. Check column '",unique(col_dupl),"' in final data set for entries separated by ';'. \n"))
  
  # fill habitat column
  alldat <- alldat %>%
    group_by(taxonID) %>%
    mutate(
      habitat = habitat[which.max(nchar(habitat, keepNA = FALSE))])
  
  ## output #############################################

  ## sort by taxonomic tree
  fulltaxalist <- read.table(file.path("Output",paste0(outputfilename,"_",version,"_","FullTaxaList.csv")),stringsAsFactors = F,header=T)
  fulltaxalist <- unique(fulltaxalist[,c("taxonID","kingdom","phylum","class","order","family")])
  alldat <- merge(alldat,fulltaxalist,by="taxonID",all.x=T)
  
  oo <- order(alldat$location,alldat$kingdom,alldat$phylum,alldat$class,alldat$order,alldat$family,alldat$scientificName)

  alldat <- alldat[oo,]
  
  ## identify columns for output
  all_addit_cols <- paste(FileInfo$Column_additional,collapse="; ")
  all_addit_cols <- unlist(strsplit(all_addit_cols,"; "))
  all_addit_cols <- all_addit_cols[all_addit_cols!="NA"]
  
  # identify columns for output depending on desired spatial aggregation
  key_columns <- c("location", "locationID", "taxon", "scientificName", "taxonID", "eventDate", 
                     "habitat", "occurrenceStatus", "establishmentMeans", "degreeOfEstablishment", 
                     "pathway", "datasetName", "bibliographicCitation")

  
  columns_out <- colnames(alldat)[colnames(alldat) %in% c(key_columns, all_addit_cols)]
  
  ## Reorder columns: prioritize key columns, then additional columns
  ind <- match(key_columns, columns_out)
  columns_out <- unique(c(columns_out[ind[which(!is.na(ind))]], columns_out[which(is.na(ind))], all_addit_cols))
  
  ## final cleaning
  columns_out <- columns_out[!is.na(columns_out)]
  
  alldat_out <- alldat[,columns_out]
  alldat_out[is.na(alldat_out)] <- ""
  alldat_out[alldat_out=="NA"] <- ""
  
  ## set final names for the columns
  colnames(alldat_out)[colnames(alldat_out) == "FirstRecord_orig"] <- "eventDate_orig"
  
  ## remove scientificName from main output file (only provide in the taxonomic table)
  alldat_out <- subset(alldat_out, select= -scientificName)
  
  ## write main output of workflow
  setDT(alldat_out) # Convert output to a data.table if it isn’t already

  ## Remove additional characters "\" from the reference column
  alldat_out$bibliographicCitation <- gsub('["\']', '', alldat_out$bibliographicCitation)

  ## store final output file
  write.table(alldat_out,file.path("Output",paste(outputfilename,"_",version,".csv",sep="")),row.names=F)

  ## ending line
  cat(paste("\n Successfully established version",version,"of",outputfilename,"file. \n"))
  
  ## delete intermediate results if selected ########
  if (!output) {
    interm_res <- list.files(file.path("Output","Intermediate"))
    interm_res <- interm_res[grep("Standard",interm_res)]
    file.remove(paste0("Output",interm_res))
  }
}
