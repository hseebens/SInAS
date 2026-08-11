####### SInAS workflow: Integration and standardisation of alien species data ###########
##
## Step 2c: Standardisation of taxon names
## check and replace species names using 'rgibf' GBIF taxonomy
##
## Hanno Seebens, Giessen, 19.02.2025
#########################################################################################

CheckGBIFTaxParallel <- function(taxon_names=NULL,
                                     column_name_taxa=NULL){
  
  ## check input variable
  if (is.null(taxon_names)){
    stop("No taxon names provided.")
  } else if (is.character(taxon_names)){ # check if input file is a vector
    dat <- as.data.frame(taxon_names)
    colnames(dat) <- "taxon"
  } else if (is.data.frame(taxon_names)){ # check if input file is a data.frame
    dat <- taxon_names
  } else {
    stop("Cannot coerce data into data.frame. Please provide a data.frame or vector as input.")
  }
  
  if (!is.null(column_name_taxa)){ # check if column name of taxa provided
    colnames(dat)[colnames(dat)==column_name_taxa] <- "taxon" # rename to standard column name
  }
  if (all(colnames(dat)!="taxon")){ # check if column "taxon_orig" can be found
    stop("No column with taxon names found. Please specify in column_name_taxa.")
  }
  
  if (any(colnames(dat)=="kingdom_user")){
    taxlist_lifeform <- unique(dat[,c("taxon","kingdom_user")])
    taxlist <- taxlist_lifeform$taxon
  } else if (any(colnames(dat)=="Author")){
    taxlist <- unique(paste(dat$taxon,dat$Author))
  } else {
    taxlist <- unique(dat$taxon)
  }
  n_taxa <- length(taxlist)
  cat(" - The number of taxa to be processed is", n_taxa, "\n")

  ## get taxonomic information from GBIF #######################################
  
  cat(" - Retrieve information from GBIF...\n")
  # Rprof()
  # check: Acentropus niveus; Acaena spp (multiple accepted names); Acacia fimbriata GDon; Acacia nilotica (L) Del subsp nilotica
  
  ##############################################################################
  
  ## assign taxonomic information to records ###################################
  ## Helper functions to avoid redundancy in the code ##########################
  ## function to fill cells related to taxonomy
  fill_taxonomy <- function(src, dst, i = 1) {
    fields <- c("species", "genus", "family", "class", "order", "phylum", "kingdom")
    for (f in fields) {
      if (f %in% colnames(src)) {
        dst[[f]] <- src[[f]][i]
      }
    }
    dst
  }
  
  ## function to fill cells related to GBIF matches
  fill_gbif_core <- function(src, dst, i = 1) {
    dst$GBIFstatus     <- src$status[i]
    dst$GBIFmatchType  <- src$matchType[i]
    dst$GBIFtaxonRank  <- src$rank[i]
    dst$GBIFusageKey   <- src$usageKey[i]
    dst$GBIFconfidence   <- src$confidence[i]
    dst
  }
 
  cores=detectCores()
  cl <- makeCluster(cores[1]-1) # -1 to avoid overloading your computer
  registerDoParallel(cl)

  all_out <- foreach(i=1:n_taxa, .packages=c("rgbif")) %dopar% {
  # for (i in 1:length(db_all)){ #
  
    gbif_entry <- name_backbone_verbose(taxlist[i], strict=F) # check for names and synonyms
    # gbif_entry <- name_backbone_verbose("uronema marinum", strict=F) # check for names and synonyms
    
    originalNameUsage <- attr(gbif_entry, "args")$name  
    
    out_entry <- data.frame(originalNameUsage, taxon=NA, scientificName=NA, 
                            species=NA, genus=NA, family=NA, class=NA, order=NA, phylum=NA, kingdom=NA, 
                            GBIFusageKey=NA, GBIFtaxonRank=NA, GBIFstatus= "NoMatch", 
                            GBIFmatchType=NA, workflowNote=NA, GBIFstatus_Synonym=NA, GBIFconfidence=NA)

    # select species name and download taxonomy
    db <- gbif_entry$data
    alternatives <- gbif_entry$alternatives
    
    if ("status"%in%colnames(db)){
      if ((db$status=="ACCEPTED" | db$status=="DOUBTFUL") & db$matchType=="EXACT"  & "canonicalName"%in%colnames(db)) { 
        
        ##########################################################################
        ### EXACT MATCHES: select accepted or doubtful names and exact matches ###
        ##########################################################################
        
        ## define criterion for selection
        # criterion <- db$matchType=="EXACT"
        
        out_entry$taxon      <- db$canonicalName[1]
        out_entry$scientificName <- db$scientificName[1]
        
        ## fill cells related to GBIF matches
        out_entry <- fill_gbif_core(db, out_entry)
        
        ## fill cells related to taxonomy
        out_entry <- fill_taxonomy(db, out_entry)
        
        out_entry$workflowNote <- "case 1"
        
        return(out_entry)
        
      } else if ((db$status=="SYNONYM" | db$status=="HETEROTYPIC_SYNONYM") & "species"%in%colnames(db)) { # select synonyms
        
        ################################################################################      
        ## SYNONYMS ####################################################################
        ################################################################################      
        
        if ("acceptedUsageKey"%in%colnames(db)) {  ## try to get author name of synonym (not provided in 'db')(works only for species)
          #if (db$rank=="SPECIES")
          
          synonym_info <- name_usage(key=db$acceptedUsageKey)[[2]]
          
          db2_all <- name_backbone_verbose(synonym_info$scientificName) #$data

          if (db2_all$data$matchType!="NONE"){ # check if data set exists
            db2 <- db2_all$data
          } else {
            db2 <- db2_all$alternatives
          }
            
          ## fill cells related to GBIF matches
          out_entry <- fill_gbif_core(db, out_entry) # original match
          
          ## if entries in new match exists, fill cells
          if (any((db2$status=="ACCEPTED" | db2$status=="DOUBTFUL") & db2$matchType=="EXACT")){
            
            criterion <- (db2$status=="ACCEPTED" | db2$status=="DOUBTFUL") & db2$matchType=="EXACT" 
            
            if (sum(criterion)>1){ # select entry with largest confidence
              ind <- max(db2$confidence[criterion])
              criterion <- (db2$status=="ACCEPTED" | db2$status=="DOUBTFUL") & db2$matchType=="EXACT" & db2$rank=="SPECIES" & db2$confidence==ind
            }
            
            ## fill cells related to taxonomy
            out_entry <- fill_taxonomy(db2[criterion,], out_entry) # new match
            out_entry$GBIFstatus_Synonym <- db2$status[criterion] # include status of synonym
            
            out_entry$taxon    <- db2[criterion, ]$canonicalName
            out_entry$scientificName    <- db2[criterion, ]$scientificName

            out_entry$workflowNote        <- "case 2a"
            
            return(out_entry)
            
          } else { # if not, use species level information from original match
            
            out_entry$taxon    <- db[(db$status=="SYNONYM"), ]$species # scientificName is from the original synonym
            
            out_entry$workflowNote        <- "case 2b"
            
            return(out_entry)
            
          }

          out_entry$workflowNote        <- "case 2c"
          
          return(out_entry)
          
        }
      } else if (any(db$matchType=="FUZZY" & db$confidence>50 & "canonicalName"%in%colnames(db))) { 
        
        ##########################################################################      
        ## FUZZY MATCHES #########################################################
        ##########################################################################    

        out_entry$taxon      <- db$canonicalName[1]
        out_entry$scientificName <- db$scientificName[1]
        
        ## fill cells related to GBIF matches
        out_entry <- fill_gbif_core(db, out_entry)
        ## fill cells related to taxonomy
        out_entry <- fill_taxonomy(db, out_entry)
        
        out_entry$workflowNote <- "case 3"
        
        return(out_entry)
        
      } else if (db$matchType=="HIGHERRANK") { # if nothing found at that taxonomic level 
        
        ########################################################################      
        ## HIGHERRANK: check if lower rank information exists ##################
        ########################################################################    

        # if (any(alternatives$matchType=="EXACT")){ # check for exact matches in alternatives
        #   criterion <- min(which(alternatives$matchType=="EXACT")) # take the first records, which seemed to be the best match (based on experience)
        #   if (alternatives$status[criterion]=="ACCEPTED"){
        if (nrow(alternatives)>0){
          # criterion <- alternatives$matchType=="EXACT" & alternatives$status=="ACCEPTED" #& which(alternatives$confidence==max(alternatives$confidence, na.rm=T))
          criterion <- 1 # first row seems to match best if confidence is sufficiently high
          if (alternatives$confidence[criterion]>60){
            if (alternatives$status[criterion]=="ACCEPTED"){
              out_entry$taxon      <- alternatives[criterion,]$canonicalName[1]
              out_entry$scientificName <- alternatives[criterion,]$scientificName[1]
              
              ## fill cells related to GBIF matches
              out_entry <- fill_gbif_core(alternatives[criterion,] , out_entry)
              ## fill cells related to taxonomy
              out_entry <- fill_taxonomy(alternatives[criterion,], out_entry)
              
              # check if taxon rank is correct (sometimes GBIF reports SPECIES rather than GENUS in alternatives)
              db2_all <- name_backbone_verbose(out_entry$scientificName) # get new match with new taxon name
              if (db2_all$data$matchType!="NONE"){
                if (db2_all$data$rank!=out_entry$GBIFtaxonRank){
                  out_entry$GBIFtaxonRank <- db2_all$data$rank # overwrite rank if deviating
                }
              }
              if (db2_all$data$matchType=="NONE"){
                if (db2_all$alternatives$rank[1]!=out_entry$GBIFtaxonRank){
                  out_entry$GBIFtaxonRank <- db2_all$alternatives$rank[1] # overwrite rank if deviating
                }
              }
              out_entry$workflowNote <- "case 4a"
              return(out_entry)
            }
            
            ## check if information of synonyms is provided in alternatives
            if (alternatives$status[criterion]=="SYNONYM" | alternatives$status[criterion]=="HETEROTYPIC_SYNONYM"){
              # criterion <- alternatives$status=="SYNONYM" & alternatives$matchType=="EXACT" #& which(alternatives$confidence==max(alternatives$confidence, na.rm=T))
              criterion <- 1

              # check if taxon rank is correct (sometimes GBIF reports SPECIES rather than GENUS in alternatives)
              db2_all <- name_backbone_verbose(alternatives$species[criterion]) # get new match with new taxon name
              
              if (db2_all$data$matchType!="NONE"){
                
                if ("species"%in%colnames(db2_all$data)){
                  
                  if (db2_all$data$matchType=="EXACT" & db2_all$data$status=="ACCEPTED"){
                    
                    out_entry$taxon      <- db2_all$data$canonicalName # take the first entry which seemingly is the best match (independent of confidence level)
                    out_entry$scientificName <- db2_all$data$scientificName
                    
                    ## fill cells related to GBIF matches
                    out_entry <- fill_gbif_core(db2_all$data, out_entry)
                    ## fill cells related to taxonomy
                    out_entry <- fill_taxonomy(db2_all$data, out_entry)
                    
                    out_entry$GBIFstatus <- "SYNONYM"
                    out_entry$GBIFstatus_Synonym <- db2_all$data$status
                    out_entry$workflowNote <- "case 4b"
                    return(out_entry)
                  }
                }
              }
            }
          }
        } else {

          ########################################################################
          ## ACCEPT ALSO HIGHERRANK if canonicalName exists ######################
          ########################################################################
          # eg. subspec to species level
          
          ## define criterion for selection
          # criterion <- (db$status=="ACCEPTED" | db$status=="DOUBTFUL") & db$matchType=="HIGHERRANK"
          
          # out_entry$taxon      <- db$canonicalName[1]
          # out_entry$scientificName <- db$scientificName[1]
          
          # ## fill cells related to GBIF matches
          # out_entry <- fill_gbif_core(db, out_entry)
          ## fill cells related to taxonomy
          out_entry <- fill_taxonomy(db, out_entry)
          out_entry$workflowNote <- "case 4c"
          return(out_entry)
        }
      }      
    }
    if (nrow(alternatives)>0){
      
      criterion <- 1
      
      if (alternatives$confidence[criterion]>60){
        
        if (alternatives$status[criterion]=="ACCEPTED"){
          
          out_entry$taxon <- alternatives[criterion,]$canonicalName
          out_entry$scientificName <- alternatives[criterion,]$scientificName
          
          ## fill cells related to GBIF matches
          out_entry <- fill_gbif_core(alternatives[criterion,] , out_entry)
          ## fill cells related to taxonomy
          out_entry <- fill_taxonomy(alternatives[criterion,], out_entry)
          
          out_entry$workflowNote <- "case 5a"
          
          return(out_entry)
          
        }
        if (alternatives$status[criterion]=="SYNONYM"){
          
          # check if taxon rank is correct (sometimes GBIF reports SPECIES rather than GENUS in alternatives)
          db2_all <- name_backbone_verbose(alternatives$species[criterion]) # get new match with new taxon name
          
          if (db2_all$data$matchType!="NONE"){
            
            if ("species"%in%colnames(db2_all$data)){
              
              if (db2_all$data$matchType=="EXACT" & db2_all$data$status=="ACCEPTED"){
                
                out_entry$taxon      <- db2_all$data$canonicalName # take the first entry which seemingly is the best match (independent of confidence level)
                out_entry$scientificName <- db2_all$data$scientificName
                
                ## fill cells related to GBIF matches
                out_entry <- fill_gbif_core(db2_all$data, out_entry)
                ## fill cells related to taxonomy
                out_entry <- fill_taxonomy(db2_all$data, out_entry)
                out_entry$GBIFstatus <- "SYNONYM"
                out_entry$GBIFstatus_Synonym <- db2_all$data$status
                out_entry$workflowNote <- "case 5b"
                return(out_entry)
              }
            }
          }
        }
      }
    } 
    if ("status"%in%colnames(db)){
      
      out_entry$taxon = originalNameUsage
      out_entry$GBIFstatus = if (!is.null(db$status)) as.character(db$status[1]) else NA_character_
      out_entry$GBIFmatchType = if (!is.null(db$matchType)) as.character(db$matchType[1]) else NA_character_
      
      if (db$status=="ACCEPTED"){
        ## fill cells related to GBIF matches
        out_entry <- fill_gbif_core(db , out_entry)
      }
      
      ## fill cells related to taxonomy
      out_entry <- fill_taxonomy(db, out_entry)
      
      out_entry$workflowNote <- "case 6"
      
      return(out_entry)
      
    } else {
      
      out_entry$workflowNote <- "case 7"
      
      return(out_entry)
      
    }
  }
  ## stop cluster  
  stopCluster(cl)


  all_out_dt <- rbindlist(all_out, fill=TRUE)
  missing <- subset(all_out_dt, is.na(taxon))
  
  # Rprof(NULL)    ## Turn off the profiler
  # summaryRprof()

  n_taxa_matched <- sum(!is.na(all_out_dt$taxon))
  n_synonyms <- sum(all_out_dt$GBIFstatus=="SYNONYM", na.rm=T)

  out <- list()
  out[[1]] <- all_out_dt
  out[[2]] <- missing
  
  # cat("Step 2a completed: species names have been standardized across the GBIF backbone taxonomy\n")
  cat(" -", n_synonyms, "synonyms identified\n")
  cat(" -", n_taxa_matched, "(",round(n_taxa_matched/n_taxa,2)*100, "% )", "taxon names retrieved from GBIF\n")
  cat(" - For", n_taxa-n_taxa_matched, "taxa no match was found and the taxon names provided by the user were taken")

  return(out)
}
