
####### SInAS workflow: Integration and standardisation of alien species data ###########
##
## Step 2c: Standardisation of species names using the GBIF backbone taxonomy
##
## Species names are standardised according to the GBIF backbone taxonomy. The protocol 
## to access GBIF and treat results is implemented in CheckGBIFTax.r.
## Script requires internet connection.
##
## Hanno Seebens, Giessen, 25.09.2025
#########################################################################################


StandardiseTaxonNames <- function(FileInfo=NULL){

  ## load data #########################################################
  
  ## identify input datasets based on file name "StandardColumns_....csv"
  allfiles <- list.files(file.path("Output","Intermediate"))
  inputfiles_all <- allfiles[grep("Step3_StandardLocationNames_",allfiles)]
  
  inputfiles <- vector()
  for (i in 1:nrow(FileInfo)){
    inputfiles <- c(inputfiles,paste("Step3_StandardLocationNames_",FileInfo[i,"Dataset_brief_name"],".csv",sep=""))
  }
  inputfiles <- inputfiles[!is.na(inputfiles)]
  

  ## loop over all data sets ################################################
  
  fullspeclist <- data.frame()
  for (i in 1:length(inputfiles)){ # loop over inputfiles 
    
    dat <- read.table(file.path("Output","Intermediate",paste0(inputfiles[i])),header=T,stringsAsFactors = F)
    
    dat <- dat[,!colnames(dat)%in%c("location_orig")]
    
    # remove white space #######################################
    dat$originalNameUsage <- gsub("  "," ",dat$originalNameUsage)
    dat$originalNameUsage <- gsub("^\\s+|\\s+$", "",dat$originalNameUsage) # trim leading and trailing whitespace
    dat$originalNameUsage <- gsub("[$\xc2\xa0]", " ",dat$originalNameUsage) # replace weird white space with recognised white space
    dat$originalNameUsage <- gsub("  "," ",dat$originalNameUsage)
    dat$originalNameUsage <- gsub("\n"," ",dat$originalNameUsage)
    
    dat <- dat[!is.na(dat$originalNameUsage),]
    dat <- dat[dat$originalNameUsage!="",]
    
    #### standardise taxon names against GBIF backbone taxonomy ###

    cat(paste0("\n    Working on ",FileInfo[i,"Dataset_brief_name"],"... \n"))
    dat2 <- CheckGBIFTaxParallel(taxon_names = dat$originalNameUsage)

    ## output #####################################################
    
    DB <- unique(dat2[[1]])
    mismatches <-  unique(dat2[[2]][order(dat2[[2]]$taxon)])
    
    ## merge standardized names to OG df
    DB <- merge(dat, DB, by="originalNameUsage", all=T)
    
    ## export full species list with original species names and names assigned by GBIF for checking
    fullspeclist <- rbind(fullspeclist,unique(DB[,c("originalNameUsage","taxon","scientificName","GBIFstatus","GBIFstatus_Synonym","GBIFtaxonRank","GBIFusageKey","workflowNote","species","genus","family","order","class","phylum","kingdom")]))
    
    DB <- unique(DB) # remove duplicates
    DB$GBIFstatus[is.na(DB$GBIFstatus)] <- "NoMatch"
    DB <- DB[,!colnames(DB)%in%c("GBIFstatus","GBIFconfidence","GBIFmatchType","GBIFtaxonRank","GBIFusageKey","workflowNote","GBIFstatus_Synonym","species","genus","family","class","order","phylum","kingdom")]
    
    ## output main database with standardised taxon names
    write.table(DB,file.path("Output","Intermediate",paste0("Step4_StandardTaxonNames_",FileInfo[i,"Dataset_brief_name"],".csv")),row.names=F)
    
    oo <- order(mismatches$taxon)
    mismatches <- unique(mismatches[oo,])
    
    ## output unmatched names for cross-checking
    write.table(mismatches,file.path("Output","Check",paste0("Missing_Taxa_",FileInfo[i,"Dataset_brief_name"],".csv")),row.names=F,col.names=F)
  }
  
  ## Prepare a full list of taxon names across all data sets incl. IDs #########
  
  oo <- order(fullspeclist$kingdom,fullspeclist$phylum,fullspeclist$class,fullspeclist$order,fullspeclist$taxon)
  fullspeclist <- unique(fullspeclist[oo,])
  
  ## assign taxon ID unique to individual taxa #############
  ## identify unique taxa (obtained from GBIF)
  fullspeclist$sequence <- 1:nrow(fullspeclist)
  uni_taxa <- unique(fullspeclist$scientificName)
  uni_taxa <- data.frame(scientificName=uni_taxa[!is.na(uni_taxa)],stringsAsFactors=F)
  uni_taxa$taxonID <- 1:nrow(uni_taxa)

  ## define major taxonomic groups ###################
  fullspeclist$taxaGroup <- NA
  fullspeclist$taxaGroup[fullspeclist$originalNameUsage%in%subset(fullspeclist, phylum%in%c("Foraminifera","Cercozoa","Ciliophora","Ochrophyta","Oomycota","Myzozoa","Peronosporea", "Bigyra"))$originalNameUsage] <- "SAR"
  fullspeclist$taxaGroup[fullspeclist$originalNameUsage%in%subset(fullspeclist, genus%in%c("Plasmodium"))$originalNameUsage] <- "SAR"
  fullspeclist$taxaGroup[fullspeclist$originalNameUsage%in%subset(fullspeclist,class%in%"Mammalia")$originalNameUsage] <- "Mammals"
  fullspeclist$taxaGroup[fullspeclist$originalNameUsage%in%subset(fullspeclist,class%in%"Aves")$originalNameUsage] <- "Birds"
  fullspeclist$taxaGroup[fullspeclist$originalNameUsage%in%subset(fullspeclist,class%in%c("Cephalaspidomorphi","Actinopterygii","Elasmobranchii","Sarcopterygii", "Petromyzonti", "Chondrichthyes", "Myxini"))$originalNameUsage] <- "Fishes"
  fullspeclist$taxaGroup[fullspeclist$originalNameUsage%in%subset(fullspeclist,order%in%c("Polypteriformes", "Acipenseriformes", "Lepisosteiformes", "Amiiformes", "Osteoglossiformes", "Hiodontiformes", 
                                                                                       "Elopiformes", "Albuliformes", "Notacanthiformes", "Anguilliformes", "Saccopharyngiformes", "Clupeiformes", "Ceratodontiformes",
                                                                                       "Gonorynchiformes", "Cypriniformes", "Characiformes", "Gymnotiformes", "Siluriformes", "Salmoniformes", "Esociformes", 
                                                                                       "Osmeriformes", "Ateleopodiformes", "Stomiiformes", "Aulopiformes", "Myctophiformes", "Lampriformes", "Polymixiiformes", 
                                                                                       "Percopsiformes", "Batrachoidiformes", "Lophiiformes", "Gadiformes", "Ophidiiformes", "Mugiliformes", "Atheriniformes", 
                                                                                       "Beloniformes", "Cetomimiformes", "Cyprinodontiformes", "Stephanoberyciformes", "Beryciformes", "Zeiformes", 
                                                                                       "Gobiesociformes", "Gasterosteiformes", "Syngnathiformes", "Synbranchiformes", "Tetraodontiformes", "Pleuronectiformes", 
                                                                                       "Scorpaeniformes", "Perciformes"))$originalNameUsage] <- "Fishes"
  fullspeclist$taxaGroup[fullspeclist$originalNameUsage%in%subset(fullspeclist,class%in%c("Reptilia", "Testudines","Squamata", "Crocodylia"))$originalNameUsage] <- "Reptiles"
  fullspeclist$taxaGroup[fullspeclist$originalNameUsage%in%subset(fullspeclist,class%in%"Amphibia")$originalNameUsage] <- "Amphibians"
  fullspeclist$taxaGroup[fullspeclist$originalNameUsage%in%subset(fullspeclist, class%in%c("Insecta"))$originalNameUsage] <- "Insects"
  fullspeclist$taxaGroup[fullspeclist$originalNameUsage%in%subset(fullspeclist, class%in%c("Arachnida", "Pycnogonida"))$originalNameUsage] <- "Arachnids"
  fullspeclist$taxaGroup[fullspeclist$originalNameUsage%in%subset(fullspeclist, class%in%c("Collembola", "Chilopoda", "Diplopoda", "Diplura", "Merostomata", "Pauropoda", "Protura", "Symphyla"))$originalNameUsage] <- "Other arthropods"
  fullspeclist$taxaGroup[fullspeclist$originalNameUsage%in%subset(fullspeclist,class%in%c("Branchiopoda","Hexanauplia","Maxillopoda","Ostracoda","Malacostraca", "Copepoda"))$originalNameUsage] <- "Crustaceans"
  fullspeclist$taxaGroup[fullspeclist$originalNameUsage%in%subset(fullspeclist,family%in%c("Elminiidae"))$originalNameUsage] <- "Crustaceans"
  fullspeclist$taxaGroup[fullspeclist$originalNameUsage%in%subset(fullspeclist,phylum%in%"Mollusca")$originalNameUsage] <- "Molluscs"
  fullspeclist$taxaGroup[fullspeclist$originalNameUsage %in% subset(fullspeclist, phylum %in% "Tracheophyta")$originalNameUsage] <- "Vascular plants"
  fullspeclist$taxaGroup[fullspeclist$originalNameUsage%in%subset(fullspeclist,phylum%in%c("Bryophyta","Anthocerotophyta", "Marchantiophyta"))$originalNameUsage] <- "Bryophytes"
  fullspeclist$taxaGroup[fullspeclist$originalNameUsage%in%subset(fullspeclist,phylum%in%c("Rhodophyta","Chlorophyta","Charophyta","Cryptophyta","Haptophyta"))$originalNameUsage] <- "Algae"
  fullspeclist$taxaGroup[fullspeclist$originalNameUsage%in%subset(fullspeclist,phylum%in%c("Ascomycota", "Dothideomycetes", "Sordariomycetes", "Chytridiomycota","Basidiomycota","Microsporidia","Zygomycota", "Entomophthoromycota"))$originalNameUsage] <- "Fungi"
  fullspeclist$taxaGroup[fullspeclist$originalNameUsage%in%subset(fullspeclist,phylum%in%c("Actinobacteria","Chlamydiae","Cyanobacteria","Firmicutes","Proteobacteria"))$originalNameUsage |
                           fullspeclist$class == "Ichthyosporea"] <- "Bacteria and protozoans"
  fullspeclist$taxaGroup[fullspeclist$originalNameUsage%in%subset(fullspeclist,kingdom%in%c("Bacteria", "Protozoa","Euglenozoa"))$originalNameUsage] <- "Bacteria and protozoans"
  fullspeclist$taxaGroup[fullspeclist$originalNameUsage%in%subset(fullspeclist,kingdom%in%c("Viruses"))$originalNameUsage] <- "Viruses"
  fullspeclist$taxaGroup[fullspeclist$originalNameUsage%in%subset(fullspeclist,phylum%in%c("Annelida", "Nematoda", "Platyhelminthes", "Sipuncula", "Nemertea", "Onychophora", "Acanthocephala"))$originalNameUsage] <- "Annelids, nematodes, platyhelminthes, and other worms"
  fullspeclist$taxaGroup[fullspeclist$originalNameUsage%in%subset(fullspeclist,phylum%in%c("Bryozoa", "Entoprocta", "Chaetognatha", "Cnidaria", "Ctenophora", "Echinodermata", "Phoronida", "Porifera", "Rotifera", "Xenacoelomorpha","Brachiopoda"))$originalNameUsage] <- "Other aquatic animals"
  fullspeclist$taxaGroup[fullspeclist$originalNameUsage%in%subset(fullspeclist,class%in%"Ascidiacea")$originalNameUsage] <- "Other aquatic animals"

  
  
  ## merge taxonID with full taxa list
  fullspeclist_2 <- merge(fullspeclist,uni_taxa,by="scientificName",all=T)
  fullspeclist_2$taxonID[is.na(fullspeclist_2$taxonID)] <- 1:length(which(is.na(fullspeclist_2$taxonID))) +  max(fullspeclist_2$taxonID,na.rm=T)
  
  fullspeclist_2 <- fullspeclist_2[order(fullspeclist_2$sequence),]
  fullspeclist_2 <- fullspeclist_2[,-which(colnames(fullspeclist_2)=="sequence")]
  #fullspeclist_2 <- fullspeclist_2[, !c("sequence"), with = FALSE]
  
  
  write.table(fullspeclist_2,file.path("Output",paste0(outputfilename,"_",version,"_","FullTaxaList.csv")),row.names=F)
              
  
  ## add taxon ID to data sets ##########
  
  taxon_id <- unique(fullspeclist_2[,c("taxonID","originalNameUsage")])
  for (i in 1:length(inputfiles)){ # loop over inputfiles 
    dat <- read.table(file.path("Output","Intermediate",paste0("Step4_StandardTaxonNames_",FileInfo[i,"Dataset_brief_name"],".csv")),header=T,stringsAsFactors = F)
    dat <- merge(dat,taxon_id,by="originalNameUsage",all.x=T)
    
    dat <- dat[!is.na(dat$taxon) & dat$taxon != "", ] # remove taxa with no match in GBIF from the intermediate output file
    
    ## save intermediate output all data sets 
    write.table(dat,file.path("Output","Intermediate",paste0("Step4_StandardTaxonNames_",FileInfo[i,"Dataset_brief_name"],".csv")),row.names=F)
  }  
}
