
####### SInAS workflow: Integration and standardisation of alien species data ###########
##
## Step 1: Prepare databases of alien taxon distribution and first records
## as input datasets to create a merged database
## 
## Hanno Seebens, Giessen, 25.09.2025
################################################################################


PrepareDatasets <- function(FileInfo=NULL){
  
  ## create output folder #####
  if (!file.exists("Output")){
    dir.create("Output")
    dir.create(file.path("Output","Intermediate"))
    dir.create(file.path("Output","Check"))
  }
  
  ######## Load data sets and standardise columns names for each ######
  
  for (i in 1:nrow(FileInfo)){#
    
    ## load data set
    data_name <- FileInfo[i,"File_name_to_load"]
    if (!file.exists(file.path("Inputfiles",data_name))) stop(paste0("File ’",data_name,"’ does not exist!"))
    
    if (grepl("\\.xlsx$",data_name)){ # xlsx file
      dat <- try(read.xlsx(file.path("Inputfiles",data_name),sheet=1),silent=T)
    } 
    
    if (grepl("\\.csv$",data_name)){ # csv file
      # dat <- try(read.csv(file.path("Inputfiles",data_name)),silent=T)
      dat <- try(as.data.frame(fread(file.path("Inputfiles",data_name)),stringsAsFactors=F),silent=T)
      if (dim(dat)[2]<3)  dat <- try(as.data.frame(fread(file.path("Inputfiles",data_name)),stringsAsFactors=F),silent=T)
    }
    
    if (class(dat)=="try-error") stop(paste0("File ’",data_name,"’ should be csv or xlsx format"))

    ## store original column names
    col_names_import <- colnames(dat)

    ## check and rename column names if column is provided
    all_column_names <- vector()
    
    ## Column recordID
    if (!is.na(FileInfo[i,"Column_recordID"]) & FileInfo[i,"Column_recordID"]!=""){
      col_recordID <- FileInfo[i,"Column_recordID"]
      colnames(dat)[col_names_import==col_recordID] <- paste("recordID",FileInfo[i,"Dataset_brief_name"],sep="_")
      all_column_names <- c(all_column_names,paste("recordID",FileInfo[i,"Dataset_brief_name"],sep="_"))
    }
    
    ## Column taxon
    if (!is.na(FileInfo[i,"Column_taxon"]) & FileInfo[i,"Column_taxon"]!=""){
      col_spec_names <- FileInfo[i,"Column_taxon"]
      all_column_names <- col_spec_names
      if (is.na(col_spec_names)) stop(paste("Column with taxon names not found in",FileInfo[i,"Dataset_brief_name"],"file!"))
      if (!is.na(FileInfo[i,"Column_author"]) & FileInfo[i,"Column_author"]!=""){
        col_author <- FileInfo[i,"Column_author"]
        # all_column_names <- c(all_column_names,"Author")
        dat[,col_spec_names] <- paste(dat[,col_spec_names],dat[,col_author]) # add author to taxon name
        dat[,col_spec_names] <- gsub(" NA","",dat[,col_spec_names]) # remove missing author names
      }
    }
    
    ## Column scientificName
    if (!is.na(FileInfo[i,"Column_scientificName"]) & FileInfo[i,"Column_scientificName"]!=""){
      col_spec_names <- FileInfo[i,"Column_scientificName"]
      if (is.na(col_spec_names)) stop(paste("Column with taxon names not found in",FileInfo[i,"Dataset_brief_name"],"file!"))
      all_column_names <- col_spec_names
    }

    ## Column location
    col_reg_names <- FileInfo[i,"Column_location"]
    if (is.na(col_reg_names)) stop(paste("Column with location names not found in",FileInfo[i,"Dataset_brief_name"],"file!"))
    all_column_names <- c(all_column_names,col_reg_names)

    ## Column kingdom
    if (!is.na(FileInfo[i,"Column_kingdom"]) & FileInfo[i,"Column_kingdom"]!=""){
      col_kingdom <- FileInfo[i,"Column_kingdom"]
      all_column_names <- c(all_column_names,col_kingdom)
    }

    ## Column ISO codes of countries
    if (!is.na(FileInfo[i,"Column_country_ISO"]) & FileInfo[i,"Column_country_ISO"]!=""){
      col_country_code <- FileInfo[i,"Column_country_ISO"]
      all_column_names <- c(all_column_names,col_country_code)
    }
    
    ## column eventDate 1
    if (!is.na(FileInfo[i,"Column_eventDate1"]) & FileInfo[i,"Column_eventDate1"]!=""){
      col_eventDate_1 <- FileInfo[i,"Column_eventDate1"]
      all_column_names <- c(all_column_names,col_eventDate_1)
    }
    
    ## column eventDate 2
    if (!is.na(FileInfo[i,"Column_eventDate2"]) & FileInfo[i,"Column_eventDate2"]!=""){
      col_eventDate_2 <- FileInfo[i,"Column_eventDate2"]
      all_column_names <- c(all_column_names,col_eventDate_2)
    }
    
    ## column establishmentMeans
    if (!is.na(FileInfo[i,"Column_establishmentMeans"]) & FileInfo[i,"Column_establishmentMeans"]!=""){
      col_establishmentMeans <- FileInfo[i,"Column_establishmentMeans"]
      colnames(dat)[col_names_import==col_establishmentMeans] <- "establishmentMeans"
      all_column_names <- c(all_column_names,"establishmentMeans")
      dat$establishmentMeans <- tolower(dat$establishmentMeans)
    }
    
    ## column occurrenceStatus
    if (!is.na(FileInfo[i,"Column_occurrenceStatus"]) & FileInfo[i,"Column_occurrenceStatus"]!=""){
      col_occurrenceStatus <- FileInfo[i,"Column_occurrenceStatus"]
      if (col_establishmentMeans==col_occurrenceStatus){ # check if same column has been assigned before in establishmentMeans
        dat$occurrenceStatus <- dat$establishmentMeans # if yes, duplicate column
      } else {
        colnames(dat)[col_names_import==col_occurrenceStatus] <- "occurrenceStatus"
      }
      all_column_names <- c(all_column_names,"occurrenceStatus")
      dat$occurrenceStatus <- tolower(dat$occurrenceStatus)
    }
    
    ## column degreeOfEstablishment
    if (!is.na(FileInfo[i,"Column_degreeOfEstablishment"]) & FileInfo[i,"Column_degreeOfEstablishment"]!=""){
      col_degreeOfEstablishment <- FileInfo[i,"Column_degreeOfEstablishment"]
      if (col_establishmentMeans==col_degreeOfEstablishment){ # check if same column has been assigned before in establishmentMeans
        dat$degreeOfEstablishment <- dat$establishmentMeans # if yes, duplicate column
      } else if (col_establishmentMeans==col_occurrenceStatus){
        dat$degreeOfEstablishment <- dat$occurrenceStatus # if yes, duplicate column
      } else {
        colnames(dat)[col_names_import==col_degreeOfEstablishment] <- "degreeOfEstablishment"
      }
      all_column_names <- c(all_column_names,"degreeOfEstablishment")
      dat$degreeOfEstablishment <- tolower(dat$degreeOfEstablishment)
    }
    
    ## column pathway
    if (!is.na(FileInfo[i,"Column_pathway"]) & FileInfo[i,"Column_pathway"]!=""){
      col_pathway <- FileInfo[i,"Column_pathway"]
      all_column_names <- c(all_column_names,col_pathway)
    }
    
    ## column habitat
    if (!is.na(FileInfo[i,"Column_habitat"]) & FileInfo[i,"Column_habitat"]!=""){
      col_habitat <- FileInfo[i,"Column_habitat"]
      all_column_names <- c(all_column_names,col_habitat)
    }
    
    ## column bibliographicCitation
    if (!is.na(FileInfo[i,"Column_bibliographicCitation"]) & FileInfo[i,"Column_bibliographicCitation"]!=""){
      col_bibliographicCitation <- FileInfo[i,"Column_bibliographicCitation"]
      all_column_names <- c(all_column_names,col_bibliographicCitation)
    }

    ## additional columns to keep throughout the following processing
    if (!is.na(FileInfo[i,"Column_additional"]) & FileInfo[i,"Column_additional"]!=""){
      col_additional <- FileInfo[i,"Column_additional"]
      addit_cols <- unlist(strsplit(col_additional,"; "))
      all_column_names <- c(all_column_names,colnames(dat)[pmatch(addit_cols,colnames(dat))])
    }
    
    ## keep required, optional and additional columns
    dat_out <- dat[,all_column_names]
    dat_out[is.null(dat_out)] <- ""
    dat_out[is.na(dat_out)] <- ""
    
    ## standardise column names
    col_names_import <- colnames(dat_out)
    if (exists("col_spec_names")) colnames(dat_out)[col_names_import==col_spec_names] <- "verbatimTaxonRank"
    if (exists("col_reg_names")) colnames(dat_out)[col_names_import==col_reg_names] <- "location_orig"
    if (exists("col_kingdom")) colnames(dat_out)[col_names_import==col_kingdom] <- "Kingdom_user"
    if (exists("col_country_code")) colnames(dat_out)[col_names_import==col_country_code] <- "Country_ISO"
    if (exists("col_eventDate_1")) colnames(dat_out)[col_names_import==col_eventDate_1] <- "eventDate"
    if (exists("col_eventDate_2")) colnames(dat_out)[col_names_import==col_eventDate_2] <- "eventDate2"
    if (exists("col_habitat")) colnames(dat_out)[col_names_import==col_habitat] <- "habitat"
    if (exists("col_pathway")) colnames(dat_out)[col_names_import==col_pathway] <- "pathway"
    if (exists("col_bibliographicCitation")) colnames(dat_out)[col_names_import==col_bibliographicCitation] <- "bibliographicCitation"
    
    if (exists("col_habitat")) dat$habitat <- tolower(dat$habitat) # for easier matching later

    
    ## Prepare output and clean ################################################
    
    options(warn=-1)
    rm(col_spec_names,col_reg_names,col_kingdom,col_country_code,col_eventDate_1,
       col_eventDate_2,col_establishmentMeans,col_occurrenceStatus,
       col_habitat,col_bibliographicCitation,col_establishmentMeans,col_occurrenceStatus,
       col_degreeOfEstablishment,col_pathway)
    options(warn=1)
    
    ## remove rows with missing taxon and region names
    dat_out <- dat_out[!dat_out$location_orig=="",]
    dat_out <- dat_out[!dat_out$verbatimTaxonRank=="",]
    
    dat_out$Taxon_group <- FileInfo[i,"Taxon_group"]
    
    colnames(dat_out) <- gsub("\\.+","_",colnames(dat_out))
    dat_out$verbatimTaxonRank <- gsub("\"","",dat_out$verbatimTaxonRank) # remove additional quotes to avoid difficulties with export
    dat_out$verbatimTaxonRank <- gsub("\\\\","",dat_out$verbatimTaxonRank) # remove back shlashes

    dat_out <- unique(dat_out) # remove duplicates
    
    # output mis-matches for reference
    write.table(dat_out,file.path("Output","Intermediate",paste("Step1_StandardColumns_",FileInfo[i,"Dataset_brief_name"],".csv",sep="")),row.names = F)
  }
}
