
####### SInAS workflow: Integration and standardisation of alien species data ###########
##
## Step 3: Merging of standardised databases of alien species data
##
## Manuela Gómez-Suárez, Hanno Seebens, Giessen, 25.09.2025
#########################################################################################


MergeDatabases <- function(FileInfo = NULL,
                           version = NULL,
                           outputfilename = NULL,
                           output = NULL) {
  
  ## 1. Create a list with the standardized file names #########################
  inputfiles <- paste0("Step5_StandardIntroYear_",
                       na.omit(FileInfo$Dataset_brief_name), ".csv")
  
  ## Columns to drop in final dataset
  drop_cols <- c("originalNameUsage", "location_orig", "Kingdom_user",
                 "Country_ISO", "ISO3", "eventDate_orig", "Taxon_group")
  
  ## Key columns to merge in final dataset
  by_cols <- c("location", "locationID", "taxon", "scientificName",
               "taxonID", "establishmentMeans")
  
  ## 2. Load, clean, and pre-process each file #################################
  datasets <- lapply(seq_along(inputfiles), function(i) {
    
    ## read file
    
    dat <- as.data.table(
      read.table(file.path("Output", "Intermediate", inputfiles[i]),
                 header = TRUE, stringsAsFactors = FALSE)
    )
    
    ## Drop unwanted columns and add source metadata
    dat[, (intersect(drop_cols, names(dat))) := NULL]
    dat[, datasetName           := FileInfo[i, 1L]]
    dat[, bibliographicCitation := FileInfo[i, "Column_bibliographicCitation"]]
    
    ## Drop records with no location match
    dat <- dat[!is.na(locationID)]
    
    ## Resolve absent-present conflicts WITHIN each source ##################
    if ("occurrenceStatus" %in% names(dat)) {
      
      resolve_abs <- function(sub_dt) {
        ## Returns rows to REMOVE: "absent" rows where a "present" also exists
        sub_dt[
          sub_dt[, .I[all(c("absent", "present") %in% occurrenceStatus)
                      & occurrenceStatus == "absent"],
                 by = .(location, locationID, taxonID, taxon,
                        establishmentMeans)]$V1
        ]
      }
      
      # resolve conflicts individually for introduced and for native records 
      
      abs_int <- dat[establishmentMeans %in% c("uncertain", "introduced"),
                     resolve_abs(.SD)]
      abs_nat <- dat[establishmentMeans == "native", resolve_abs(.SD)]
      
      remove_key <- c("location", "locationID", "taxonID", "taxon",
                      "establishmentMeans", "occurrenceStatus")
      if (nrow(abs_int) > 0L)
        dat <- dat[!abs_int, on = remove_key]
      if (nrow(abs_nat) > 0L)
        dat <- dat[!abs_nat, on = remove_key]
    }
    
    dat
  })
  
  ##  3. Stack all datasets at once ############################################
  alldat <- rbindlist(datasets, use.names = TRUE, fill = TRUE)
  
  ## 4. Collapse introduced + uncertain records into one group #################
  alldat[, dist_group := fcase(
    establishmentMeans %in% c("introduced", "uncertain"), "a",
    establishmentMeans == "native",                        "b",
    establishmentMeans == "vagrant",                       "c",
    default = NA_character_
  )]
  
  alldat[, eventDate := as.character(eventDate)]
  
  ## Concatenate non-key columns within each group, then restore establishmentMeans
  non_key <- setdiff(names(alldat), c(by_cols, "dist_group"))
  
  alldat <- alldat[,
                   c(.(establishmentMeans = paste(sort(unique(establishmentMeans)), collapse = "; ")),
                     lapply(.SD, function(x) paste(unique(x), collapse = "; "))),
                   by = .(location, locationID, taxon, scientificName, taxonID, dist_group),
                   .SDcols = non_key
  ]
  
  ## 5. Row merge: collapse duplicated taxon-location-status rows ##############
  alldat[, occurrenceStatus := ifelse(
    occurrenceStatus == "absent; present", "present", occurrenceStatus
  )]
  alldat[, occurrenceStatus := ifelse(
    occurrenceStatus == "present; absent", "present", occurrenceStatus
  )]
  
  alldat <- alldat[,
                   lapply(.SD, function(x) paste(x, collapse = "; ")),
                   by = .(taxon, location, locationID, taxonID, scientificName, establishmentMeans)
  ]
  
  ## 6. Clean all non-key columns to remove dup. strings #######################
  key_set  <- c("taxon", "taxonID", "location", "locationID",
                "scientificName", "establishmentMeans")
  val_cols <- setdiff(names(alldat), key_set)
  
  ## Convert value columns to a matrix for fast vectorised gsub
  m <- as.matrix(alldat[, val_cols, with = FALSE])
  
  ## Sort and remove duplicates between semicolons
  m[] <- vapply(m, function(cell) {
    parts <- unique(trimws(strsplit(cell, "; ", fixed = TRUE)[[1L]]))
    parts <- parts[nzchar(parts) & parts != "NA"]
    if (length(parts) == 0L) "" else paste(sort(parts), collapse = "; ")
  }, character(1L))
  
  alldat[, (val_cols) := as.data.table(m)]
  
  ## 7. Select earliest eventDate where multiple dates are available ###########
  multi_date <- grepl("; ", alldat$eventDate, fixed = TRUE)
  if (any(multi_date)) {
    alldat[multi_date, eventDate := vapply(
      strsplit(eventDate, "; ", fixed = TRUE),
      function(s) as.character(min(as.integer(s), na.rm = TRUE)),
      character(1L)
    )]
  }
  
  ## 8. Resolve duplicated rows ################################################
  alldat <- unique(alldat)
  
  ## 9. Resolve remaining taxon-location-establishmentMeans duplicates #########
  dup_key <- intersect(
    c("taxon", "location", "establishmentMeans",
      "taxonID", "locationID", "scientificName"),
    names(alldat)
  )
  merge_cols <- setdiff(names(alldat), dup_key)
  
  alldat <- alldat[,
                   {
                     if (.N == 1L) {
                       .SD
                     } else {
                       res <- lapply(names(.SD), function(cn) {
                         vals <- .SD[[cn]]
                         non_na <- vals[!is.na(vals) & vals != ""]
                         if (length(non_na) == 0L) return(NA_character_)
                         if (cn == "eventDate") {
                           as.character(min(as.integer(non_na), na.rm = TRUE))
                         } else {
                           paste(sort(unique(non_na)), collapse = "; ")
                         }
                       })
                       setNames(as.data.table(res), names(.SD))
                     }
                   },
                   by = dup_key,
                   .SDcols = merge_cols
  ]
  
  ## 10. Fill habitat per taxon ################################################
  alldat[, habitat := habitat[which.max(nchar(habitat, keepNA = FALSE))],
         by = taxonID]
  
  ## 11. Join taxonomic tree and sort ##########################################
  fulltaxalist <- as.data.table(
    read.table(
      file.path("Output", paste0(outputfilename, "_", version, "_FullTaxaList.csv")),
      header = TRUE, stringsAsFactors = FALSE
    )
  )[, .(taxonID, kingdom, phylum, class, order, family)]
  fulltaxalist <- unique(fulltaxalist)
  
  alldat <- fulltaxalist[alldat, on = "taxonID"]
  setorder(alldat, location, kingdom, phylum, class, order, family, scientificName)
  
  ##  12. Select and reorder output columns ####################################
  all_addit_cols <- unique(unlist(strsplit(
    paste(FileInfo$Column_additional, collapse = "; "), "; "
  )))
  all_addit_cols <- all_addit_cols[all_addit_cols != "NA"]
  
  key_columns <- c("location", "locationID", "taxon", "scientificName",
                   "taxonID", "eventDate", "habitat", "occurrenceStatus",
                   "establishmentMeans", "degreeOfEstablishment",
                   "pathway", "datasetName", "bibliographicCitation")
  
  columns_out <- intersect(c(key_columns, all_addit_cols), names(alldat))
  
  alldat_out <- alldat[, columns_out, with = FALSE]
  alldat_out[is.na(alldat_out)] <- ""
  alldat_out[alldat_out == "NA"] <- ""
  
  setnames(alldat_out, "FirstRecord_orig", "eventDate_orig",
           skip_absent = TRUE)
  
  alldat_out[, scientificName := NULL]
  
  ## 13. Final cleanup of citation and write output ############################
  alldat_out[, bibliographicCitation :=
               gsub('["\']', '', bibliographicCitation)]
  
  write.table(alldat_out,
              file.path("Output", paste0(outputfilename, "_", version, ".csv")),
              row.names = FALSE)
  
  cat(paste("\n Successfully established version", version,
            "of", outputfilename, "file. \n"))
  
  ## 14. Optionally remove intermediate files ##################################
  if (!output) {
    interm_res <- list.files(file.path("Output", "Intermediate"))
    interm_res <- grep("Standard", interm_res, value = TRUE)
    file.remove(file.path("Output", "Intermediate", interm_res))
  }
}
