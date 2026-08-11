
####### SInAS workflow: Integration and standardisation of alien species data ###########
##
## Step 2b: Standardisation of location names using a translation table "AllLocations.xlsx"
##
## Manuela Gómez-Suárez, Hanno Seebens, Giessen, 17.04.2026
#########################################################################################


StandardiseLocationNames <- function(FileInfo = NULL) {
  
  ## 1. Load translation tables ################################################
  regions_raw <- as.data.table(
    read.xlsx(file.path("Config", "AllLocations.xlsx"), sheet = "location", na.strings = "")
  )[, .(locationID, location, location_var)]
  regions_raw[, location_var   := tolower(location_var)]
  regions_raw[, location_lower := tolower(location)]
  
  subregions_raw <- as.data.table(
    read.xlsx(file.path("Config", "AllLocations.xlsx"), sheet = "stateProvince", na.strings = "")
  )[, .(locationID, location, gadm1_name, gadm1_var)]
  subregions_raw[, gadm1_var   := tolower(gadm1_var)]
  subregions_raw[, gadm1_lower := tolower(gadm1_name)]
  
  ## Add compound name variants: "gadm1_name, country" and "gadm1_name (country)"
  subregions_raw[, combined_name    := tolower(paste(gadm1_name, location, sep = ", "))]
  subregions_raw[, parenthesis_name := tolower(paste0(gadm1_name, " (", location, ")"))]
  subregions_raw[, gadm1_var := fcase(
    is.na(gadm1_var) | gadm1_var == "",
    paste(combined_name, parenthesis_name, sep = "; "),
    default = paste(gadm1_var, combined_name, parenthesis_name, sep = "; ")
  )]
  subregions_raw[, c("combined_name", "parenthesis_name") := NULL]
  
  ## Detect ambiguous subregion base names (shared across countries)
  base_names <- gsub("\\s*\\(.*?\\)", "", subregions_raw$gadm1_name)
  dup <- unique(base_names[duplicated(base_names)])
  
  ## 2. Build flat key → value lookup tables ###################################
  ## One key for each location variation
  expand_keys <- function(dt, key_col, value_cols) {
    keys_list <- strsplit(dt[[key_col]], ";", fixed = TRUE)
    lens      <- lengths(keys_list)
    out       <- dt[rep(seq_len(nrow(dt)), lens), value_cols, with = FALSE]
    out[, key := trimws(unlist(keys_list))]
    out[nzchar(key)]          # drop any empty strings
  }
  
  ## Regions lookup: direct country names + location_var variants
  regions_lut <- unique(
    rbind(
      regions_raw[, .(key = location_lower, locationID, location)],
      regions_raw[!is.na(location_var),
                  expand_keys(.SD, "location_var", c("locationID", "location"))]
    ),
    by = "key"
  )
  setkey(regions_lut, key)
  
  ## Subregions lookup: direct gadm1 names + gadm1_var variants
  subregions_lut <- unique(
    rbind(
      subregions_raw[, .(key = gadm1_lower, locationID, location, gadm1_name)],
      subregions_raw[nzchar(gadm1_var),
                     expand_keys(.SD, "gadm1_var", c("locationID", "location", "gadm1_name"))]
    ),
    by = "key"
  )
  setkey(subregions_lut, key)
  
  ## Country-level table for final locationID re-assignment
  ## (locationID always reflects country level, not subregion level)
  regions_clean <- unique(regions_raw[, .(location, locationID)], by = "location")
  setkey(regions_clean, location)
  
  ## 3. Input file list ########################################################
  inputfiles <- paste0("Step2_StandardTerms_",
                       na.omit(FileInfo$Dataset_brief_name),
                       ".csv")
  
  ## 4. Process files in parallel ##############################################
  n_cores <- min(length(inputfiles), max(1L, detectCores() - 1L))
  if (n_cores > 1L) {
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
    on.exit(stopCluster(cl), add = TRUE)
    
    
    clusterExport(
      cl,
      c("regions_lut", "subregions_lut", "regions_clean", 
        "dup", "FileInfo", "inputfiles"),
      envir = environment()
    )
    
    `%exec%` <- `%dopar%`
  } else {
    `%exec%` <- `%do%`
  }

  
  foreach(
    i         = seq_along(inputfiles),
    .packages = "data.table"
  ) %exec% {
    
    dat <- as.data.table(
      read.table(file.path("Output", "Intermediate", inputfiles[i]),
                 header = TRUE, stringsAsFactors = FALSE)
    )
    
    ## Normalise the location string ###########################################
    dat[, loc_key := tolower(location_orig)]
    dat[, loc_key := gsub("\u00a0|\u00c2", " ", loc_key)]
    dat[, loc_key := trimws(loc_key)]
    dat[, loc_key := gsub("  +", " ",  loc_key)]
    dat[, loc_key := gsub(" \\(the\\)$", "", loc_key)]
    
    ## Initialise output columns (ensures clean state regardless of input)
    dat[, c("locationID", "location", "stateProvince") :=
          .(NA_integer_, NA_character_, NA_character_)]
    
    ## Two-pass keyed join #####################################################
    ## Pass 1: subregion match — sets stateProvince, location, locationID
    dat[subregions_lut, on = .(loc_key = key),
        `:=`(locationID    = i.locationID,
             location      = i.location,
             stateProvince = i.gadm1_name)]
    
    ## Pass 2: region match only for rows still unmatched 
    dat[is.na(locationID),
        `:=`(locationID = regions_lut[.SD, on = .(key = loc_key), x.locationID],
             location   = regions_lut[.SD, on = .(key = loc_key), x.location])]
    
    ## ── Warn about ambiguous non-matched subregion names #####################
    ambig <- unique(dat$location_orig[dat$location_orig %in% dup])
    if (length(ambig) > 0L)
      message(sprintf(
        "\n    Warning: Unresolved terms in %s. Location name(s) matching multiple subregions: %s.\n    Add country in parentheses, e.g. 'Amazonas (Colombia)'.\n",
        FileInfo[i, "Dataset_brief_name"],
        paste(ambig, collapse = ", ")
      ))
    
    ## Remove duplicates #######################################################
    dat <- unique(dat)
    
    ## Re-assign locationID from country-level table ###########################
    dat[, locationID := NULL]
    dat[regions_clean, on = "location", locationID := i.locationID]
    
    ## Write output and drop working columns ###################################
    drop_cols <- intersect(c("stateProvince", "loc_key"), names(dat))
    write.table(
      dat[, -drop_cols, with = FALSE],
      file.path("Output", "Intermediate",
                paste0("Step3_StandardLocationNames_",
                       FileInfo[i, "Dataset_brief_name"], ".csv")),
      row.names = FALSE
    )
    
    ## Export unmatched locations ##############################################
    missing_locs <- sort(unique(dat$location_orig[is.na(dat$locationID)]))
    if (length(missing_locs) > 0L)
      write.table(
        missing_locs,
        file.path("Output", "Check",
                  paste0("Missing_Locations_",
                         FileInfo[i, "Dataset_brief_name"], ".csv")),
        row.names = FALSE, col.names = FALSE
      )
    
    invisible(NULL)
  }
  
  ## 5. Post-process: aggregate changed location names #########################
  reg_names <- rbindlist(lapply(seq_along(inputfiles), function(i) {
    dat <- as.data.table(
      read.table(
        file.path("Output", "Intermediate",
                  paste0("Step3_StandardLocationNames_",
                         FileInfo[i, "Dataset_brief_name"], ".csv")),
        header = TRUE, stringsAsFactors = FALSE
      )
    )
    data.table(location      = dat$location,
               location_orig = dat$location_orig,
               datasetName   = FileInfo[i, 1L])
  }))
  
  reg_names <- reg_names[location != location_orig]
  reg_names <- unique(reg_names[order(location)])
  reg_names[regions_clean, on = "location", locationID := i.locationID]
  
  write.table(reg_names,
              file.path("Output", "Translated_LocationNames.csv"),
              row.names = FALSE)
}