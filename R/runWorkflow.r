####### SInAS workflow: Integration and standardisation of alien species data ###########
##
## This script executes the SInAS workflow to standardise and integrate alien species 
## occurrence data. All required files and a manual can be found on GitHub:
## https://github.com/hseebens/SInAS
##
## Run the whole workflow by e.g. copy-pasting "source("runWorkflow.r")" into the R terminal.
##
## Manuela Gómez-Suárez, Hanno Seebens, Gießen, 02.07.2025
#########################################################################################

graphics.off()
rm(list=ls())

## required libraries
library(rgbif) # for checking names, records and taxonomy; note: usage of rgbif may cause warnings like "Unknown or uninitalised column: " which is a bug. Can be ignored.
library(openxlsx)
library(data.table)
library(tidyverse)

## option for storing the intermediate and final output
outputfilename <- "SInAS" # name of final output file

version <- "3.0" # which version of the database are you going to produce? this will be attached to the end of 'outputfilename'

output <- T # shall intermediate results be stored to disk? (may overwrite existing files!)

################################################################################
### load functions #############################################################
source(file.path("R","PrepareDatasets.r")) # preparing example data sets as input files
source(file.path("R","StandardiseTaxonNames.r")) # standardising taxon names, requires GBIF connection, takes some time...
source(file.path("R","OverwriteTaxonNames.r")) # replace taxon names with user-defined ones
source(file.path("R","StandardiseLocationNames.r")) # standardising location names
source(file.path("R","StandardiseTerms.r")) # standardising location names
source(file.path("R","GeteventDate.r")) # standardising location names
source(file.path("R","MergeDatabases.r")) # combine data sets
source(file.path("R","CheckGBIFTax.r")) #function to check taxon names using GBIF taxonomy
################################################################################


################################################################################
######## Load data set table ###################################################

FileInfo <- read.xlsx(file.path("Config","DatabaseInfo.xlsx"),sheet=1)
if (nrow(FileInfo)==0) stop("No database information provided. Add information to Config/DatabaseInfo.xlsx.")

## load databases, extract required information and harmonise taxon names...
cat("\n Step 1 Preparation of provided data sets \n")
PrepareDatasets(FileInfo)

## load databases, extract required information and harmonise taxon names...
cat("\n Step 2 Standardisation of terminology \n")
StandardiseTerms(FileInfo)

## harmonise taxon names...
cat("\n Step 3 Standardisation of location names \n")
StandardiseLocationNames(FileInfo)

## load databases, extract required information and harmonise taxon names...
cat("\n Step 4 Standardisation of taxon names \n")
StandardiseTaxonNames(FileInfo)
OverwriteTaxonNames(FileInfo) # user-defined taxon names

## standardise first records....
cat("\n Step 5 Standardisation of eventDate \n")
GeteventDate(FileInfo)

## merge databases...
cat("\n Step 6 Merging databases \n")
MergeDatabases(FileInfo,version,outputfilename,output)
