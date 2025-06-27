#########################################################################################
####### SInAS workflow: Integration and standardisation of alien species data ###########
#########################################################################################
##
## Hanno Seebens, Frankfurt, 25.06.2025
#########################################################################################

The R scripts contain the implementation of a workflow to standardise and merge 
databases of alien species occurrences and years of first record. The workflow is 
described in detail in Seebens, H., 
D. A. Clarke, Q. Groom, J. R. U. Wilson, E. García-Berthou, I. Kühn, M. Roigé, 
S. Pagad, F. Essl, J. Vicente, M. Winter, and M. McGeoch. 2020. A workflow for 
standardising and integrating alien species distribution data. NeoBiota 59:39–59.

The R scripts can be used and modified freely as long as the work is properly cited
with the aforementioned citation.

In brief, the workflow consists of the following five steps:
1. : Preparation of column names of alien taxon databases  
2a.: Standardisation of terminologies
2b.: Standardisation of location names
2c.: Standardisation of taxon names
2d.: Standardisation of event dates 
3. : Merging all standardised databases

Input: 
Information about databases has to be provided in DatabaseInfo.xlsx.
Modification of location names, taxon names, terminologies and rules to treat first records can be
done in UserDefinedSpeciesNames.xlsx, AllLocations.xlsx,
Guidelines_eventDates.xlsx and five translation tables for pathway, habitat, 
occurrence status, degree of establishment and establishment means.
Note that only the first sheet of the Excel file is read in. Others are ignored.

Output: 
A standardised masterfile built from all databases and a full list of taxon names with 
further taxonomic information.
Several data sets will be exported by the workflow depending on the degree of 
matching e.g., missing location names, unresolved terms, missing taxon names. These 
files can be used for cross-checking and further refinement of the original databases
and the translation tables.

