## SInAS workflow (vs. 2.0): Integration and standardisation of alien species occurrence data
Author: Hanno Seebens, Manuela Gómez-Suárez, Giessen-Frankfurt, 25.06.2025




The R scripts contain the implementation of a workflow to standardise and integrate 
databases of alien species occurrences and years of first record. The workflow is 
described in detail in the manual 'Manual SInAS workflow.pdf' and in Seebens, H., 
D. A. Clarke, Q. Groom, J. R. U. Wilson, E. García-Berthou, I. Kühn, M. Roigé, 
S. Pagad, F. Essl, J. Vicente, M. Winter, and M. McGeoch. 2020. A workflow for 
standardising and integrating alien species distribution data. NeoBiota 59:39–59.

The R scripts can be used and modified freely as long as the work is properly cited
with the aforementioned citation of the scientific paper.

In brief, the workflow consists of the following steps:
1 : Preparation of column names of alien taxon databases  
2a: Standardisation of terminologies
2b: Standardisation of location names
2c: Standardisation of taxon names
2d: Standardisation of event dates 
3 : Merging all standardised databases

Input: 
Information about databases has to be provided in DatabaseInfo.xlsx.
The workflow can be adjusted by modifying input tables providing information about
location names, taxon names, terminologies and rules to treat first records can be
done in UserDefinedTaxonNames.xlsx, AllLocations.xlsx, 
Guidelines_eventDates.xlsx and five translation tables for pathway, habitat, 
occurrence status, degree of establishment and establishment means, respectively.
Note that only the first sheet of the Excel file is read in. Others are ignored.

Output: 
A standardised merged database built from all databases and a full list of taxon names with 
further taxonomic information will be exported.
Several data sets will be exported in addition by the workflow depending on the degree of 
matching e.g., missing location names, unresolved terms, missing taxon names. These 
files can be used for cross-checking and further refinement of the original databases
and the translation tables.

Please consult the manual and the scientific paper for further reference.

Pleaes cite the scientific paper when using the workflow.

