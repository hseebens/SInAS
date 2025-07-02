## SInAS workflow (vs. 3.0): Integration and standardisation of alien species occurrence data
Author: Hanno Seebens, Manuela Gómez-Suárez, Giessen-Frankfurt, 01.07.2025


The workflow of Standardising and Integrating Alien Species distributional data (SInAS) 
has been developed to harmonise distributional data of alien species including their native
ranges. The workflows has the capability to handle data provided as regional lists (i.e.
checklists) of alien species and is fexible enough to adjust the spatial delineation 
according to the user's need. Original records will be assigned to the selected 
geographic classification if possible. If not, lists of mis-matches will be provided. Information
about taxa, locations, invasion status, habitat and year of first records will be standardised
while applying the workflow following Darwin Core terminology.

In brief, the workflow consists of the following steps:
1 : Preparation of column names of alien taxon databases  
2a: Standardisation of terminologies
2b: Standardisation of location names
2c: Standardisation of taxon names
2d: Standardisation of event dates 
3 : Merging all standardised databases

Input: 
Information about databases to be integrated has to be provided in a table format.
The workflow can be adjusted by modifying configuration files (e.g. translation tables). These
files provide the user the opportunity to adjust the configurations of standardisation 
according to the need of the respective analysis.

Output: 
A standardised merged dataset built from all provided databases and a full list of taxon names with 
further taxonomic information will be exported.
Several data sets will be exported in addition by the workflow depending on the degree of 
matching e.g. missing location names, unresolved terms or missing taxon names. These 
files can be used for cross-checking and further refinement of the original databases
and the translation tables.

Please consult the manual (Manual SInAS workflow.pdf) and the following scientific papers 
for further information. Pleaes cite the respestive scientific paper when using the workflow.

- Gómez-Suárez, M., Laeseke, P., and Seebens, H. (submitted) A global dataset of native and
alien distributions of alien species

- Seebens, H., D. A. Clarke, Q. Groom, J. R. U. Wilson, E. García-Berthou, I. Kühn, M. Roigé, 
S. Pagad, F. Essl, J. Vicente, M. Winter, and M. McGeoch. 2020. A workflow for 
standardising and integrating alien species distribution data. NeoBiota 59:39–59.



The R scripts can be used and modified freely as long as the work is properly cited
with the aforementioned citation of the scientific paper.
