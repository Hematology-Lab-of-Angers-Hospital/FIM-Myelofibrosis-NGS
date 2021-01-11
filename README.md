# Repository information :
# Genomic analysis of primary and secondary myelofibrosis redefines the prognostic relevance of ASXL1 mutations: a French Intergroup of Myeloproliferative neoplasms (FIM) study
## Author of code : Laboratory of Hematology, Angers Universitary Hospital





# Data/
## Database of variant
Data/Data_Variant_result_479MF.csv
## Classification of gene categories for circos
Groupe_gene.csv
## Gene of panel NGS
List_gene_panel.csv
## Code for figure of Genomics
 Code to generate figure for description of the molecular landscape



# Download
git clone https://github.com/Hematology-Lab-of-Angers-Hospital/FIM-Myelofibrosis-NGS

# Launch generation of figure
python3 Code/Figure_genomics.py


## Dependencies
## Software
python 3.6.9 (2020-10-8)
R 3.4.4 (2018-03-15)
# Module 
## R
### dplyr
install.packages("dplyr")
library(dplyr)
### circlize
install.packages("circlize") 
library(circlize)
### RColorBrewer
install.packages("RColorBrewer")
library("RColorBrewer")

## Launch with Python 3.6.9
os
# Treat Data
numpy==1.18.2
pandas==0.24.2
# Library for figure
matplot==lib3.2.1
seaborn==0.9.0

### To install python library
#### Requirement
pip install -r requirement.txt
