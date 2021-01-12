# Repository information :
# Genomic analysis of primary and secondary myelofibrosis redefines the prognostic relevance of ASXL1 mutations: a French Intergroup of Myeloproliferative neoplasms (FIM) study
#### Author of code : Laboratory of Hematology, Angers Universitary Hospital

Luque Paz D. et al. Blood Advances. in press

The folder contains the clinical and molecular data of the cohort.
## :question: Contacts
For any question about the data, feel free to contact Valérie Ugo (valerie.ugo@chu-angers.fr), Damien Luque Paz (damien.luquepaz@univ-angers.fr), and Maxime Renard for Bioinformatic part (maxime.renard@chu-angers.fr). Code for statistics can be share upon request.

## :open_file_folder: Data/
### Database of all mutations in the cohort
Data_Variants.csv 

Classification of mutation:
 A=Pathogenic
 B=Likely pathogenic
 C=Variant of unknow significance
### Clinical data of the 479 patients of the cohort
Data_Cohort.csv

Follow-up is reported as days
### Classification of gene categories for circos
Groupe_gene.csv
### Gene of the NGS panel
List_gene_panel.csv
## :open_file_folder: Code/
 Code to generate figures for description of the molecular landscape



## Download
```shell
git clone https://github.com/Hematology-Lab-of-Angers-Hospital/FIM-Myelofibrosis-NGS
```
## Launch generation of figure
```python
python3 Code/Figure_genomics.py
```
## Dependencies
python 3.6.9 (2020-10-8)
R 3.4.4 (2018-03-15)
### R
-dplyr

install.packages("dplyr")

library(dplyr)

-circlize

install.packages("circlize") 

library(circlize)

-RColorBrewer

install.packages("RColorBrewer")

library("RColorBrewer")

## Versionning
numpy==1.18.2

pandas==0.24.2

matplot==lib3.2.1

seaborn==0.9.0

### To install python library
#### Requirement
```shell
pip install -r requirement.txt
```
