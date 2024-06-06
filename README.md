# BSI-F-Plasmids
Data analysis and visualisation for manuscript "F plasmids in Escherichia coli bloodstream infections" by Reid C.J., Cummins M.L. and Djordjevic S.P.

## Overview
The following repository allows conscientious readers of the manuscript to reproduce the data processing, statistics and figures presented in the paper.

It comprises two directories: __`scripts`__ and __`data`__ (the contents of which should be self explanatory) and generates a __`data/outputs`__ folder with __`figures`__ and __`data`__ subdirectories.


## Installation
### Software requirements
These scripts are currently functional on mac OS Ventura 13.2.1 using RStudio 2022.07.1+554 and R version 4.2.2. We cannot guarantee they will work on other distributions of R or RStudio, but they should.

### Packages required
You will need to install the following packages and versions to work with the scripts:
- rstatix_0.7.0     
- plasmidmapR_0.1.0  
- abricateR_0.1.1    
- RColorBrewer_1.1-3
- ggpubr_0.4.0       
- reshape2_1.4.4     
- lemon_0.4.5        
- forcats_0.5.2      
- stringr_1.4.1     
- dplyr_1.0.10     
- purrr_0.3.4        
- readr_2.1.2        
- tidyr_1.2.1        
- tibble_3.1.8       
- ggplot2_3.3.6      
- tidyverse_1.3.2  

### plasmidmapR and abricateR
Please see [plasmidmapR](https://github.com/maxlcummins/plasmidmapR) and [abricateR](https://github.com/maxlcummins/abricateR) repositories to install these packages.

## Usage
Clone this repository
```
git clone [https://github.com/CJREID/ST372.git](https://github.com/CJREID/BSI-F-Plasmids.git)
cd BSI-F-Plasmids
pwd
```
Open the BSI_F_Plasmids.R script in a text editor or RStudio and set the variable __`wrkdir`__ on line 16 to the output of __`pwd`__ above and save the script.

Run the BSI_F_Plasmids.R script and your __`outputs`__ folder will populate with figures and tables.

## Outputs
### Figures
1. Figure 1. F plasmid groups and RSTs in the collection
2. Figure 2. Counts of F RSTs by F plasmid group for the top ten E. coli STs
3. Figure 3. Relationship between HPI carriage, ColV carriage and the presence of iron acquisition systems stacked by E. coli ST
4. Figure 4. Comparison of total ARG and VAG carriage across F plasmid groups, ST frequency and common STs

### Supplementary
#### Tables
1. Table S1. Metadata, accession numbers and gene screening results for 4711 *E. coli* genomes used in the study

#### Figures
1. Fig S1. F plasmids groups in the collection
2. Fig S2. Proportional distribution of phylogroups by ST frequency and F plasmid group
3. Fig S3. F RSTs in the collection
4. Fig S4. Counts of F RSTs by F plasmid group for other common E. coli STs
