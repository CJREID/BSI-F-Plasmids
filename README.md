# BSI-F-Plasmids
Data analysis and visualisation for manuscript "Characterisation of major F plasmid clusters in bloodstream isolates of Escherichia coli" by Reid C.J., Cummins M.L. and Djordjevic S.P.

## Overview
The following repository allows conscientious readers of the manuscript to reproduce the data processing, statistics and figures presented in the paper.

It comprises two directories: __`scripts`__ and __`data`__ (the contents of which should be self explanatory) and generates a __`data/outputs`__ folder with __`figures`__ and __`data`__ subdirectories.


## Installation
### Software requirements
These scripts are currently functional on mac OS Sonoma 14.5 using RStudio 2023.12.1+402 and R version 4.3.2. We cannot guarantee they will work on other distributions of R or RStudio, but they should.

### Packages required
You will need to install the following packages and versions to work with the scripts:
- tidytree_0.4.6
- abricateR_0.1.2
- dplyr_1.1.4
- pheatmap_1.0.12
- RColorBrewer_1.1-3
- purrr_1.0.2
- data.table_1.15.0
- ggpubr_0.6.0
- readr_2.1.5
- magrittr_2.0.3
- reshape2_1.4.4
- tidyr_1.3.1
- scales_1.3.0
- lemon_0.4.7
- tibble_3.2.1
- ggplot2_3.4.4
- forcats_1.0.0
- tidyverse_2.0.0
- plasmidmapR_0.1.0
- stringr_1.5.1


### plasmidmapR and abricateR
Please see [plasmidmapR](https://github.com/maxlcummins/plasmidmapR) and [abricateR](https://github.com/maxlcummins/abricateR) repositories to install these packages.

## Usage
Clone this repository
```
git clone https://github.com/CJREID/BSI-F-Plasmids.git
cd BSI-F-Plasmids
pwd
```
Open the BSI_F_Plasmids.R script in a text editor or RStudio and set the variable __`wrkdir`__ on line 16 to the output of __`pwd`__ above and save the script.

Run the BSI_F_Plasmids.R script and your __`outputs`__ folder will populate with figures and tables form the publication.
