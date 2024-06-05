####-----PACKAGES-----####
# Load required packages
library(tidyverse)
library(ggplot2)
library(lemon)
library(reshape2)
library(ggpubr)
library(RColorBrewer)
library(abricateR)
library(plasmidmapR)
library(rstatix)
library(ggforce)
library(scales)

####-----SETUP-----####
# SET WORKING DIRECTORY
## YOU NEED TO CHANGE THIS TO THE DIRECTORY YOU HAVE CLONED THE GITHUB REPOSITORY TO
wrkdir <- "/Volumes/126451/WORK/Projects/Sepsis_Plasmids/BSI_F_Plasmids/"
setwd(wrkdir)

# Define 'not in' function for subsetting
'%notin%' <- Negate('%in%')

# Make output directory
if (dir.exists("outputs")){
} else {
  dir.create("outputs")
  dir.create("outputs/data")
  dir.create("outputs/figures")
}
####-----FUNCTIONS-----####
# Function to quick display a custom pallette
see_palette <- function(palette){
  require(ggplot2)
  palette_df <- data.frame(
    labels = names(palette),
    colors = as.character(palette),
    n = rep(1/length(palette), length(palette))
  )
  
  # Plot
  ggplot(palette_df, aes(x = "Palette", y = n, fill =colors)) +
    geom_col() +
    geom_text(aes(label = labels), position = position_stack(vjust = 0.5)) +
    labs(fill = "labels") +
    scale_fill_identity() +
    scale_y_continuous(expand = c(0, 0)) +
    theme_classic() +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank())
}

####-----IMPORT DATA-----####
#Get input filenames
infiles <-c(list.files("data/full_summaries/", pattern = "\\.txt"))

#Read in and name data.frames with whatever comes before .txt
for (file in infiles){
  indir <- c("data/full_summaries")
  f <- read_delim(paste(indir, file, sep = "/"), delim = "\t", col_names = TRUE, trim_ws = TRUE)
  assign(paste(substr(file, 1, nchar(file)-4), sep = ""), f)
}

# Filter sequences with expected genome size for E. coli
assembly_ok <- assembly_stats %>% filter(total_length >= 4500000, total_length <= 6500000)

# Get rid of genomes with no ST
mlst_ok <- mlst %>% mutate(ST = gsub("\\s.*", "", ST)) %>% filter(ST != "-",scheme == "ecoli")

# Join to only include genomes of normal size with legitimate ST in the analysis
meta <- mlst_ok %>% filter(name %in% assembly_ok$name) %>% select(Name = name, ST)

####-----PROCESS ABRICATE DATA-----####
#Load paths to files needed for abricateR
abricate_path <- "data/full_summaries/genotype.txt"
pointfinder_path <- "data/full_summaries/pointfinder.txt"
pMLST_data <- "data/full_summaries/pMLST.txt"

# Provide output names
output_name <- "sepsis"
output_dir <- "processed_R"

# Run abricateR 
abricateR::abricateR(
  abricate_in = abricate_path,
  output = output_name,
  identity = 95,
  length = 95,
  output_directory = output_dir,
  writecsv = FALSE,
  pointfinder_data = pointfinder_path,
  pMLST_data = pMLST_data
)

####-----PROCESS METADATA-----####
# Select ColV data to add to meta
colv <- sepsis_simple_summary_N95L95 %>% select(Name = name, ColV, IncF_RST)

# Join plasmid meta and drop STs
meta <- left_join(meta, colv)

# Categorise ColV
meta <- meta %>% mutate(ColV = case_when(ColV == "0" ~ "No", ColV == "1" ~ "Yes"))

# Join phylogroup information
meta <- left_join(meta, phylogroup, by = c("Name"="name"))

# Split STs into groups/ranks based on ST frequency
# common_ST <- meta %>% group_by(ST) %>% filter(n() >= 30) %>% ungroup %>% select(Name) 
# intermediate_ST <- meta %>% group_by(ST) %>% filter(n() >=5, n() <30) %>% ungroup %>% select(Name) 
# rare_ST <- meta %>% group_by(ST) %>% filter(n() < 5) %>% ungroup %>% select(Name)
# 
# # Designate ST ranks
# common_ST$ST_rank <- "Common"
# intermediate_ST$ST_rank <- "Intermediate"
# rare_ST$ST_rank <- "Rare"
# 
# # Stick them back together
# st_rank <- rbind(common_ST, intermediate_ST) %>% select(Name, ST_rank)
# st_rank <- rbind(st_rank, rare_ST) %>% select(Name, ST_rank)
# 
# # Join ST ranks to meta
# meta <- left_join(meta, st_rank)

top_20_ST <- meta %>% 
  group_by(ST) %>% 
  tally(sort=TRUE) %>% 
  slice_head(n=20) %>% 
  pull(ST)

# Create a simple ST column where non-common STs are grouped as 'other' 
meta <- meta %>% mutate(ST_simple = ifelse(ST %notin% top_20_ST, "Other", as.character(ST)))

# Fix the format of the F RST data
meta <- meta %>% mutate(IncF_RST = paste(IncF_RST, str_match(meta$IncF_RST, "C[0-9]{1}"), sep = ":"),
                        IncF_RST = gsub("C[0-9]{1}:", "", IncF_RST),
                        IncF_RST = gsub("^", "F-:", IncF_RST), 
                        IncF_RST = gsub("F-:F", "F", IncF_RST),
                        IncF_RST = gsub(":NA", "", IncF_RST))

# Get the top 15 F types
top15_F <- meta %>% 
  group_by(IncF_RST) %>% 
  tally %>% 
  arrange(desc(n)) %>% 
  slice(1:15) %>% 
  pull(IncF_RST)

# Create column where F types not in the top 15 are designated "Other"
meta <- meta %>% mutate(IncF_simple = ifelse(IncF_RST %notin% top15_F, "Other", as.character(IncF_RST)))

# Extract working names for collection
working_names <- as.vector(meta$Name)

# Get gene data for senB categorisation
senB_cjrABC <- sepsis.N95.L95.PASS %>% 
  filter(grepl("senB", GENE) | grepl("cjrA", GENE) | grepl("cjrB", GENE) |grepl("cjrC", GENE)) %>%
  select(name, SEQUENCE, GENE, perc_coverage, perc_identity) %>% 
  pivot_wider(id_cols = c(name, SEQUENCE), names_from = GENE, values_from = perc_identity) %>% 
  drop_na %>%
  select(-vfdb_senB)

# Assign senB status
senB_cjrABC <- senB_cjrABC %>% mutate(senB_cjrABC = rep("Yes", nrow(senB_cjrABC))) %>% select(name, senB_cjrABC)

# Join to meta
meta <- left_join(meta, senB_cjrABC, by= c("Name"="name")) %>% mutate(senB_cjrABC = replace_na(senB_cjrABC, "No"))

# Categorise sequences as ColV, senB+, or senB- F plasmid carriers (or null)
meta <- meta %>% mutate(`Plasmid_Markers` = case_when(ColV == "No" & IncF_RST != "F-:A-:B-" & senB_cjrABC == "Yes" ~ "senB",
                                                      ColV == "No" & IncF_RST != "F-:A-:B-" & senB_cjrABC == "No" ~ "No Plasmid Markers",
                                                      ColV == "Yes" & IncF_RST != "F-:A-:B-" ~ "ColV",
                                                      IncF_RST == "F-:A-:B-" ~ "No F Plasmid"))

# Re-order F group factor
meta <- meta %>% mutate(`Plasmid_Markers` = as.factor(`Plasmid_Markers`),
                        `Plasmid_Markers` = fct_relevel(`Plasmid_Markers`, "ColV", "senB", "No Plasmid Markers", "No F Plasmid"))

# Vector of ST order
ST_fct_relevel <- meta %>% group_by(ST) %>% tally(sort=TRUE) %>% pull(ST)

# Order ST factor by ST frequency
meta <- meta %>% mutate(ST =as.factor(ST),
                        ST = fct_relevel(ST, ST_fct_relevel))

# Meta table with added count so we can do grouped proportional plots
plot_meta <- meta %>% group_by_all() %>% tally(name = "Count")

# Define levels of ST again for some reason idk
# ST_levels <- plot_meta %>% filter(ST_rank == "Common") %>% group_by(ST) %>% tally(sort=TRUE) %>% pull(ST)

# IncF_RST levels
incf_levels <- meta %>% group_by(IncF_simple) %>% tally(sort=TRUE) %>% pull(IncF_simple)
incf_levels <- c(incf_levels[-(1:2)], incf_levels[1:2])

meta <- meta %>% mutate(IncF_simple = factor(IncF_simple, levels = incf_levels))

# Make count labels for supp plot - not using this anymore?
# rst_lab_inter <- plot_meta %>% filter(ST_rank == "Intermediate") %>% group_by(ST) %>% tally %>% mutate(IncF_simple=NA)

####-----PROCESS pUTI89 ALIGNMENT-----####
# Get tree path and abricate alignment data
path_to_abricate <- "data/full_summaries/pUTI89.tab"
plasrefname <- "pUTI89"

# Get reference length in a very convoluted way...
pUTI89_ref_length <- read_delim(path_to_abricate, n_max = 1) %>% 
  mutate(COVERAGE = as.integer(gsub("^.*\\/", "", COVERAGE))) %>% pull(COVERAGE)

# Use plasmid_mapR to determine binned plasmid coverage at 90% nucleotide ID
plasmid_mapR(path_to_abricate = path_to_abricate,
             plasmid_reference_name = plasrefname,
             output_directory = NULL,
             min_hit_id = 90,
             min_hit_length = 0.5,
             writecsv = FALSE)

# Generate dataframe with a summary statistic for how many bins were covered at >90%ID
pUTI89_ID <- as.data.frame(rowSums(pUTI89_binned_hits))

# Add names column
pUTI89_ID$Name <- rownames(pUTI89_binned_hits)

# Rename columns
colnames(pUTI89_ID) <- c("pUTI89_ID","Name")

# Convert value to a percentage
pUTI89_ID$pUTI89_ID <- round((pUTI89_ID$pUTI89_ID/pUTI89_ref_length) * 100)

# Add pUTI89 IDs to meta
meta <- left_join(meta, pUTI89_ID) %>% mutate(pUTI89_ID = as.numeric(replace_na(pUTI89_ID, 0)))

#
### pCERC4
# Get tree path and abricate alignment data
path_to_abricate <- "data/full_summaries/pCERC4.tab"
plasrefname <- "pCERC4"

# Get ref length
pCERC4_ref_length <- read_delim(path_to_abricate, n_max = 1) %>% 
  mutate(COVERAGE = as.integer(gsub("^.*\\/", "", COVERAGE))) %>% pull(COVERAGE)

# Use plasmid_mapR to determine binned plasmid coverage at 90% nucleotide ID
plasmid_mapR(path_to_abricate = path_to_abricate,
             plasmid_reference_name = plasrefname,
             output_directory = NULL,
             min_hit_id = 90,
             min_hit_length = 0.5,
             writecsv = FALSE)

# Generate dataframe with a summary statistic for how many bins were covered at >90%ID
pCERC4_ID <- as.data.frame(rowSums(pCERC4_binned_hits))

# Add names column
pCERC4_ID$Name <- rownames(pCERC4_binned_hits)

# Rename columns
colnames(pCERC4_ID) <- c("pCERC4_ID","Name")

# Convert value to a percentage
pCERC4_ID$pCERC4_ID <- round((pCERC4_ID$pCERC4_ID/pCERC4_ref_length) * 100)

# Add pCERC4 IDs to meta
meta <- left_join(meta, pCERC4_ID) %>% mutate(pCERC4_ID = as.numeric(replace_na(pCERC4_ID, 0)))


# Plot pUTI89 ID for sequences designated senB+
puti89 <- meta %>% group_by(`Plasmid_Markers`)
ggplot(puti89 %>% filter(`Plasmid_Markers`=="senB") %>% 
         group_by(pUTI89_ID) %>% 
         summarise(count = n()), aes(x=count, y = pUTI89_ID)) +
  geom_col()

# 87.2% of senB+ sequences align to >=70% of pUTI89 at >=90% nucleotide ID

# Plot pCERC4 ID for sequences designated ColV+
pcerc4 <- meta %>% group_by(`Plasmid_Markers`)
ggplot(pcerc4 %>% filter(`Plasmid_Markers`=="ColV") %>% 
         group_by(pCERC4_ID) %>% 
         summarise(count = n()), aes(x=pCERC4_ID, y = count)) +
  geom_col()



####-----MGE-CLUSTER METADATA-----####
F_type_colours <- c("ColV" = "#17e6ae", "senB" = "#E7298A", "No Plasmid Markers" = "#dbdbdb", "No F Plasmid" = "#f0efe6")

# Common STs
ST_vars <- unique(meta$ST_simple)
ST_clrs <- colorRampPalette(brewer.pal(12, "Paired"))(length(unique(meta$ST_simple)))
names(ST_clrs) <- sort(ST_vars)
ST_clrs["Other"] <- "#dbdbdb"

mge_data <- read_delim("data/mge-cluster/results/F_only/mge-cluster_prediction_minclust-80.csv", delim = ",", col_names = TRUE, trim_ws = TRUE) %>%
  mutate(Sample_Name = gsub("_plasmidcontigs", "", Sample_Name)) %>% 
  mutate(Standard_Cluster = ifelse(Standard_Cluster == "-1", "pc_null", as.character(Standard_Cluster)),
         Standard_Cluster = ifelse(Standard_Cluster == "-", "pc_null", as.character(Standard_Cluster)),
         Standard_Cluster = ifelse(Standard_Cluster != "pc_null", paste0("pc_", Standard_Cluster), as.character(Standard_Cluster)))

# Get the names of clusters that represent >5% of clustered genomes
clust_gtr_100 <- mge_data %>% 
  group_by(Standard_Cluster) %>% 
  tally %>% 
  filter(n>0.05*(nrow(mge_data %>% filter(Standard_Cluster != "pc_null")))) %>% 
  pull(Standard_Cluster)

mge_data <- mge_data %>% 
  mutate(Standard_Cluster_simple = ifelse(Standard_Cluster %notin% clust_gtr_100, "Other", as.character(Standard_Cluster)))

mge_meta <- left_join(meta, mge_data, by= c("Name"= "Sample_Name")) %>%
  select(Name, ST, ST_simple, tsne1D, tsne2D, Standard_Cluster, Standard_Cluster_simple, IncF_simple, IncF_RST, `Plasmid_Markers`)

mge_vars <- mge_meta %>% filter(!is.na(Standard_Cluster)) %>% pull(Standard_Cluster) %>% unique()
mge_clrs <- colorRampPalette(brewer.pal(12, "Set3"))(length(mge_vars))
names(mge_clrs) <- sort(mge_vars)
mge_clrs["pc_null"] <- "#b5b3b3"
mge_clrs["None"] <- "#d9d9d9"

mge_plot_meta <- mge_meta %>% filter(!is.na(tsne1D),
                                     tsne1D != "-") %>% 
  mutate(tsne1D = as.numeric(tsne1D), tsne2D = as.numeric(tsne2D))

simp_mge_clrs <- mge_clrs[names(mge_clrs) %in% clust_gtr_100]
simp_mge_clrs["Other"] <- "#d9d9d9"

simp_mge_clrs2 <- mge_clrs[names(mge_clrs) %in% clust_gtr_100]
simp_mge_clrs2["Other"] <- "#d9d9d9"
simp_mge_clrs2["None"] <- "#faf9f0"

# Add to meta
meta <- left_join(meta, mge_meta %>% select(Name, Standard_Cluster, Standard_Cluster_simple)) %>% 
  mutate(Standard_Cluster_simple = replace_na(Standard_Cluster_simple, "None")) %>% 
  mutate(Standard_Cluster_simple = factor(Standard_Cluster_simple, levels = names(simp_mge_clrs2)))

####-----INTER-TYPE RSTS-----####
# have a look at RSTs that are present in different Plasmid_Markerss
type_incf <- meta %>% group_by(Name, `Plasmid_Markers`, IncF_RST) %>% summarise() %>% pivot_wider(names_from = `Plasmid_Markers`, values_from = IncF_RST)

colv_incs <- type_incf$ColV %>% na.omit() %>% unique()
senb_incs <- type_incf$`senB` %>% na.omit() %>% unique()
other_incs <- type_incf$`Other` %>% na.omit() %>% unique()

a <- intersect(colv_incs, senb_incs)

b <- intersect(senb_incs, other_incs)

c <- intersect(colv_incs, other_incs)

y <- intersect(b,c)
intersect(y,a)

x <- union(a,c)
union(x, b)

unique(meta$IncF_RST) %>% length()

####------COLOURS-----####
# Common STs
ST_vars <- unique(meta$ST_simple)
ST_clrs <- colorRampPalette(brewer.pal(12, "Paired"))(length(unique(meta$ST_simple)))
names(ST_clrs) <- sort(ST_vars)
ST_clrs["Other"] <- "#dbdbdb"

# ST Ranks
rank_vars <- unique(meta$ST_rank)
rank_clrs <- brewer.pal(3, "Spectral")
names(rank_clrs) <- rank_vars

# Phylogroup colours
phylogroup_vars <- unique(meta$phylogroup)[1:11]
phylogroup_clrs <- colorRampPalette(brewer.pal(11, "Set3"))(length(unique(meta$phylogroup)[1:11]))
names(phylogroup_clrs) <- sort(phylogroup_vars)

# F RST Colours
F_vars <- unique(meta$IncF_simple)
F_clrs <- colorRampPalette(brewer.pal(12, "Paired"))(length(unique(meta$IncF_simple)))
names(F_clrs) <- F_vars 
F_clrs["Other"] <- "#dbdbdb"
F_clrs["F-:A-:B-"] <- "#FFFFFF"

# Alternative F colours
alt_F_clrs <- c("#AA4499", "#795e86", "#9b7121", "#4c94c0", "#DDCC77", "#f2c8b3", "#CC6677","#44AA99", 
                "#e91c04", "#523cbd", "#9a8106", "#117733", "#88CCEE", "#73e3a7", "#882255", "#9e8e96")

set.seed(12)
alt_F_clrs <- sample(alt_F_clrs, length(unique(meta$IncF_simple)))
ordered_F <- meta %>% group_by(IncF_simple) %>% tally(sort=TRUE) %>% pull(IncF_simple)
ordered_F <- c(ordered_F[-1:-2], ordered_F[2:1])
names(alt_F_clrs) <- ordered_F
alt_F_clrs["F-:A-:B-"] <- "#f0efe6"
alt_F_clrs["Other"] <- "#d3d2d9"


# see_palette(alt_F_clrs)

# Colours for categorised F types
#F_type_colours <- c("ColV" = "#17e6ae", "senB" = "#26e9ff", "Other" = "#dec3fa", "None" = "#dbdbdb")
F_type_colours <- c("ColV" = "#17e6ae", "senB" = "#E7298A", "No Plasmid Markers" = "#26e9ff", "No F Plasmid" = "#f0efe6")
F_type_colours <- c("ColV" = "#17e6ae", "senB" = "#E7298A", "No Plasmid Markers" = "grey", "No F Plasmid" = "#f0efe6")

# ColV
colv_clrs <- c("ColV-Pos" = "#17e6ae", "ColV-Neg" = 'white')

colv_clrs_meta <- c("Yes" = "#17e6ae", "No" = '#dbdbdb')

# pUTI89
puti89_clrs <- c("pUTI89-Pos" = "#E7298A", "pUTI89-Neg" = 'white')

# senB
senB_clrs <- c("senB" = "#26e9ff", "senB-Neg" = 'white')

# Combine for tree
tree_vars <- c(ST_clrs, F_type_colours, F_clrs)

# Select phylogroup data for tips
phylo_tips <- meta %>% select(Name, phylogroup)

####-----MOBSUITE DATA-----####
pmob <- read_csv('data/mobsuite/plasmid_concatenated_contig_report.csv')

pmob <- pmob %>% 
  filter(molecule_type == "plasmid") %>% 
  select("Name", "primary_cluster_id", "secondary_cluster_id", "rep_type" = "rep_type(s)", "relaxase_type" = "relaxase_type(s)", "mpf_type", "orit_type" = "orit_type(s)")

pmob_filt <- pmob %>%
  group_by(Name, primary_cluster_id) %>%
  summarize_all(~paste(unique(na.omit(.)), collapse = ',')) %>% 
  ungroup() %>% 
  filter(grepl("IncF", rep_type), grepl("MOBF", relaxase_type), grepl("MOBF", orit_type))

pmob_filt_no_orit <-pmob %>%
  group_by(Name, primary_cluster_id) %>%
  summarize_all(~paste(unique(na.omit(.)), collapse = ',')) %>% 
  ungroup() %>% 
  filter(grepl("IncF", rep_type), grepl("MOBF", relaxase_type))

pmob_meta <- left_join(meta %>% select(Name, IncF_RST), pmob_filt)

# Get the top 15 primary clusters
top15_pclust <- pmob_filt %>% 
  group_by(primary_cluster_id) %>% 
  tally %>% 
  arrange(desc(n)) %>% 
  slice(1:15) %>% 
  pull(primary_cluster_id)

# Top 15 secondary clusters
top15_sclust <- pmob_filt %>% 
  group_by(secondary_cluster_id) %>% 
  tally %>% 
  arrange(desc(n)) %>% 
  slice(1:15) %>% 
  pull(secondary_cluster_id)

pmob_filt <- pmob_filt %>% 
  mutate(primary_cluster_id = ifelse(is.na(primary_cluster_id), "None", as.character(primary_cluster_id)),
         secondary_cluster_id = ifelse(is.na(secondary_cluster_id), "None", as.character(secondary_cluster_id)),
         primary_cluster_simple = ifelse(primary_cluster_id %in% top15_pclust, as.character(primary_cluster_id), "Other"),
         secondary_cluster_simple = ifelse(secondary_cluster_id %in% top15_sclust, as.character(secondary_cluster_id), "Other"))

####------GENE PROCESSING-----####
# Create 'geno_meta' df with all meta and gene screening data
geno_meta <- left_join(meta %>% select(Name, ST, ColV, IncF_RST, IncF_simple, phylogroup, ST_simple, `Plasmid_Markers`, Standard_Cluster, Standard_Cluster_simple), 
                       sepsis_simple_summary_N95L95 %>% select(-ColV, -starts_with("Inc")), by = c("Name" = "name"))

# Define gene columns as those that are integers
gene_cols <- names(geno_meta %>% select(where(is.integer)))

# Recode multiple hits as a single hit
geno_meta <- geno_meta %>% mutate(across(where(is.integer), ~ gsub(2, 1, .x)))
geno_meta <- geno_meta %>% mutate(across(where(is.integer), ~ gsub(3, 1, .x)))
geno_meta <- geno_meta %>% mutate(across(where(is.integer), ~ gsub(4, 1, .x)))

# Convert back to integer (gsub makes everything character class)
geno_meta <- geno_meta %>% mutate(across(all_of(gene_cols), as.integer))

# Extract metadata column names for later use with the geno_meta dataframe
meta_cols <- c("Name", "ST", "ColV", "IncF_RST", "IncF_simple", "phylogroup", "ST_simple", "Plasmid_Markers", "Standard_Cluster", "Standard_Cluster_simple")

# senB data
senB <- geno_meta %>% select(Name, senB = vfdb_senB) %>% mutate(senB = case_when(senB == "1" ~ "senB", senB == "0" ~ "senB-Neg"))

# Split ABRicate gene hits into their functional groups
# Get all the hits from CARD database and intI1 and intI2. Fix all the messy names. Filter out genes present in >90% of isolates. These are housekeeping
# genes that sometimes mutate to confer AMR phenotypes but we are only concerned with acquired resistance genes
args <- geno_meta %>%
  select(where(is.character), where(is.factor), where( ~ is.integer(.x) && sum(.x) >0 && sum(.x) <= .95*nrow(geno_meta))) %>%
  select(all_of(meta_cols), starts_with("card"), contains("intI")) %>%
  rename_with(~ gsub("card_", "", .x, fixed = TRUE)) %>%
  rename_with(~ gsub("Escherichia_coli_", "", .x, fixed = TRUE)) %>%
  rename_with(~ gsub("_beta-lactamase", "", .x, fixed = TRUE)) %>%
  rename_with(~ gsub("(", "_", .x, fixed = TRUE)) %>%
  rename_with(~ gsub(")", "", .x, fixed = TRUE)) %>%
  rename_with(~ gsub("PC1__", "PC1_", .x, fixed = TRUE)) %>%
  rename_with(~ gsub("EC_custom_intI1.*", "intI1", .x)) %>%
  rename_with(~ gsub("EC_custom_intI2.*", "intI2", .x)) %>%
  rename_with(~ gsub("Shigella_flexneri_chloramphenicol_acetyltransferase", "catA1", .x, fixed = TRUE)) %>% 
  rename_with(~ gsub("Klebsiella_pneumoniae_", "", .x, fixed = TRUE)) %>%
  select(sort(names(.)), -ugd) %>%
  relocate(all_of(meta_cols), contains("intI"))

# Calculate total carriage of each gene in the collection
arg_totals <- t(args %>% summarise(across(where(is.integer), sum)))

arg_totals <- as_tibble(arg_totals, rownames = "Gene", .name_repair = "minimal") %>% 
  rename(Total = 2) %>% 
  mutate(Percentage = round(Total/nrow(geno_meta)*100, 2))

# Total ARGs
args2 <- args %>% select(-intI1, -intI2) %>% 
  mutate(`Total ARGs` = rowSums(across(where(is.integer))), 
         `Plasmid_Markers` = factor(`Plasmid_Markers`, levels = c("ColV", "senB", "Other", "None"))) %>% 
  select(all_of(meta_cols), `Total ARGs`)

# Get hits from VFDB
vags <- geno_meta %>% 
  select(all_of(meta_cols), starts_with("vfdb")) %>%
  rename_with(~ gsub("vfdb_", "", .x))

# Select additional virulence genes from our custom database
custom_vags <- geno_meta %>% 
  select(Name, contains("EC_custom")) %>%
  rename_with(~ gsub("EC_custom_", "", .x, fixed = TRUE)) %>%
  mutate(usp = as.integer(rowSums(select(.,starts_with("usp"))))) %>%
  select(-starts_with("usp_")) %>%
  mutate(eitA = as.integer(rowSums(select(.,starts_with("eitA"))))) %>%
  select(-starts_with("eitA_")) %>% 
  select(Name,
         starts_with(c("cba", "cbi", "cjr", 
                       "cva", "cvi","eit",
                       "fecA", "hek", "hyl",
                       "iha","iss", "merA",
                       "ompT", "silA",
                       "terA", "traT", "usp"))) %>%
  rename_with(~ gsub("_[A-Z]{1,2}.*", "", .x)) %>%
  rename_with(~ gsub("_pAPEC-O1-ColBM", "", .x, fixed =TRUE)) %>%
  rename_with(~ gsub("_pUTI89", "", .x, fixed =TRUE)) %>%
  rename_with(~ gsub("_type3", "", .x, fixed =TRUE)) %>%
  rename_with(~ gsub("_|VFG1539", "", .x, fixed =TRUE)) %>%
  rename_with(~ gsub("_chromosomal", "_1", .x, fixed =TRUE)) %>%
  rename_with(~ gsub("_episomal", "_2", .x, fixed =TRUE))

# Join VFDB and custom hits and filter genes in less than 5% of strains
vags <- left_join(vags, custom_vags) %>% select(sort(names(.))) %>%
  relocate(all_of(meta_cols)) %>% 
  select(where(is.character), where(is.factor), where( ~ is.integer(.x) && sum(.x) >= .05*nrow(vags)))

# Total VAGs
vags2 <- vags %>% 
  mutate(`Total VAGs` = rowSums(across(where(is.integer))),
         `Plasmid_Markers` = factor(`Plasmid_Markers`, levels = c("ColV", "senB", "Other", "None"))) %>% 
  select(all_of(meta_cols), `Total VAGs`)

# Infer HPI carriage by presence of fyuA and irp2
HPI <- vags %>% mutate(HPI = case_when(fyuA == "1" & irp2 == "1" ~ "1", TRUE ~ "0")) %>% select(Name, HPI)

# Add HPI to meta
meta <- left_join(meta, HPI)

#####################################
####---- FIGURE 1 MGE CLUSTER ----####
#####################################
#### FIG 1A - RESULTS MINCLUST-80 ####
f <- read_delim("data/mge-cluster/results/F_only/mge-cluster_prediction_minclust-80.csv", delim = ",", col_names = TRUE, trim_ws = TRUE) %>%
  mutate(Sample_Name = gsub("_plasmidcontigs", "", Sample_Name)) %>% 
  mutate(Standard_Cluster = ifelse(Standard_Cluster == "-1", "pc_null", as.character(Standard_Cluster)),
         Standard_Cluster = ifelse(Standard_Cluster == "-", "pc_null", as.character(Standard_Cluster)),
         Standard_Cluster = ifelse(Standard_Cluster != "pc_null", paste0("pc_", Standard_Cluster), as.character(Standard_Cluster)))

f_meta <- mge_meta %>% 
  filter(!is.na(tsne1D),
         tsne1D != "-") %>% 
  mutate(tsne1D = as.numeric(tsne1D), tsne2D = as.numeric(tsne2D))

# Define colours for specific clusters 
# mge-clust colours
mge_vars2 <- unique(f_meta$Standard_Cluster)
mge_clrs2 <- colorRampPalette(brewer.pal(12, "Set3"))(length(unique(f_meta$Standard_Cluster)))
names(mge_clrs2) <- sort(mge_vars2)
mge_clrs2["pc_null"] <- "#dbdbdb"

# Make plots and assign their names
## mge-clusters
f_meta_clusters <- f_meta %>% filter(Standard_Cluster != "pc_null")
f_meta_other <- f_meta %>% filter(Standard_Cluster == "pc_null")

##
## pc_null is causing the whole thing to go weird when you coerce tsne columns to pc_null
mgeplot <- ggplot() +
  geom_point(data = f_meta_other,
             aes(x = tsne1D, 
                 y = tsne2D,
                 color = Standard_Cluster), size = .3) +
  scale_color_manual(values = mge_clrs2,
                     name = "F Plasmid Cluster")+
  #xlim(-40, 45) +
  #ylim(-40, 30) +
  geom_point(data = f_meta_clusters,
             aes(x = tsne1D, 
                 y = tsne2D,
                 color = Standard_Cluster), size =.3) +
  scale_color_manual(values = mge_clrs2, 
                     name = "F Plasmid Cluster") +
  geom_mark_ellipse(data = f_meta_clusters, aes(x = tsne1D, 
                                                y = tsne2D,
                                                color = Standard_Cluster), 
                    expand = unit(0.5,"mm"),
                    show.legend = FALSE)+
  #xlim(-20,45) +
  #ylim(-45, 5) +
  theme_classic() +
  theme(axis.title = element_text(size = 8),
        axis.text = element_text(size=6),
        axis.line = element_blank(),
        legend.position = c(0.1,.18),
        legend.justification = "center",
        legend.key.size = unit(5, "mm"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6, hjust = 0.5),
        legend.direction = "vertical",
        plot.margin = margin(.2,2,.2,2, "cm"),
        panel.border = element_rect(color = "black", 
                                    fill = NA, 
                                    linewidth = .5),
        legend.spacing.x = unit(0, 'cm')) +
  guides(color = guide_legend(ncol=2, byrow=TRUE, override.aes = list(size=1)))

#### FIG 1B F GROUP CLUSTER ####
F_plot <- ggplot(f_meta,
                 aes(x=tsne1D, 
                     y =tsne2D)) + 
  geom_point(aes(color = `Plasmid_Markers`), size =.1) +
  geom_mark_ellipse(data = f_meta_clusters, aes(x = tsne1D, 
                                                y = tsne2D,
                                                color = Standard_Cluster), 
                    expand = unit(0.5,"mm"),
                    show.legend = FALSE)+
  scale_color_manual(values = c(mge_clrs2, F_type_colours), breaks= c(names(F_type_colours)),
                     name = "Plasmid Markers") +
  #xlim(-20,45) +
  #ylim(-45, 5) +
  theme_classic() +
  theme(axis.text = element_text(size=6),
        axis.title = element_text(size = 8),
        axis.line = element_blank(),
        legend.position = "none",
        legend.justification = "center",
        legend.key.size = unit(2.5, "mm"),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 5, hjust = 0.5),
        legend.direction = "vertical",
        plot.margin = margin(.5,.1,.7,.1, "cm"),
        panel.border = element_rect(color = "black", 
                                    fill = NA, 
                                    linewidth = .5)) +
  guides(color = guide_legend(nrow=1, byrow=TRUE, override.aes = list(size=1)))

#### FIG 1C RST CLUSTER ####
rst_plot <- ggplot(f_meta,
                   aes(x=tsne1D, 
                       y =tsne2D)) + 
  geom_point(aes(color = IncF_simple), size =.1) +
  geom_mark_ellipse(data = f_meta_clusters, aes(x = tsne1D, 
                                                y = tsne2D,
                                                color = Standard_Cluster), 
                    expand = unit(0.5,"mm"),
                    show.legend = FALSE)+
  scale_color_manual(values = c(mge_clrs2, alt_F_clrs), breaks= c(names(F_clrs)),
                     name = "IncF RST") +
  #xlim(-20,45) +
  #ylim(-45, 5) +
  theme_classic() +
  theme(axis.title = element_text(size = 8),
        axis.text.x = element_text(size=6),
        axis.text.y = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        legend.justification = "center",
        legend.key.size = unit(2.5, "mm"),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 5, hjust = 0.5),
        legend.direction = "vertical",
        plot.margin = margin(.5,.1,.7,.1, "cm"),
        panel.border = element_rect(color = "black", 
                                    fill = NA, 
                                    linewidth = .5)) +
  guides(color = guide_legend(nrow=3, byrow=TRUE, override.aes = list(size=1)))

#### FIG 1D ST CLUSTER ####
ST_plot <- ggplot(f_meta,
                  aes(x=tsne1D, 
                      y =tsne2D)) + 
  geom_point(aes(color = ST_simple), size =.1) +
  geom_mark_ellipse(data = f_meta_clusters, aes(x = tsne1D, 
                                                y = tsne2D,
                                                color = Standard_Cluster), 
                    expand = unit(0.5,"mm"),
                    show.legend = FALSE)+
  scale_color_manual(values = c(ST_clrs, mge_clrs2), breaks = c(names(ST_clrs)),
                     name = "ST") +
  #xlim(-20,45) +
  #ylim(-45, 5) +
  theme_classic() +
  theme(axis.text.x = element_text(size=6),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 8),
        axis.title.y = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none",
        legend.justification = "center",
        legend.key.size = unit(2.5, "mm"),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 5, hjust = 0.5),
        legend.direction = "vertical",
        plot.margin = margin(.5,.1,.7,.1, "cm"),
        panel.border = element_rect(color = "black", 
                                    fill = NA, 
                                    linewidth = .5)) +
  guides(color = guide_legend(nrow=3, byrow=TRUE, override.aes = list(size=1)))

# Save plot as PDF with title as filename
#ggsave("outputs/figures/mge-cluster100.pdf", plot = final_plot)



plot_meta <- meta %>% group_by_all() %>% tally(name = "Count")

##### FIG 1E F GROUP COL ####
# F Plasmid category by F Group with labels
FSCs <- meta %>% filter(Standard_Cluster_simple != "None") %>%
  group_by(Standard_Cluster_simple, `Plasmid_Markers`) %>% tally %>% mutate(perc=n/sum(n))

fig1e <- ggplot(FSCs , aes(x = Standard_Cluster_simple, y=perc, fill=`Plasmid_Markers`)) +
  geom_col(position = "stack") +
  #geom_text(aes(label = ifelse(perc > .02, IncF_simple, "")), position = position_stack(vjust=0.5), size =1.8) +
  scale_fill_manual(name = "Plasmid Markers", values = F_type_colours,
                    guide = guide_legend(nrow=1)) +
  #scale_x_discrete(name = "Cluster")+
  scale_y_continuous(name = "Proportion", expand = c(0, 0), n.breaks = 5) +
  theme_classic()+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(color = "grey20", 
                                   size = 6, 
                                   vjust = 1,
                                   hjust = ,
                                   face = "plain",
                                   angle = ),
        axis.text.y = element_text(color = "grey20", 
                                   size = 6, 
                                   face = "plain"),
        legend.position = "bottom",
        legend.key.size = unit(2, "mm"),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.direction = "vertical")

#### LABEL COLOURS E-G ####
label_clrs <- c("pc_2" = "#D9A0AB", "pc_4" = "#EFB36E", "pc_5" = "#BBDC76", "pc_6" = "#F0D1E1", "pc_7" = "#C8A7C9")


##### FIG 1F RST COL ####
# Mutate data to make labels possible
SCs <- meta %>% filter(Standard_Cluster_simple != "None") %>%
  group_by(Standard_Cluster_simple, IncF_simple) %>% tally %>% mutate(perc=n/sum(n))

# F Plasmid category by F RST with labels
fig1f <- ggplot(SCs , aes(x = Standard_Cluster_simple, y=perc, fill=IncF_simple)) +
  geom_col(position = "stack") +
  geom_text(aes(label = ifelse(perc > .05, as.character(IncF_simple), "")), position = position_stack(vjust=0.5), size =1) +
  scale_fill_manual(name = "F Plasmid RST", values = alt_F_clrs,
                    guide = guide_legend(nrow=4)) +
  #scale_x_discrete(name = "Cluster")+
  scale_y_continuous(name = "Proportion", expand = c(0, 0), n.breaks = 5) +
  theme_classic()+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(color = "grey20", 
                                   size = 6, 
                                   vjust = 1,
                                   hjust = ,
                                   face = "plain"),
        legend.position = "bottom",
        legend.key.size = unit(2, "mm"),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.direction = "vertical")


##### FIG 1G ST COL ####
# F Plasmid category by F RST with labels
FSTs <- meta %>% filter(Standard_Cluster_simple != "None") %>%
  mutate(ST_simple = factor(ST_simple, levels = c(top_20_ST, "Other"))) %>% 
  group_by(Standard_Cluster_simple, ST_simple) %>% tally %>% mutate(perc=n/sum(n))

fig1g <- ggplot(FSTs , aes(x = Standard_Cluster_simple, y=perc, fill=ST_simple)) +
  geom_col(position = "stack") +
  geom_text(aes(label = ifelse(perc > .05, as.character(ST_simple), "")), position = position_stack(vjust=0.5), size =1.8) +
  scale_fill_manual(name = "ST", values = ST_clrs,
                    guide = guide_legend(nrow=4)) +
 #scale_x_discrete(name = "Cluster")+
  scale_y_continuous(name = "Proportion", expand = c(0, 0), n.breaks = 5) +
  theme_classic()+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(color = "grey20", 
                                   size = 6, 
                                   vjust = 1,
                                   hjust = ,
                                   face = "plain",
                                   angle = ),
        legend.position = "bottom",
        legend.key.size = unit(2, "mm"),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.direction = "vertical")

##### SAVE FIG 1 ####
fig1 <- ggarrange(mgeplot,
                        ggarrange(F_plot,rst_plot,ST_plot,
                                  labels =c("b)","c)", "d)"),
                                  ncol =3,
                                  font.label = list(size=8),
                                  align = "v"),
                        ggarrange(fig1e,fig1f,fig1g,
                                  labels =c("e)","f)", "g)"),
                                  ncol = 3,
                                  font.label = list(size=8),
                                  align = "hv"),
                        ncol =1, nrow =3, 
                        labels =c("a)"),
                        common.legend = FALSE,
                        heights = c(1.5,1,1),
                        font.label = list(size=8)) 


ggsave("outputs/figures/fig1.clusters.pdf", 
       plot = fig1,
       width = 210,
       height =297,
       units="mm")

#####################################
####----- FIGURE 2 MAJOR STS ----####
#####################################
#### FIGURE 2 ####
# f2_meta <- meta %>% group_by(ST) %>% filter(n() >= 65, Standard_Cluster_simple != "None") %>% ungroup() %>% 
#  mutate(ST = factor(paste0("ST",ST), levels = paste0("ST", levels(fct_infreq(meta$ST)))))
f2_meta <- meta %>%
  filter(ST %in% (meta %>%
                    count(ST, sort = TRUE) %>%
                    slice_head(n = 20) %>%
                    pull(ST))) %>%
  mutate(ST = factor(paste0("ST", ST), levels = paste0("ST", levels(
    fct_infreq(meta$ST)
  ))))

alt_F_clrs2 <- alt_F_clrs[-(length(alt_F_clrs) -1)]

fig2 <- ggplot(f2_meta, aes(x = Standard_Cluster_simple)) +
  geom_bar(aes(fill = `IncF_simple`)) +
  facet_wrap(
    ~ ST,
    nrow = 4,
    strip.position = "top",
    scales = "free_y"
  ) +
  scale_fill_manual(values = alt_F_clrs,
                    name = "F Plasmid RST",
                    guide = guide_legend(ncol = 1)) +
  scale_x_discrete(name = "Plasmid Cluster") +
  scale_y_continuous(name = "Count",
                     expand = expansion(mult = c(0, 0.02)),
                    n.breaks = 6,
                     labels = scales::number_format(accuracy = 1, trim = TRUE)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      color = "grey20",
      size = 9,
      vjust = 1,
      hjust = 1,
      face = "plain",
      angle = 40
    ),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size=14),
    axis.line = element_blank(),
    legend.position = "bottom",
    legend.key.size = unit(5.5, "mm"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.direction = "horizontal",
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.x = element_text(size = 12),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = .7
    )
  ) +
  guides(fill = guide_legend(nrow = 2, byrow=TRUE, override.aes = list(size=1)))

#### SAVE FIGURE 2 ####
ggsave("fig2.commonSTs.clusters.RSTs.pdf", 
       fig2, 
       path = "outputs/figures", 
       device = "pdf", 
       width = 270,
       height = 180,
       unit ="mm")

# g <- ggplot_gtable(ggplot_build(p))
# stripb <- which(grepl('strip-t', g$layout$name))
# fills <- colorRampPalette(brewer.pal(10, "Set3"))(10)
# k <- 1
# for (i in stripb) {
#   j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
#   g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
#   k <- k+1
# }
# grid::grid.draw(g)

#####################################
#### FIGURE S1 MGE-CLUSTER MODEL ####
#####################################
#### FIGURE S1 ####
#Get input filenames
model <- read_delim("data/mge-cluster/models/F-plasmid-model-minclust-80_results.csv", delim = ",", col_names = TRUE, trim_ws = TRUE) %>% 
  mutate(Plasmid_Cluster = ifelse(Standard_Cluster == "-1", "pc_null", as.character(Standard_Cluster)),
         Plasmid_Cluster = ifelse(Plasmid_Cluster != "pc_null", paste0("pc_", Plasmid_Cluster), as.character(Plasmid_Cluster)))

# Define colours for specific clusters 
# mge-clust colours
extra_clrs <- c("pc_1" = "#FB8072", "pc_9" = "#80B1D3")
mge_clrs <- c(mge_clrs2, extra_clrs)
mge_clrs["pc_null"] <- "#cfad95"

model_clusters <- model %>% filter(Plasmid_Cluster != "pc_null")
model_other <- model %>% filter(Plasmid_Cluster == "pc_null")

figS1 <- ggplot() + 
  geom_point(data = model_other,
             aes(x=tsne1D, 
                 y =tsne2D,
                 color = Plasmid_Cluster), 
             size = 1) +
  scale_color_manual(values = mge_clrs, name = "Plasmid Cluster") +
  geom_point(data = model_clusters,
             aes(x=tsne1D, 
             y =tsne2D,
             color = Plasmid_Cluster), 
             size = 1) +
  scale_color_manual(values = mge_clrs, name = "Plasmid Cluster") +
  geom_mark_ellipse(data = model_clusters, 
                    aes(x = tsne1D, 
                    y = tsne2D,
                    color = Plasmid_Cluster), 
                    expand = unit(0.5,"mm"),
                    show.legend = FALSE) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.key.size = unit(10, "mm"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size =10),
        panel.border = element_rect(color = "black", 
                                    fill = NA, 
                                    linewidth = .5)) +
  guides(color = guide_legend(nrow =2)) 

#### SAVE FIGURE S1 ####
ggsave("figS1.mge-model.pdf", 
       figS1, 
       path = "outputs/figures", 
       device = "pdf", 
       width = 297,
       height = 210,
       unit ="mm")

#####################################
####---- FIGURE S2 COMMON STs ----####
#####################################
# Meta table with added count so we can do grouped proportional plots
plot_meta <- meta %>% group_by_all() %>% tally(name = "Count")

# Levels to plot STs in collection frequency order
ST_levels <- meta %>% filter(ST_simple != "Other") %>% group_by(ST_simple) %>% tally(sort=TRUE) %>% pull(ST_simple)

#### FIGURE S2A ####
# Absolute common ST by F category
figS2a <- ggplot(meta %>% filter(ST_simple != "Other"), aes(x = fct_infreq(factor(ST_simple)))) +
  geom_bar(aes(fill = Standard_Cluster_simple), position="stack") +
  scale_fill_manual(values = simp_mge_clrs2,
                    guide = guide_legend(title.position = "left"),
                    name = "Cluster") +
  scale_x_discrete(name = "ST")+
  scale_y_continuous(name = "Count", expand = c(0, 0), limits = c(0,800), n.breaks = 9) +
  theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", 
                                   size = 6, 
                                   vjust = 1,
                                   hjust = ,
                                   face = "plain",
                                   angle = ),
        axis.title.x = element_blank(),
        legend.position = "right",
        legend.key.size = unit(2.5, "mm"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        legend.direction = "horizontal")

#### FIGURE S2B ####
# Make count labels for below plot
# f_labs_common <- plot_meta %>% filter(ST_rank == "Common") %>% group_by(ST) %>% tally %>% mutate(Standard_Cluster_simple=NA)

# Proportional common ST by F category
figS2b <- ggplot(plot_meta %>% filter(ST_simple != "Other") %>% group_by(ST_simple, Standard_Cluster_simple) %>% tally %>% mutate(perc=n/sum(n)), 
                aes(x = factor(ST_simple, levels = ST_levels), y = perc, fill =Standard_Cluster_simple)) +
  geom_col()+
  #geom_text(data=f_labs_common,aes(label=n, y = rep(0.97,length(unique(meta %>% filter(ST_rank == "Common") %>% pull(ST)))))) +
  scale_fill_manual(values = simp_mge_clrs2,
                    guide = guide_legend(title.position = "left")) +
  scale_x_discrete(name = "ST")+
  scale_y_continuous(name = "Proportion", expand = c(0, 0), limits = c(0, 1.01), n.breaks = 6) +
  theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", 
                                   size = 6, 
                                   vjust = 1,
                                   hjust = ,
                                   face = "plain",
                                   angle = ), 
        legend.position = "right",
        legend.key.size = unit(2.5, "mm"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        legend.direction = "horizontal")

#### FIGURE S2C ####
# Absolute common ST by F RST
figS2c <- ggplot(meta %>% filter(ST_simple != "Other"), aes(x = fct_infreq(factor(ST_simple)))) +
  geom_bar(aes(fill = IncF_simple)) +
  scale_fill_manual(values = alt_F_clrs, name = "F Plasmid RST",
                    guide = guide_legend(nrow=2)) +
  scale_x_discrete(name = "ST")+
  scale_y_continuous(name = "Count", expand = c(0, 0), limits=c(0,800), n.breaks = 9) +
  theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", 
                                   size = 6,
                                   vjust = ,
                                   face = "plain"),
        axis.title.x = element_blank(),
        legend.position = "none",
        legend.key.size = unit(2.5, "mm"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        legend.direction = "horizontal",
  )

#### FIGURE S2D ####
# Proportional common ST by F RST
figS2d <- ggplot(plot_meta %>% filter(ST_simple != "Other") %>% group_by(ST, IncF_simple) %>% tally %>% mutate(perc=n/sum(n)), 
                aes(x = factor(ST, levels = ST_levels), y = perc, fill = IncF_simple)) +
  geom_col()+
  scale_fill_manual(values = alt_F_clrs, name = "F Plasmid RST",
                    guide = guide_legend(nrow=2)) +
  scale_x_discrete(name = "ST") +
  scale_y_continuous(name = "Proportion", expand = c(0, 0), n.breaks = 6) +
  theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", 
                                   size = 6,
                                   vjust = ,
                                   face = "plain"),
        legend.position = "none",
        legend.key.size = unit(2.5, "mm"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        legend.direction = "horizontal")

#### COMPILE AND SAVE FIGURE S2 ####
figS2ab <- ggarrange(plotlist= list(figS2a, figS2b),
                    ncol =1, nrow =2, labels =c("a)", "b)"),
                    font.label = list(size= 10),
                    common.legend = TRUE,
                    legend = "bottom")

figS2cd <- ggarrange(plotlist= list(figS2c, figS2d),
                    ncol =1, nrow =2, labels =c("c)", "d)"),
                    font.label = list(size= 10),
                    common.legend = TRUE,
                    legend = "bottom")

figS2 <- ggarrange(plotlist = list(figS2ab, figS2cd),
                  ncol =2,
                  legend = "bottom")

ggsave("figS2.ST.clusters.RSTs.pdf", 
       figS2, 
       path = "outputs/figures", 
       device = "pdf", 
       width = 297,
       height = 210,
       unit ="mm")

#####################################
####---- FIGURE S3 STs MARKERS----####
#####################################
#### FIGURE S3 ####
figS3 <- ggplot(f2_meta, aes(x = Standard_Cluster_simple)) +
  geom_bar(aes(fill = Plasmid_Markers)) +
  facet_wrap(
    ~ ST,
    nrow = 4,
    strip.position = "top",
    scales = "free_y"
  ) +
  scale_fill_manual(values = F_type_colours,
                    name = "Plasmid Markers",
                    guide = guide_legend(ncol = 1)) +
  scale_x_discrete(name = "Plasmid Cluster") +
  scale_y_continuous(name = "Count",
                     expand = expansion(mult = c(0, 0.02)),
                     n.breaks = 6,
                     labels = scales::number_format(accuracy = 1, trim = TRUE)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      color = "grey20",
      size = 9,
      vjust = 1,
      hjust = 1,
      face = "plain",
      angle = 30
    ),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size =14),
    axis.line = element_blank(),
    legend.position = "bottom",
    legend.key.size = unit(5.5, "mm"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.direction = "horizontal",
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.x = element_text(size = 12),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = .7
    )
  ) + guides(fill = guide_legend(ncol=4))

#### SAVE FIGURE S3 ####
ggsave("figS3.commonSTs.clusters.markers.pdf", 
       figS3, 
       path = "outputs/figures", 
       device = "pdf", 
       width = 270,
       height = 180,
       unit ="mm")

#####################################
####---- FIGURE S4 STs MARKER RSTs----####
#####################################
#### FIGURE S4 ####
figS4 <- ggplot(f2_meta, aes(x = Plasmid_Markers)) +
  geom_bar(aes(fill = IncF_simple)) +
  facet_wrap(
    ~ ST,
    nrow = 4,
    strip.position = "top",
    scales = "free_y"
  ) +
  scale_fill_manual(values = alt_F_clrs,
                    name = "F Plasmid RST",
                    guide = guide_legend(ncol = 1)) +
  scale_x_discrete(name = "Plasmid Markers") +
  scale_y_continuous(name = "Count",
                     expand = expansion(mult = c(0, 0.02)),
                     n.breaks = 9,
                     labels = scales::number_format(accuracy = 1, trim = TRUE)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      color = "grey20",
      size = 8,
      vjust = 1,
      hjust = 1,
      face = "plain",
      angle = 30
    ),
    axis.title.x = element_text(size = 14),
    axis.line = element_blank(),
    legend.position = "bottom",
    legend.key.size = unit(5.5, "mm"),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 9),
    legend.direction = "horizontal",
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.x = element_text(size = 12),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = .7
    )
  )+ guides(fill = guide_legend(nrow=2))

#### SAVE FIGURE S4 ####
ggsave("figS4.commonSTs.markers.RSTs.pdf", 
       figS4, 
       path = "outputs/figures", 
       device = "pdf", 
       width = 270,
       height = 180,
       unit ="mm")

#####################################
####     SUPP AND SOURCE DATA    ####
#####################################
# Get input filenames
sourceinfiles <-c(list.files("data/supp_data/", pattern = "\\.csv"))

# Read in accessions for all the collections
for (file in sourceinfiles){
  indir <- c("data/supp_data")
  f <- read_delim(paste(indir, file, sep = "/"), delim = ",", col_names = TRUE, trim_ws = TRUE)
  assign(paste(substr(file, 1, nchar(file)-4), sep = ""), f)
}

# Make dataframe of all the accessions
accessions <- rbind(french,cambridge,oxford,sydney) %>% filter(SRA %in% meta$Name | Name %in% meta$Name)

# Change names so it will bind to meta
accessions <- accessions %>% mutate(Name = ifelse(is.na(SRA), Name, SRA))

# Get relevant colnames
acc_nms <- names(accessions)

# Join meta and mge-cluster tsne coordinates
suppdata1 <- left_join(meta, mge_meta %>% select(Name, tsne1D, tsne2D))

# Join meta and accessions
suppdata1 <- left_join(suppdata1, accessions) %>% select(all_of(acc_nms), everything())

# Get gene screening data
genes <- left_join(args, vags) %>% select(Name, 10:ncol(.))

# Join gene screening to meta and accessions
suppdata1 <- left_join(suppdata1, genes)

# Write it out as Supplementary Data 1
write_csv(suppdata1, "outputs/data/SupplementaryData1.csv")

# Add the list of PLSDB accessions
# Accessions
pls_acc <- read_tsv("data/mge-cluster/models/plasmid_accessions.txt", col_names=c("Accession"))

# PLSDB metadata
pls_meta <- read_tsv("data/mge-cluster/models/plsdb_filtered_F.tsv")

# Join them
plsdb_data <- left_join(pls_acc, pls_meta, by = c("Accession"="NUCCORE_ACC"))

# Save as SuppData2
write_csv(plsdb_data, "outputs/data/SupplementaryData2.csv")      

rm(pls_acc, pls_meta, plsdb_data)
