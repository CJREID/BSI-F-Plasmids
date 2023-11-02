####-----PACKAGES-----####
# Load required packages
library(tidyverse)
library(ggplot2)
library(lemon)
library(reshape2)
library(ggpubr)
library(pheatmap)
library(ggtree)
library(RColorBrewer)
library(abricateR)
library(plasmidmapR)
library(paletteer)
library(viridis)
library(rstatix)

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
common_ST <- meta %>% group_by(ST) %>% filter(n() >= 30) %>% ungroup %>% select(Name) 
intermediate_ST <- meta %>% group_by(ST) %>% filter(n() >=5, n() <30) %>% ungroup %>% select(Name) 
rare_ST <- meta %>% group_by(ST) %>% filter(n() < 5) %>% ungroup %>% select(Name)


# Designate ST ranks
common_ST$ST_rank <- "Common"
intermediate_ST$ST_rank <- "Intermediate"
rare_ST$ST_rank <- "Rare"

# Stick them back together
st_rank <- rbind(common_ST, intermediate_ST) %>% select(Name, ST_rank)
st_rank <- rbind(st_rank, rare_ST) %>% select(Name, ST_rank)

# Join ST ranks to meta
meta <- left_join(meta, st_rank)

# Create a simple ST column where non-common STs are grouped as 'other' 
meta <- meta %>% mutate(ST_simple = case_when(ST_rank == "Common" ~ ST, TRUE ~ "Other"))

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

# Extract metadata column names for later use with the geno_meta dataframe
meta_cols <- c("Name", "ST", "ColV", "IncF_RST", "IncF_simple", "phylogroup", "ST_rank", "ST_simple", "F Plasmid Group")

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
meta <- meta %>% mutate(`F Plasmid Group` = case_when(ColV == "No" & IncF_RST != "F-:A-:B-" & senB_cjrABC == "Yes" ~ "senB",
                                                     ColV == "No" & IncF_RST != "F-:A-:B-" & senB_cjrABC == "No" ~ "Other",
                                                    ColV == "Yes" & IncF_RST != "F-:A-:B-" ~ "ColV",
                                                    IncF_RST == "F-:A-:B-" ~ "None"))

# Re-order F group factor
meta <- meta %>% mutate(`F Plasmid Group` = as.factor(`F Plasmid Group`),
                        `F Plasmid Group` = fct_relevel(`F Plasmid Group`, "None", "ColV", "senB", "Other"))

# Vector of ST order
ST_fct_relevel <- meta %>% group_by(ST) %>% tally(sort=TRUE) %>% pull(ST)

# Order ST factor by ST frequency
meta <- meta %>% mutate(ST =as.factor(ST),
                        ST = fct_relevel(ST, ST_fct_relevel))

# Create 'geno_meta' df with all meta and gene screening data
geno_meta <- left_join(meta %>% select(Name, ST, ColV, IncF_RST, IncF_simple, phylogroup, ST_rank, ST_simple, `F Plasmid Group`), 
                       sepsis_simple_summary_N95L95 %>% select(-ColV, -starts_with("Inc")), by = c("Name" = "name"))

# Define gene columns as those that are integers
gene_cols <- names(geno_meta %>% select(where(is.integer)))

# Recode multiple hits as a single hit
geno_meta <- geno_meta %>% mutate(across(where(is.integer), ~ gsub(2, 1, .x)))
geno_meta <- geno_meta %>% mutate(across(where(is.integer), ~ gsub(3, 1, .x)))
geno_meta <- geno_meta %>% mutate(across(where(is.integer), ~ gsub(4, 1, .x)))

# Convert back to integer (gsub makes everything character class)
geno_meta <- geno_meta %>% mutate(across(all_of(gene_cols), as.integer))

# Meta table with added count so we can do grouped proportional plots
plot_meta <- meta %>% group_by_all() %>% tally(name = "Count")

# Define levels of ST again for some reason idk
ST_levels <- plot_meta %>% filter(ST_rank == "Common") %>% group_by(ST) %>% tally(sort=TRUE) %>% pull(ST)

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
puti89 <- meta %>% group_by(`F Plasmid Group`)
ggplot(puti89 %>% filter(`F Plasmid Group`=="senB") %>% 
         group_by(pUTI89_ID) %>% 
         summarise(count = n()), aes(x=count, y = pUTI89_ID)) +
  geom_col()

# 87.2% of senB+ sequences align to >=70% of pUTI89 at >=90% nucleotide ID

# Plot pCERC4 ID for sequences designated ColV+
pcerc4 <- meta %>% group_by(`F Plasmid Group`)
ggplot(pcerc4 %>% filter(`F Plasmid Group`=="ColV") %>% 
         group_by(pCERC4_ID) %>% 
         summarise(count = n()), aes(x=pCERC4_ID, y = count)) +
  geom_col()

####------GENE PROCESSING-----####
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
         `F Plasmid Group` = factor(`F Plasmid Group`, levels = c("ColV", "senB", "Other", "None"))) %>% 
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
         `F Plasmid Group` = factor(`F Plasmid Group`, levels = c("ColV", "senB", "Other", "None"))) %>% 
  select(all_of(meta_cols), `Total VAGs`)

####-----INTER-TYPE RSTS-----####
# have a look at RSTs that are present in different F plasmid groups
type_incf <- meta %>% group_by(Name, `F Plasmid Group`, IncF_RST) %>% summarise() %>% pivot_wider(names_from = `F Plasmid Group`, values_from = IncF_RST)

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
F_type_colours <- c("ColV" = "#17e6ae", "senB" = "#E7298A", "Other" = "#26e9ff", "None" = "#dbdbdb")

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

# senB data
senB <- geno_meta %>% select(Name, senB = vfdb_senB) %>% mutate(senB = case_when(senB == "1" ~ "senB", senB == "0" ~ "senB-Neg"))

#####################################
####---- FIGURE 1 COMMON STs ----####
#####################################
#### FIGURE 1A ####
# Absolute common ST by F category
fig1a <- ggplot(meta %>% filter(ST_rank == "Common"), aes(x = fct_infreq(factor(ST)))) +
  geom_bar(aes(fill = `F Plasmid Group`), position="stack") +
  scale_fill_manual(values = F_type_colours,
                    guide = guide_legend(title.position = "left")) +
  scale_x_discrete(name = "ST")+
  scale_y_continuous(name = "Count", expand = c(0, 0), limits = c(0,800), n.breaks = 9) +
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
        legend.text = element_text(size = 6),
        legend.direction = "vertical")

#### FIGURE 1B ####
# Make count labels for below plot
f_labs_common <- plot_meta %>% filter(ST_rank == "Common") %>% group_by(ST) %>% tally %>% mutate(`F Plasmid Group`=NA)

# Proportional common ST by F category
fig1b <- ggplot(plot_meta %>% filter(ST_rank == "Common") %>% group_by(ST, `F Plasmid Group`) %>% tally %>% mutate(perc=n/sum(n)), 
                aes(x = factor(ST, levels =ST_levels), y = perc, fill =`F Plasmid Group`)) +
  geom_col()+
  #geom_text(data=f_labs_common,aes(label=n, y = rep(0.97,length(unique(meta %>% filter(ST_rank == "Common") %>% pull(ST)))))) +
  scale_fill_manual(values = F_type_colours,
                    guide = guide_legend(title.position = "left")) +
  scale_x_discrete(name = "ST")+
  scale_y_continuous(name = "Proportion", expand = c(0, 0), limits = c(0, 1), n.breaks = 6) +
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
        legend.text = element_text(size = 6),
        legend.direction = "vertical")

#### FIGURE 1C ####
# Make count labels for below plot
f_labs_rares <- plot_meta %>% filter(ST_rank == "Rare") %>%
  group_by(ST_rank) %>% tally %>% mutate(`F Plasmid Group`=NA)

f_labs_all <- plot_meta %>%
  group_by(ST_rank) %>% tally %>% mutate(`F Plasmid Group`=NA)

phylo_labs <- plot_meta %>%
  group_by(ST_rank) %>% tally %>% mutate(phylogroup=NA)

# Proportional F plasmid type by ST rank
fig1c <- ggplot(plot_meta %>% group_by(ST_rank, `F Plasmid Group`) %>% tally %>% mutate(perc = n/sum(n)), 
                aes(x = fct_infreq(factor(ST_rank)), y = perc, fill = `F Plasmid Group`)) +
  geom_col() +
  geom_text(data=f_labs_all,aes(label=paste0("n=",n), y = rep(0.97,3)), size =3) +
  scale_fill_manual(values = F_type_colours,
                    guide = guide_legend(title.position = "left")) +
  scale_x_discrete(name = "ST Frequency", labels = c("Common (>29)", "Intermediate (5-29)", "Rare (1-4)"))+
  scale_y_continuous(name = "Proportion", expand = c(0, 0), limits = c(0, 1), n.breaks = 6) +
  theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", 
                                   size = 8, 
                                   vjust = 1,
                                   hjust = ,
                                   face = "plain",
                                   angle = ),
        legend.position = "right",
        legend.key.size = unit(2.5, "mm"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        legend.direction = "vertical")

#### FIGURE 1D ####
# Absolute common ST by F RST
fig1d <- ggplot(meta %>% filter(ST_rank == "Common"), aes(x = fct_infreq(factor(ST)))) +
  geom_bar(aes(fill = IncF_simple)) +
  scale_fill_manual(values = alt_F_clrs, name = "F Plasmid RST",
                    guide = guide_legend(nrow=4)) +
  scale_x_discrete(name = "ST")+
  scale_y_continuous(name = "Count", expand = c(0, 0), limits=c(0,800), n.breaks = 9) +
  theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", 
                                   size = 6,
                                   vjust = ,
                                   face = "plain"),
        legend.position = "right",
        legend.key.size = unit(2.5, "mm"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        legend.direction = "horizontal",
        )

#### FIGURE 1E ####
# Make count labels for below plot
rst_lab_common <- plot_meta %>% filter(ST_rank == "Common") %>% group_by(ST) %>% tally %>% mutate(IncF_simple=NA)

# Proportional common ST by F RST
fig1e <- ggplot(plot_meta %>% filter(ST_rank == "Common") %>% group_by(ST, IncF_simple) %>% tally %>% mutate(perc=n/sum(n)), 
                aes(x = factor(ST, levels = ST_levels), y = perc, fill = IncF_simple)) +
  geom_col()+
  # geom_text(data=rst_lab_common, 
   #         aes(label=paste0("n=",n), 
    #        y = rep(0.97,length(unique(meta %>% filter(ST_rank == "Common") %>% pull(ST)))))) +
  scale_fill_manual(values = alt_F_clrs, name = "F Plasmid RST",
                    guide = guide_legend(nrow=4)) +
  scale_x_discrete(name = "ST") +
  scale_y_continuous(name = "Proportion", expand = c(0, 0), n.breaks = 6) +
  theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", 
                                   size = 6,
                                   vjust = ,
                                   face = "plain"),
        legend.position = "right",
        legend.key.size = unit(2.5, "mm"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        legend.direction = "horizontal")

#### FIGURE 1F ####
# Mutate data to make labels possible
Fs <- meta %>% filter(`F Plasmid Group` != "None") %>%
  group_by(`F Plasmid Group`, IncF_simple) %>% tally %>% mutate(perc=n/sum(n))

# F Plasmid category by F RST with labels
fig1f <- ggplot(Fs , aes(x = `F Plasmid Group`, y=perc, fill=IncF_simple)) +
  geom_col(position = "stack") +
  geom_text(aes(label = ifelse(n > 40, IncF_simple, "")), position = position_stack(vjust=0.5), size =1.8) +
  scale_fill_manual(name = "F Plasmid RST", values = alt_F_clrs,
                    guide = guide_legend(nrow=4)) +
  scale_x_discrete(name = "F Plasmid Group", labels = c("ColV (n=924)", "senB (n=1290)", "Other (n=1249)"))+
  scale_y_continuous(name = "Count", expand = c(0, 0), n.breaks = 10) +
  theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", 
                                   size = 9, 
                                   vjust = 1,
                                   hjust = ,
                                   face = "plain",
                                   angle = ), 
        legend.position = "right",
        legend.key.size = unit(2.5, "mm"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        legend.direction = "horizontal")

#### COMPILE AND SAVE FIGURE 1 ####
fig1abc <- ggarrange(plotlist= list(fig1a, fig1b, fig1c),
                     ncol =1, nrow =3, labels =c("a)", "b)","c)"),
                     common.legend = TRUE,
                     legend = "bottom")

fig1def <- ggarrange(plotlist= list(fig1d, fig1e, fig1f),
                     ncol =1, nrow =3, labels =c("d)", "e)","f)"),
                     common.legend = TRUE,
                     legend = "bottom")

fig1 <- ggarrange(plotlist = list(fig1abc, fig1def),
                  ncol=2, legend = "bottom", align = "v")

ggsave("fig1.commonSTs.pdf", 
       fig1, 
       path = "outputs/figures", 
       device = "pdf", 
       width = 297,
       height = 210,
       unit ="mm")

#####################################
####----- FIGURE 2 MAJOR STS ----####
#####################################
#### FIGURE 2 ####
f2_meta <- meta %>% group_by(ST) %>% filter(n() >= 65, `F Plasmid Group` != "None") %>% ungroup() %>% 
  mutate(ST = factor(paste0("ST",ST), levels = paste0("ST", levels(fct_infreq(meta$ST)))))
          
alt_F_clrs2 <- alt_F_clrs[-(length(alt_F_clrs) -1)]

fig2 <- ggplot(f2_meta, aes(x = `F Plasmid Group`)) +
  geom_bar(aes(fill = `IncF_simple`)) +
  facet_rep_wrap(~ ST, nrow = 2, strip.position="top", scales = "free_y", repeat.tick.labels = c("left","bottom")) +
  scale_fill_manual(values = alt_F_clrs2, name = "F Plasmid RST", guide = guide_legend(ncol=1)) +
  scale_x_discrete(name = "F Plasmid Group")+
  scale_y_continuous(name = "Count", expand = expansion(mult = c(0, 0.02)), n.breaks = 9) +
  theme_classic()+
  theme(axis.text.x = element_text(color = "grey20",
                                   size = 12, 
                                   vjust = 1,
                                   hjust = ,
                                   face = "plain",
                                   angle = ),
        axis.title.x = element_text(
          size=14
        ),
        legend.position = "right",
        legend.key.size = unit(5, "mm"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.direction = "vertical",
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text.x = element_text(size = 12),
        panel.border = element_rect(colour = "black", fill = NA, size = .7)
        )

#### SAVE FIGURE 2 ####
ggsave("fig2.STs_plasmids.pdf", 
       fig2, 
       path = "outputs/figures", 
       device = "pdf", 
       width = 297,
       height = 210,
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
#### FIGURE 3 - COLV, HPI, iuc + iro ####
#####################################
#### FIGURE 3 ####
# Infer HPI carriage by presence of fyuA and irp2
HPI <- vags %>% mutate(HPI = case_when(fyuA == "1" & irp2 == "1" ~ "1", TRUE ~ "0")) %>% select(Name, HPI)

# Join ColV, HPI and metadata
colv_gene_meta <- left_join(meta, sepsis.colv.N95.L95 %>% select(-ColV), by = c("Name" = "name")) %>% 
  left_join(., HPI, by = c("Name")) %>% 
  mutate(HPI = as.double(HPI)) %>% 
  mutate(across(where(is.double), ~replace_na(.,0)))

# Filter to select top 15 STs 
ST_filt <- colv_gene_meta %>% group_by(ST) %>%
  tally %>% arrange(desc(n)) %>% slice(1:15) %>% pull(ST)

# Reorder by ST frequency but put Other at the end
ST_order <- colv_gene_meta %>%
  mutate(ST = if_else(ST %notin% ST_filt, "Other", as.character(ST))) %>%
  group_by(ST) %>% 
  tally %>%
  arrange(desc(n)) %>% 
  pull(ST)

ST_order <- c(ST_order[-1], ST_order[1])

# Table for plot of ColV, HPI and operon presence/absence
colv_st_salmo_aero_noHPI <- colv_gene_meta %>%
  mutate(ST = if_else(ST %notin% ST_filt, "Other", as.character(ST))) %>% 
  group_by(ST, HPI, ColV, grp2_iroBCDEN, grp3_iucABCD_iutA) %>% 
  tally %>% 
  select(ST, HPI, ColV, salmo = grp2_iroBCDEN, aero = grp3_iucABCD_iutA, count=n) %>% 
  ungroup %>% 
  mutate(`Salmochelin/Aerobactin` = case_when(
    salmo == 0 & aero == 0 ~ "None",
    salmo == 1 & aero == 0 ~ "Salmochelin",
    salmo == 0 & aero == 1 ~ "Aerobactin",
    salmo == 1 & aero == 1 ~ "Both"
  ),
  ColV_factor = factor(ColV, levels = c("Yes", "No")),
  HPI_factor = factor(ifelse(HPI == 1, "HPI+", "HPI-"), levels = c("HPI+", "HPI-")))

# Alternative plot for senB plasmids
hpi_senB_plot <- colv_gene_meta %>%
  mutate(ST = if_else(ST %notin% ST_filt, "Other", as.character(ST))) %>% 
  group_by(ST, HPI, `F Plasmid Group`, grp2_iroBCDEN, grp3_iucABCD_iutA,senB_cjrABC) %>% 
  tally %>% 
  select(ST, HPI, senB_cjrABC, salmo = grp2_iroBCDEN, aero = grp3_iucABCD_iutA, count=n, `F Plasmid Group`) %>% 
  ungroup %>% 
  mutate(`Salmochelin/Aerobactin` = case_when(
    salmo == 0 & aero == 0 ~ "None",
    salmo == 1 & aero == 0 ~ "Salmochelin",
    salmo == 0 & aero == 1 ~ "Aerobactin",
    salmo == 1 & aero == 1 ~ "Both"
  ),
  Ftype = factor(`F Plasmid Group`, levels = c("ColV", "senB", "Other", "None")),
  HPI_factor = factor(ifelse(HPI == 1, "HPI+", "HPI-"), levels = c("HPI+", "HPI-")))

# Another table including all F groups
hpi_senB_F_plot <- colv_gene_meta %>%
  mutate(ST = if_else(ST %notin% ST_filt, "Other", as.character(ST))) %>% 
  group_by(ST, IncF_simple, HPI, `F Plasmid Group`, grp2_iroBCDEN, grp3_iucABCD_iutA,senB_cjrABC) %>% 
  tally %>% 
  select(ST,IncF_simple, HPI, senB_cjrABC, salmo = grp2_iroBCDEN, aero = grp3_iucABCD_iutA, count=n, `F Plasmid Group`) %>% 
  ungroup %>% 
  mutate(`Salmochelin/Aerobactin` = case_when(
    salmo == 0 & aero == 0 ~ "None",
    salmo == 1 & aero == 0 ~ "Salmochelin",
    salmo == 0 & aero == 1 ~ "Aerobactin",
    salmo == 1 & aero == 1 ~ "Both"
  ),
  Ftype = factor(`F Plasmid Group`, levels = c("ColV", "senB", "Other", "None")),
  HPI_factor = factor(ifelse(HPI == 1, "HPI+", "HPI-"), levels = c("HPI+", "HPI-")))

# Colours
iron_ST_clrs <- c("69" = "#AA4499", "131" = "#795e86","95" = "#9b7121","10" = "#4c94c0","73" = "#DDCC77",
                  "58" = "#f2c8b3","88" = "#CC6677","12" = "#44AA99", "127" = "#e91c04","141" = "#332288",
                  "144" = "#9a8106","38" = "#117733","405" = "#88CCEE","59" = "#73e3a7","80" = "#882255",
                  "Other" = "#c9bfc4")

# set.seed(3)
# iron_ST_clrs2 <- sample(iron_ST_clrs2, length(unique(colv_st_salmo_aero_noHPI$ST)))

# Plot iron system carriage in HPI+ strains by ColV status
fig3 <- ggplot(colv_st_salmo_aero_noHPI, aes(x = fct_relevel(as.factor(`Salmochelin/Aerobactin`), 
                                                             c("Both", "Aerobactin", "Salmochelin", "None")), 
                                             y=count)) + 
  geom_col(aes(fill = ST)) +
  scale_fill_manual(values = iron_ST_clrs, breaks = ST_order, name = "ST") +
  facet_grid(HPI_factor ~ ColV_factor,
             labeller = labeller(ColV_factor = as_labeller(c("Yes" = "ColV+", "No"="ColV-"))),
             switch = "y",
             scales = "fixed")+
  scale_x_discrete(name = "Iron System") +
  scale_y_continuous(
    position = "right",
    name = "Count",
    expand = expansion(mult = c(0, 0.10)),
    n.breaks = 4
  ) +
  theme_classic() +
  theme(
    strip.placement = "outside",
    axis.line.x = element_line(colour="black"),
    axis.text.x = element_text(
      size = 12
    ),
    axis.title.x = element_text(
      size = 16,
    ),
    axis.title.y = element_text(
      size = 16,
    ),
    axis.text.y = element_text(
      size = 12
    ),
    strip.text = element_text(
      size = 16
    ),
    legend.position = "right",
    legend.key.size = unit(6, "mm"),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 10),
    legend.direction = "vertical",
    panel.border = element_rect(colour = "black", fill = NA, size = 1)
  )

ggsave(filename = "fig3.HPI.ColV.pdf", 
       plot = fig3,
       path = "outputs/figures", 
       device = "pdf", 
       width = 297,
       height = 210,
       unit ="mm")



#####################################
####-----FIGURE 4 - TOTAL ARGS+VAGS-----####
#####################################
#### FIGURE 4A ####
# F Plasmid Group vs Total ARGs
F_arg <- ggboxplot(args2, x = "F Plasmid Group", y = "Total ARGs",
                   color = "F Plasmid Group", 
                   order = c("ColV", "senB", "Other", "None"),
                   ylab = "Total ARGs", xlab = "F Plasmid Group") +
  scale_y_continuous(expand = c(0, 0), limits = c(0,30), n.breaks = 6) +
  scale_color_manual(values = F_type_colours) +
  theme_classic() +
  theme(axis.text.x = element_text(color = "grey20", size = 10, vjust = , face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10, vjust = ,face = "plain"),  
        axis.title.x = element_text(color = "grey20",vjust = , size = 12, face = "plain"),
        axis.title.y = element_text(color = "grey20",vjust = , size = 12, face = "plain"),
        legend.position = "none",
        legend.box = "vertical",
        legend.key.size = unit(5, "mm"),
        legend.title = element_text(size=12),
        legend.text = element_text(size = 10))

### Testing for significant differences ###
# Kruskall Wallis tests indicate that there is a difference between ST ranks and F Plasmid Groups in terms of ARG and VAG carriage rates
kruskal.test(`Total ARGs` ~ `F Plasmid Group`, data = args2)

### Plasmid ###
# Pairwise Wilcoxon tests to determine the pairwise differences between BAP groups with BH adjustment for multiple testing
# ARGs
# Mutate column name so wilcox test function doesn't get angry
args2 <- args2 %>% mutate(FPType = `F Plasmid Group`)

# Generate test statistics
arg_plasmid_type <- args2 %>% pairwise_wilcox_test(`Total ARGs` ~ FPType, p.adjust.method = "BH", detailed = TRUE)
arg_plasmid_type <- arg_plasmid_type %>% add_xy_position(x = "F Plasmid Group")
#
# Filter p-values to only display signficant comparisons between BAP2 and other groups
arg_plasmid_type <- arg_plasmid_type %>% filter(p.adj.signif != "ns")

# Replot with BH adjusted p-values on figure
fig4a <- F_arg + 
  stat_pvalue_manual(arg_plasmid_type,
                     label = "p = {p.adj}",
                     label.size = 3,
                     bracket.nudge.y = c(2, 4, 6, 8, 11)) +
  # scale_color_brewer(palette = "Dark2") +
  scale_y_continuous(limits=c(0,42)) +
  theme(legend.position = "none")

#### FIGURE 4B ####
# F Plasmid Group vs Total VAGs
F_vag <- ggboxplot(vags2, x = "F Plasmid Group", y = "Total VAGs",
                   color = "F Plasmid Group", 
                   order = c("ColV", "senB", "Other", "None"),
                   ylab = "Total VAGs", xlab = "F Plasmid Group") +
  scale_y_continuous(expand = c(0, 0), limits = c(0,150), n.breaks = 11) +
  scale_color_manual(values = F_type_colours) +
  theme_classic() +
  theme(axis.text.x = element_text(color = "grey20", size = 10, vjust = ,face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10, vjust = ,face = "plain"),  
        axis.title.x = element_text(color = "grey20",vjust = , size = 12, face = "plain"),
        axis.title.y = element_text(color = "grey20",vjust = , size = 12, face = "plain"),
        legend.position = "none",
        legend.box = "vertical",
        legend.key.size = unit(5, "mm"),
        legend.title = element_text(size=12),
        legend.text = element_text(size = 10))

### Testing for significant differences ###
# Kruskall Wallis tests indicate that there is a difference between ST ranks and F Plasmid Groups in terms of ARG and VAG carriage rates
kruskal.test(`Total VAGs` ~ `F Plasmid Group`, data = vags2)

# Pairwise Wilcoxon tests to determine the pairwise differences between BAP groups with BH adjustment for multiple testing
# VAGs
# Mutate column name so wilcox test function doesn't get angry
vags2 <- vags2 %>% mutate(FPType = `F Plasmid Group`)

# Generate test statistics
vag_plasmid_type <- vags2 %>% pairwise_wilcox_test(`Total VAGs` ~ FPType, p.adjust.method = "BH", detailed = TRUE)
vag_plasmid_type <- vag_plasmid_type %>% add_xy_position(x = "F Plasmid Group") 
#%>% mutate(xmin = case_when(group1 == "ColV" ~ 1,
#                            group1 == "None" ~ 4,
#                            group1 == "Other" ~ 3),
#           xmax = case_when(group2 == "None" ~ 4,
#                         group2 == "Other" ~ 3,
#                         group2 == "senB" ~ 2))
#
# Filter p-values to only display signficant comparisons between BAP2 and other groups
vag_plasmid_type <- vag_plasmid_type %>% filter(p.adj.signif != "ns")

# Replot with BH adjusted p-values on figure
fig4b <- F_vag + 
  stat_pvalue_manual(vag_plasmid_type,
                     label = "p = {p.adj}",
                     label.size = 3,
                     bracket.nudge.y = c(9, 25, 40, 55, 70, 90))+
  # scale_color_brewer(palette = "Dark2") +
  scale_y_continuous(limits=c(0,240)) +
  theme(legend.position = "none")

#### FIGURE 4C ####
# ST Rank vs Total ARGs
ST_arg <- ggboxplot(args2, x = "ST_rank", y = "Total ARGs",
                    color = "ST_rank", 
                    order = c("Common", "Intermediate", "Rare"),
                    ylab = "Total ARGs", xlab = "ST Frequency") +
  scale_y_continuous(expand = c(0, 0), limits = c(0,35), n.breaks = 6) +
  scale_color_brewer(palette = "Dark2") +
  theme_classic() +
  theme(axis.text.x = element_text(color = "grey20", size = 10, vjust = , face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10, vjust = ,face = "plain"),  
        axis.title.x = element_text(color = "grey20",vjust = , size = 12, face = "plain"),
        axis.title.y = element_text(color = "grey20",vjust = , size = 12, face = "plain"),
        legend.position = "none",
        legend.box = "vertical",
        legend.key.size = unit(5, "mm"),
        legend.title = element_text(size=12),
        legend.text = element_text(size = 10))

### Testing for significant differences ###
# Kruskall Wallis tests indicate that there is a difference between ST ranks and F Plasmid Groups in terms of ARG and VAG carriage rates
kruskal.test(`Total ARGs` ~ ST_rank, data = args2)

# Pairwise Wilcoxon tests to determine the pairwise differences between BAP groups with BH adjustment for multiple testing
### ST Rank ###
# Generate test statistics
arg_ST_rank <- args2 %>% mutate(ST_rank = as.character(ST_rank)) %>% 
  pairwise_wilcox_test(`Total ARGs` ~ ST_rank,p.adjust.method = "BH", detailed = TRUE)

arg_ST_rank <- arg_ST_rank %>% add_xy_position(x = "ST_rank")

# Filter p-values to only display signficant comparisons
arg_ST_rank <- arg_ST_rank %>% filter(p.adj.signif != "ns")

# Replot with BH adjusted p-values on figure
fig4c <- ST_arg + 
  stat_pvalue_manual(arg_ST_rank,
                     label = "p = {p.adj}",
                     label.size = 3,
                     bracket.nudge.y = c(.5, 3)) +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "none")

#### FIGURE 4D ####
ST_vag <- ggboxplot(vags2, x = "ST_rank", y = "Total VAGs",
                    color = "ST_rank", 
                    order = c("Common", "Intermediate", "Rare"),
                    ylab = "Total VAGs", xlab = "ST Frequency") +
  scale_y_continuous(expand = c(0, 0), limits = c(0,120), n.breaks = 6) +
  scale_color_brewer(palette = "Dark2") +
  theme_classic() +
  theme(axis.text.x = element_text(color = "grey20", size = 10, vjust = , face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10, vjust = ,face = "plain"),  
        axis.title.x = element_text(color = "grey20",vjust = , size = 12, face = "plain"),
        axis.title.y = element_text(color = "grey20",vjust = , size = 12, face = "plain"),
        legend.position = "none",
        legend.box = "vertical",
        legend.key.size = unit(5, "mm"),
        legend.title = element_text(size=12),
        legend.text = element_text(size = 10))

### Testing for significant differences ###
# Kruskall Wallis tests indicate that there is a difference between ST ranks and F Plasmid Groups in terms of ARG and VAG carriage rates
kruskal.test(`Total VAGs` ~ ST_rank, data = vags2)

# Pairwise Wilcoxon tests to determine the pairwise differences between BAP groups with BH adjustment for multiple testing
### ST Rank ###
# Generate test statistics
vag_ST_rank <- vags2 %>% mutate(ST_rank = as.character(ST_rank)) %>% 
  pairwise_wilcox_test(`Total VAGs` ~ ST_rank,p.adjust.method = "BH", detailed = TRUE)

vag_ST_rank <- vag_ST_rank %>% add_xy_position(x = "ST_rank")

# Filter p-values to only display signficant comparisons
vag_ST_rank <- vag_ST_rank %>% filter(p.adj.signif != "ns")

# Replot with BH adjusted p-values on figure
fig4d <- ST_vag + 
  stat_pvalue_manual(vag_ST_rank,
                     label = "p = {p.adj}",
                     label.size = 3,
                     bracket.nudge.y = c(4, 14, 27)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,170), n.breaks = 6) +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "none")


#### FIGURE 4E ####
# Common STs vs Total ARGs
fig4e <- ggboxplot(args2 %>% filter(ST_simple != "Other"), x = "ST_simple", y = "Total ARGs",
                   color = "ST_simple", 
                   order = c(levels(fct_infreq(factor(args2$ST)))),
                   ylab = "Total ARGs", xlab = "ST") +
  scale_y_continuous(expand = c(0, 0), limits = c(0,26), n.breaks = 6) +
  scale_color_manual(values = ST_clrs[2:length(ST_clrs)]) +
  #scale_x_discrete(labels = c(paste0("ST", levels(fct_infreq(factor(vags2$ST_simple)))[2:length(levels(fct_infreq(factor(vags2$ST_simple))))]))) +
  theme_classic() +
  theme(axis.text.x = element_text(color = "grey20", size = 10, vjust = , face = "plain", angle = 35, hjust = 1),
        axis.text.y = element_text(color = "grey20", size = 10, vjust = ,face = "plain"),  
        axis.title.x = element_text(color = "grey20",vjust = , size = 12, face = "plain"),
        axis.title.y = element_text(color = "grey20",vjust = , size = 12, face = "plain"),
        legend.position = "none",
        legend.box = "vertical",
        legend.key.size = unit(5, "mm"),
        legend.title = element_text(size=12),
        legend.text = element_text(size = 10))

#### FIGURE 4F ####
# Common STs vs Total VAGs
fig4f <- ggboxplot(vags2 %>% filter(ST_simple != "Other"), x = "ST_simple", y = "Total VAGs",
                   color = "ST_simple", 
                   order = c(levels(fct_infreq(factor(args2$ST)))),
                   ylab = "Total VAGs", xlab = "ST") +
  scale_y_continuous(expand = c(0, 0), limits = c(20,140), n.breaks = 11) +
  scale_color_manual(values = ST_clrs[2:length(ST_clrs)]) +
  #scale_x_discrete(labels = c(paste0("ST", levels(fct_infreq(factor(vags2$ST_simple)))[2:length(levels(fct_infreq(factor(vags2$ST_simple))))]))) +
  theme_classic() +
  theme(axis.text.x = element_text(color = "grey20", size = 10, vjust = , face = "plain", angle = 35, hjust =1),
        axis.text.y = element_text(color = "grey20", size = 10, vjust = ,face = "plain"),  
        axis.title.x = element_text(color = "grey20",vjust = , size = 12, face = "plain"),
        axis.title.y = element_text(color = "grey20",vjust = , size = 12, face = "plain"),
        legend.position = "none",
        legend.box = "vertical",
        legend.key.size = unit(5, "mm"),
        legend.title = element_text(size=12),
        legend.text = element_text(size = 10))

# Intermediate STs vs TotalARGs
# ggboxplot(args2 %>% filter(ST_rank == "Intermediate"), x = "ST", y = "Total ARGs",
#           order = c(levels(fct_infreq(factor(args2$ST)))),
#           ylab = "", xlab = "ST") +
#   scale_y_continuous(expand = c(0, 0), limits = c(0,30), n.breaks = 5) +
#   theme_classic() +
#   theme(axis.text.x = element_text(color = "grey20", size = 10, vjust = , hjust = .9, face = "plain", angle = 30),
#         axis.text.y = element_text(color = "grey20", size = 10, vjust = ,face = "plain"),  
#         axis.title.x = element_text(color = "grey20",vjust = , size = 12, face = "plain"),
#         axis.title.y = element_text(color = "grey20",vjust = , size = 12, face = "plain"),
#         legend.position = "none",
#         legend.box = "vertical",
#         legend.key.size = unit(5, "mm"),
#         legend.title = element_text(size=12),
#         legend.text = element_text(size = 10))


#### COMPILE FIGURE 4 ####
fig4ab <- ggarrange(plotlist= list(fig4a, fig4b),
                    ncol =2, labels =c("a)", "b)"),
                    legend = "none")

fig4cd <- ggarrange(plotlist= list(fig4c, fig4d),
                    ncol =2, labels =c("c)", "d)"),
                    legend = "none")

fig4ef<- ggarrange(plotlist= list(fig4e, fig4f),
                   ncol =2, labels =c("e)", "f)"),
                   legend = "none")

fig4 <- ggarrange(plotlist = list(fig4ab, fig4cd, fig4ef),
                  nrow =3)

ggsave("fig4.args.vags.pdf", 
       fig4, 
       path = "outputs/figures", 
       device = "pdf", 
       width = 297,
       height = 210,
       unit ="mm")

#####################################
####--- FIGURE S1 PLASMID GROUPS ---####
#####################################
#### FIGURE S1A ####
figs1a <- ggplot(plot_meta %>% group_by(Count, `F Plasmid Group`) %>% tally %>% mutate(perc = n/sum(n)), aes(x = Count, y = perc, fill = `F Plasmid Group`))+
  geom_col()+
  scale_fill_manual(values = F_type_colours, name = "F Plasmid Group",
                    guide = guide_legend(title.position = "left")) +
  scale_x_discrete(name = "All Genomes")+
  scale_y_continuous(name = "Proportion", expand = c(0, 0), limits=c(0,1), n.breaks = 11) +
  theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", 
                                   size = 6,
                                   vjust = ,
                                   face = "plain"),
        legend.position = "right",
        legend.key.size = unit(4, "mm"),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        legend.direction = "horizontal")

#### FIGURE S1B ####
figs1b <-ggplot(meta %>% filter(phylogroup != "EC_control_fail"),
                aes(x = phylogroup)) +
  geom_bar(aes(fill = `F Plasmid Group`)) +
  scale_fill_manual(values = F_type_colours, name = "F Plasmid Group",
                    guide = guide_legend(title.position = "left")) +
  scale_x_discrete(name = "Phylogroup") +
  scale_y_continuous(
    name = "Count",
    expand = c(0, 0),
    limits = c(0, 3100),
    n.breaks = 11
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      color = "grey20",
      size = 8,
      vjust = ,
      face = "plain"
    ),
    legend.position = "right",
    legend.key.size = unit(4, "mm"),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    legend.direction = "horizontal"
  )

#### FIGURE S1C ####
# SUPP FIGURE Proportional plasmid type by phylogroup
figs1c <- ggplot(
  plot_meta %>% filter(phylogroup != "EC_control_fail") %>% group_by(phylogroup, `F Plasmid Group`) %>% tally %>% mutate(perc = n /
                                                                                                                           sum(n)),
  aes(x = phylogroup, y = perc, fill = `F Plasmid Group`)
) +
  geom_col() +
  scale_fill_manual(values = F_type_colours, name = "F Plasmid Group",
                    guide = guide_legend(title.position = "left")) +
  scale_x_discrete(name = "Phylogroup") +
  scale_y_continuous(
    name = "Count",
    expand = c(0, 0),
    limits = c(0, 1),
    n.breaks = 11
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      color = "grey20",
      size = 8,
      vjust = ,
      face = "plain"
    ),
    legend.position = "right",
    legend.key.size = unit(4, "mm"),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    legend.direction = "horizontal"
  )


#### FIGURE S1D ####
# Absolute Intermediate ST by F category
figs1d <- ggplot(meta %>% filter(ST_rank == "Intermediate"), aes(x = fct_infreq(factor(ST)))) +
  geom_bar(aes(fill = `F Plasmid Group`)) +
  scale_fill_manual(values = F_type_colours) +
  scale_x_discrete(name = "ST")+
  scale_y_continuous(name = "Count", expand = c(0, 0), n.breaks = 8) +
  theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", 
                                   size = 5.5, 
                                   vjust = ,
                                   hjust = 1,
                                   face = "plain",
                                   angle = 45), 
        legend.position = "right",
        legend.key.size = unit(4, "mm"),
        legend.title = element_text(size =8),
        legend.text = element_text(size = 6),
        legend.direction = "horizontal")


#### FIGURE S1E ####
# Make count labels for below plot
f_labs_inter <- plot_meta %>% filter(ST_rank == "Intermediate") %>% group_by(ST) %>% tally %>% mutate(`F Plasmid Group`=NA)

# Proportional ST by F category
figs1e <- ggplot(plot_meta %>% filter(ST_rank == "Intermediate"), aes(x = fct_infreq(factor(ST)), y = Count)) +
  geom_bar(aes(fill = `F Plasmid Group`), position = "fill", stat = "identity") +
  # geom_text(data=f_labs_inter,aes(label=n, y = rep(0.97,length(unique(meta %>% filter(ST_rank == "Intermediate") %>% pull(ST)))))) +
  scale_fill_manual(values = F_type_colours) +
  scale_x_discrete(name = "ST")+
  scale_y_continuous(name = "Count", expand = c(0, 0), limits = c(0, 1), n.breaks = 8) +
  theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", 
                                   size = 5.5, 
                                   vjust = ,
                                   hjust = 1,
                                   face = "plain",
                                   angle = 45),
        legend.position = "right",
        legend.key.size = unit(4, "mm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.direction = "horizontal")


#### COMPILE AND SAVE FIGURE S1 ####
figs1abc <- ggarrange(plotlist= list(figs1a, figs1b, figs1c),
                      ncol =1, nrow =3, labels =c("a)", "b)", "c)"),
                      common.legend = TRUE,
                      legend = "bottom")

figs1de <- ggarrange(plotlist= list(figs1d, figs1e),
                     ncol =1, nrow =2, labels =c("d)", "e)"),
                     common.legend = TRUE,
                     legend = "none")

figs1 <- ggarrange(plotlist = list(figs1abc, figs1de),
                   ncol=2, common.legend =TRUE, legend = "bottom", align = "v")

# figs1_alt <- ggarrange(plotlist = list(figs1, figSX),
#                     ncol =1, align = "v")

ggsave("figS1.plasmid_groups.pdf", 
       figs1, 
       path = "outputs/figures", 
       device = "pdf", 
       width = 297,
       height = 210,
       unit ="mm")

#####################################
#### FIGURE S2 - PHYLOGROUP ####
#####################################

## Phylogroup ST Rank
figS2a <- ggplot(
  plot_meta %>% filter(phylogroup != "EC_control_fail") %>% group_by(ST_rank,phylogroup) %>% tally %>% mutate(perc = n /sum(n)),
  aes(x = ST_rank, y = perc, fill = phylogroup)
) +
  geom_col() +
  scale_fill_manual(values = phylogroup_clrs, name = "Phylogroup",
                    guide = guide_legend(title.position = "left")) +
  scale_x_discrete(name = "ST Frequency") +
  scale_y_continuous(
    name = "Count",
    expand = c(0, 0),
    limits = c(0, 1.0001),
    n.breaks = 11
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      color = "grey20",
      size = 9,
      vjust = ,
      face = "plain"
    ),
    legend.position = "right",
    legend.key.size = unit(4, "mm"),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    legend.direction = "horizontal"
  )

## Phylogroup Plasmid group
figS2b <- ggplot(
  plot_meta %>% filter(phylogroup != "EC_control_fail") %>% group_by(`F Plasmid Group`,phylogroup) %>% tally %>% mutate(perc = n /sum(n)),
  aes(x = `F Plasmid Group`, y = perc, fill = phylogroup)
) +
  geom_col() +
  scale_fill_manual(values = phylogroup_clrs, name = "Phylogroup",
                    guide = guide_legend(title.position = "left")) +
  scale_x_discrete(name = "F Plasmid Group", limits = c("ColV", "senB", "Other", "None")) +
  scale_y_continuous(
    name = "Count",
    expand = c(0, 0),
    limits = c(0, 1.0001),
    n.breaks = 11
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      color = "grey20",
      size = 9,
      vjust = ,
      face = "plain"
    ),
    legend.position = "right",
    legend.key.size = unit(4, "mm"),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    legend.direction = "horizontal"
  )

figS2 <- ggarrange(plotlist = list(figS2a, figS2b),
                   labels = c("a)", "b)"),
                   ncol=2, legend = "bottom", align = "v", common.legend = TRUE)

ggsave("figS2.phylogroups.pdf", 
       figS2, 
       path = "outputs/figures", 
       device = "pdf", 
       width = 297,
       height = 120,
       unit ="mm")


#####################################
#### FIGURE S3 RSTs ####
#####################################
#### FIGURE S3A ####
figs3a <- ggplot(plot_meta %>% group_by(Count, `F Plasmid Group`, IncF_simple) %>% tally %>% mutate(perc = n/nrow(meta)), aes(x = Count, y = perc, fill = IncF_simple))+
  geom_col()+
  scale_fill_manual(values = alt_F_clrs, name = "F Plasmid RST",
                    guide = guide_legend(nrow=4)) +
  scale_x_discrete(name = "All Genomes")+
  scale_y_continuous(name = "Proportion", expand = c(0, 0), limits=c(0,1), n.breaks = 11) +
  theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", 
                                   size = 8,
                                   vjust = ,
                                   face = "plain"),
        legend.position = "right",
        legend.key.size = unit(4, "mm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 6),
        legend.direction = "horizontal")

#### FIGURE S3B ####
figs3b <- ggplot(meta %>% filter(phylogroup != "EC_control_fail"),
                 aes(x = phylogroup)) +
  geom_bar(aes(fill = `IncF_simple`)) +
  scale_fill_manual(values = alt_F_clrs, name = "F Plasmid RST",
                    guide = guide_legend(nrow=4)) +
  scale_x_discrete(name = "Phylogroup") +
  scale_y_continuous(
    name = "Count",
    expand = c(0, 0),
    n.breaks = 11
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      color = "grey20",
      size = 8,
      vjust = ,
      face = "plain"
    ),
    legend.position = "right",
    legend.key.size = unit(4, "mm"),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 6),
    legend.direction = "horizontal"
  )

#### FIGURE S3C ####
figs3c <- ggplot(
  plot_meta %>% filter(phylogroup != "EC_control_fail") %>% group_by(phylogroup, `IncF_simple`) %>% tally %>% mutate(perc = n /
                                                                                                                       sum(n)),
  aes(x = phylogroup, y = perc, fill = IncF_simple)
) +
  geom_col() +
  scale_fill_manual(values = alt_F_clrs, name = "F Plasmid RST",
                    guide = guide_legend(nrow=4)) +
  scale_x_discrete(name = "Phylogroup") +
  scale_y_continuous(
    name = "Count",
    expand = c(0, 0),
    limits = c(0, 1),
    n.breaks = 11
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      color = "grey20",
      size = 8,
      vjust = ,
      face = "plain"
    ),
    legend.position = "right",
    legend.key.size = unit(4, "mm"),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 6),
    legend.direction = "horizontal"
  )
#### FIGURE S3D ####
# Absolute intermediate ST by F RST
figs3d <- ggplot(meta %>% filter(ST_rank == "Intermediate"), aes(x = fct_infreq(factor(ST)))) +
  geom_bar(aes(fill = IncF_simple)) +
  scale_fill_manual(values = alt_F_clrs, name = "F Plasmid RST") +
  scale_x_discrete(name = "ST")+
  scale_y_continuous(name = "Count", expand = c(0, 0), n.breaks = 6) +
  theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", 
                                   size = 5.5,
                                   vjust = ,
                                   hjust = 1,
                                   face = "plain",
                                   angle =45),
        legend.position = "right",
        legend.key.size = unit(4, "mm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.direction = "horizontal")

#### FIGURE S3E ####
# Proportional intermediate ST by F RST
figs3e <- ggplot(plot_meta %>% filter(ST_rank == "Intermediate"), aes(x = fct_infreq(factor(ST)), y = Count)) +
  geom_bar(aes(fill = IncF_simple), stat = "identity", position = "fill") +
  # geom_text(data=rst_lab_inter,aes(label=paste0("n=",n), y = rep(0.97,length(unique(meta %>% filter(ST_rank == "Intermediate") %>% pull(ST))))), size =2.5) +
  scale_fill_manual(values = alt_F_clrs, name = "F Plasmid RST") +
  scale_x_discrete(name = "ST")+
  scale_y_continuous(name = "Count", expand = c(0, 0), n.breaks = 6) +
  theme_classic()+
  theme(axis.text.x = element_text(color = "grey20", 
                                   size = 5.5,
                                   vjust = ,
                                   hjust = 1,
                                   face = "plain",
                                   angle = 45),
        legend.position = "right",
        legend.key.size = unit(4, "mm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.direction = "horizontal") 

#### COMPILE AND SAVE FIGURE S3 ####
figs3abc <- ggarrange(plotlist= list(figs3a, figs3b, figs3c),
                      ncol =1, labels =c("a)", "b)", "c)"),
                      common.legend = TRUE,
                      legend = "bottom")

figs3de <- ggarrange(plotlist= list(figs3d, figs3e),
                     ncol =1, labels =c("c)", "d)"),
                     legend = "none")


figs3 <- ggarrange(plotlist = list(figs3abc, figs3de),
                   ncol =2)

ggsave("figS3.RSTs.pdf", 
       figs3, 
       path = "outputs/figures", 
       device = "pdf", 
       width = 297,
       height = 210,
       unit ="mm")


#####################################
#### FIGURE S4 - EXTRA STs PLASMIDS ####
#####################################
otherST_meta <- meta %>% group_by(ST) %>% filter(n() < 65, `F Plasmid Group` != "None", ST_rank == "Common") %>% ungroup() %>% 
  mutate(ST = factor(paste0("ST",ST), levels = paste0("ST", levels(fct_infreq(meta$ST)))))

figS4 <- ggplot(otherST_meta, aes(x = `F Plasmid Group`)) +
  geom_bar(aes(fill = `IncF_simple`)) +
  facet_rep_wrap(~ ST, nrow = 2, strip.position="top", scales = "free_y", repeat.tick.labels = c("left","bottom")) +
  scale_fill_manual(values = alt_F_clrs, name = "F Plasmid RST", guide = guide_legend(ncol=1)) +
  scale_x_discrete(name = "F Plasmid Group")+
  scale_y_continuous(name = "Count", expand = c(0, 0), limits = c(0,50)) +
  theme_classic()+
  theme(axis.text.x = element_text(color = "grey20",
                                   size = 8, 
                                   vjust = 1,
                                   hjust = ,
                                   face = "plain",
                                   angle = ),
        legend.position = "right",
        legend.key.size = unit(3, "mm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.direction = "vertical",
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text.x = element_text(size = 10))


#### SAVE FIGURE S4 ####
ggsave("figS4.otherSTs_plasmids.pdf", 
       figS4, 
       path = "outputs/figures", 
       device = "pdf", 
       width = 297,
       height = 210,
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
# Join meta and accessions
suppdata1 <- left_join(meta, accessions) %>% select(acc_nms, everything())
# Get gene screening data
genes <- left_join(args, vags) %>% select(Name, 10:ncol(.))
# Join gene screening to meta and accessions
suppdata1 <- left_join(suppdata1, genes)

# Write it out as Supplementary Data 1
write_csv(suppdata1, "outputs/data/SupplementaryData1.csv")

## Create csvs for Nature Comms SourceData file by figure
# Source data - Fig. 1, 1, S1-4 (plot_meta from BSI_F_plasmids.r script)
write_csv(plot_meta, "outputs/data/src_fig1-2_S1-4.csv")

# Source data - Fig. 3
write_csv(colv_st_salmo_aero_noHPI, "outputs/data/src_fig3.csv")

# Source data - Fig. 4
src_fig4 <- left_join(args2, vags2)
write_csv(src_fig4, "outputs/data/src_fig4.csv")

## I compiled these csvs manually in Excel for the final file
