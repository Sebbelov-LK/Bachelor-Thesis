
# Full Result file

# Libraries that have been in use
library(tidyverse)
library(tidygraph)
library(ggthemes)
library(patchwork)
library(ggraph)
library(igraph)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(readxl)       
library(enrichR)      
library(ggforce)      
library(viridis)      
library(forcats)
library(readr)
library(broom)
library(Peptides)
library(rstatix)
library(ggpubr)
library(scales)
library(stringr)
library(treemapify)
library(plotly)
library(gprofiler2)


## Protein matches and cohort distinction

autoimmune_data <- read_csv("autoimmune_iedb_pepetides.csv")
donor_hlas <- read_tsv("Inhouse_donor_hlas.tsv")
donor_ligands <- read_csv("supplementary_table1.csv")
head(autoimmune_data)
head(donor_hlas)
head(donor_ligands)

### Test for the number of matches between the two csv files
### Due to the suffix of the supplementary taple there were originally no matches, meaning we had to remove those in "clean"

autoimmune_data_clean <- autoimmune_data %>%
  mutate(proteinID_clean = str_remove(`Epitope: Parent Protein Accession`, "\\.\\d+$"))

donor_ligands_clean <- donor_ligands %>%
  mutate(proteinID_clean = str_remove(Source_proteinID, "\\.\\d+$"))

# weerform an inner join on the newly created 'clean' columns
protein_matches <- inner_join(
  autoimmune_data_clean,
  donor_ligands_clean,
  by = "proteinID_clean"
)

# and print out the number of matching rows and the first few matches
cat("Number of matching rows found:", nrow(protein_matches), "\n")
print(head(protein_matches, 10))

# We found 115418 matching rows
# We here now have a dataset with 115418 observations, with nine variables. Below is established a list with rows that did not actually match anything.

autoimmune_unmatched <- anti_join(
  autoimmune_data_clean,
  donor_ligands_clean,
  by = "proteinID_clean"
)
cat("Number of autoimmune_data rows without a match:", nrow(autoimmune_unmatched), "\n")

# Rows without a match, autoimmune: 1069

donor_unmatched <- anti_join(
  donor_ligands_clean,
  autoimmune_data_clean,
  by = "proteinID_clean"
)
cat("Number of donor_ligands rows without a match:", nrow(donor_unmatched), "\n")

# Number of donor_ligands rows without a match: 361570

# Now we can take a look at the distinct ones, to see how many there are. We can also group and nest the matched rows

protein_matches_distinct <- protein_matches %>%
  distinct()

protein_matches_nested <- protein_matches_distinct %>%
  group_by(proteinID_clean) %>%
  nest() %>%
  # Adding a column to show how many matched rows each protein has
  mutate(num_matches = map_int(data, nrow))

cat("Number of distinct proteins matched:", nrow(protein_matches_nested), "\n")
print(protein_matches_nested, n = 10)

## We see that there are 47 distinct proteins from the autoimmune data that matches with the cleaned proteins from the healthy dataset

# The distinct amount of proteins that actually are in the autoimmune dataset:

distinct_autoimmune_proteins <- autoimmune_data %>%
  distinct(`Epitope: Parent Protein Accession`)


num_distinct_proteins <- nrow(distinct_autoimmune_proteins)

cat("Number of distinct proteins in the autoimmune dataset:", num_distinct_proteins, "\n")

# 116, meaning there are 69 distinct proteins with no matches

autoimmune_unmatched_distinct <- autoimmune_unmatched %>%
  distinct(`Epitope: Parent Protein Accession`)

num_unmatched_distinct <- nrow(autoimmune_unmatched_distinct)

cat("Number of distinct unmatched proteins from autoimmune:",
    num_unmatched_distinct, "\n")

# This shows us that 69 distinct proteins do not have matches

#___________________________________________________________________________________________________________________

# Protein Analysis
## Working from this we pull metadata and FASTA files from UniProt after having found all the proteins. 
## We create three files, one for the matched, one for the unmatched and one full.

matched <- read_delim("protein_data/uniprotkb_matched_proteins.csv", 
                      delim = ";",  # specify semicolon as delimiter
                      escape_double = FALSE, 
                      trim_ws = TRUE)
## Rows: 47, columns: 11

unmatched <- read_delim("protein_data/uniprotkb_unmatched_proteins.csv", 
                        delim = ";", 
                        escape_double = FALSE, 
                        trim_ws = TRUE)
## Rows: 69, columns: 11

full_data <- read_delim("protein_data/uniprotkb_full_protein_list.csv", 
                        delim = ";", 
                        escape_double = FALSE, 
                        trim_ws = TRUE)
## Rows: 116, columns: 50

# Check column names
colnames(matched)
colnames(unmatched)
colnames(full_data)

# We have taken the following columns from UniProt to the two smaller files: 
# [1] "Entry"                     "Reviewed"                 
# [3] "Entry Name"                "Protein names"            
# [5] "Gene Names"                "Organism"                 
# [7] "Length"                    "Subcellular location [CC]"
# [9] "Function [CC]"             "DNA binding"              
# [11] "Sequence similarities"    

# For the file containing all the proteins we have taken as much as 50 columns if needed for testing

# We merge and annotate them

# First we create a vector of all relevant protein IDs from the matched and unmatched data
relevant_entries <- c(matched$Entry, unmatched$Entry)

# then subset the full data to include only these proteins
full_data_subset <- full_data %>%
  filter(Entry %in% relevant_entries)

# Annotate each protein with a "match_status" column
full_data_subset <- full_data_subset %>%
  mutate(match_status = case_when(
    Entry %in% matched$Entry    ~ "Healthy & Autoimmune",
    Entry %in% unmatched$Entry  ~ "Autoimmune Only",
    TRUE                        ~ "Other"  # In case of unexpected overlap
  ))

# Let us do a quick check the distribution of match_status
table(full_data_subset$match_status)

# This is the results:
# Autoimmune Only Healthy & Autoimmune 
#             69                   47 

## Starting of with an exploratory check of length

ggplot(full_data_subset, aes(x = match_status, y = Length, fill = match_status)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.4) +
  labs(x = "", 
       y = "Protein Length", 
       fill = NULL) +
  theme_minimal() +
  theme(
    legend.text  = element_text(size = 12),
    legend.title = element_text(size = 12)
  )


## There is little difference
## Let us quickly look at means, medians, along with other statistics of length for each of the groups

full_data_subset %>%
  group_by(match_status) %>%
  summarise(
    count = n(),
    mean_length = mean(Length, na.rm = TRUE),
    median_length = median(Length, na.rm = TRUE),
    sd_length = sd(Length, na.rm = TRUE)
  )

#match_status       count mean_length median_length sd_length
#     <chr>              <int>       <dbl>         <dbl>     <dbl>
# 1 Autoimmune Only       69        644.           462      652.
# 2 Healthy & Autoimm…    47        664.           454      795.

# To assess if the difference in protein length between the two groups is statistically significant, 
# we can run a Wilcoxon rank-sum test (nonparametric) or a t-test (parametric, assuming normality):

df_two_groups <- full_data_subset %>%
  filter(match_status %in% c("Autoimmune Only", "Healthy & Autoimmune"))

wilcox.test(Length ~ match_status, data = df_two_groups)

#Wilcoxon rank sum test with continuity correction

# data:  Length by match_status
# W = 1707.5, p-value = 0.6306
# alternative hypothesis: true location shift is not equal to 0

#_____________________________________________________________________________________________________________________


## We took a bit of data from UniProt, one also being the sequence. We would therefore like to make a test regarding aa composition and net charge
## First we make sure that we actually have a sequence column.

df_seq <- full_data_subset %>%
  filter(!is.na(Sequence) & Sequence != "")

cat("Number of proteins with sequences:", nrow(df_seq), "\n")

# Number of proteins with sequences: 116 

# Perfect, now let us begin calculating the approximate net charge. We’ll define a function that:
# Counts Asp (D) and Glu (E) as -1 each.
# Counts Lys (K) and Arg (R) as +1 each.
# Gives His (H) a partial charge of +0.1 (a rough approximation at pH \~7.4, since histidine’s side chain pKa is around 6.0).
# If we see something noteworthy, we will try with more precise calculations.

calc_net_charge <- function(aa_sequence) {
  s <- toupper(aa_sequence)   # this to ensure uppercase
  nD <- str_count(s, "D")
  nE <- str_count(s, "E")
  nK <- str_count(s, "K")
  nR <- str_count(s, "R")
  nH <- str_count(s, "H")
  
  # Approx net charge ignoring partial charges except for H
  # Lys + Arg = +1 each, His = +0.1, Asp + Glu = -1 each
  net_charge <- (nK + nR + 0.1 * nH) - (nD + nE)
  
  return(net_charge)
}

# Now we apply it to all sequences
# We’ll create two new columns in df_seq:
# net_charge_approx: The numeric net charge estimate.
# charge_class: A label: “positive,” “negative,” or “neutral."

df_seq_charge <- df_seq %>%
  mutate(
    net_charge_approx = map_dbl(Sequence, calc_net_charge),
    charge_class = case_when(
      net_charge_approx > 0  ~ "positive",
      net_charge_approx < 0  ~ "negative",
      TRUE                   ~ "neutral"  # for exactly 0
    )
  )

head(df_seq_charge %>% select(Entry, net_charge_approx, charge_class, match_status))

# Entry      net_charge_approx charge_class match_status   
# <chr>                  <dbl> <chr>        <chr>          
# 1 A0A3B3IU09             0.400 positive     Autoimmune Only
# 2 A0A3B3IU43             6     positive     Autoimmune Only
# 3 A0A494C0V5             3     positive     Autoimmune Only
# 4 D6RBI3                10.2   positive     Autoimmune Only
# 5 E9PKL9                -7.9   negative     Autoimmune Only
# 6 F6XSS0                16.3   positive     Autoimmune Only

# Let us compare the charge distribution using a boxplot.

charge1_fig <- ggplot(df_seq_charge, aes(x = match_status, y = net_charge_approx, fill = match_status)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.4) +
  labs(
    x = "",
    y = "Net Charge (approx)",
    fill = NULL
  ) +
  theme_minimal() +
  theme(
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12)
  )

# For now, it doesn't appear to show us that much. Let us do a statistical test. Once again using Wilcoxon rank-sum test.

df_two_groups <- df_seq_charge %>%
  filter(match_status %in% c("Autoimmune Only", "Healthy & Autoimmune"))

wilcox.test(net_charge_approx ~ match_status, data = df_two_groups)

# Wilcoxon rank sum test with continuity correction
# data:  net_charge_approx by match_status
# W = 1546.5, p-value = 0.6752
# alternative hypothesis: true location shift is not equal to 0

# Unfortunately We do not get a statistically significant p-value. We can try again in a more refined manner using the Peptides package.

df_seq_charge_peptides <- df_seq %>%
  mutate(
    net_charge_peptides = map_dbl(Sequence, ~ charge(., pH = 7.4, pKscale = "EMBOSS"))
  )

head(df_seq_charge_peptides %>% select(Entry, net_charge_peptides))

# Entry      net_charge_peptides
# <chr>                    <dbl>
# 1 A0A3B3IU09              -0.111
# 2 A0A3B3IU43               5.82 
# 3 A0A494C0V5               2.76 
# 4 D6RBI3                   8.94 
# 5 E9PKL9                  -8.06 
# 6 F6XSS0                  15.9 

charge2_fig <- ggplot(df_seq_charge_peptides, 
       aes(x = match_status, y = net_charge_peptides, fill = match_status)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.4) +
  labs(
    x = "",
    y = "Net Charge",
    fill = NULL
  ) +
  theme_minimal() +
  theme(
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12)
  )

# The figure looks very similar and still shows us nothing

df_two_groups <- df_seq_charge_peptides %>%
  filter(match_status %in% c("Autoimmune Only", "Healthy & Autoimmune"))

wilcox.test(net_charge_peptides ~ match_status, data = df_two_groups)

# Wilcoxon rank sum test with continuity correction
# data:  net_charge_peptides by match_status
# W = 1546, p-value = 0.6732
# alternative hypothesis: true location shift is not equal to 0

# We save and print them together
combined_charge_fig <- (charge1_fig + charge2_fig) +
  plot_layout(ncol = 2, guides = "collect") & 
  theme(
    legend.position = "bottom",
    legend.text     = element_text(size = 16),
    legend.title    = element_text(size = 16)
  )

ggsave(
  filename = "results/combined_charge_plots.png",
  plot     = combined_charge_fig,
  width    = 12,    # adjust as needed (in inches)
  height   = 6,
  dpi      = 400
)

# From what I can see, the only great difference here is the median-charge between the two groups.

#_________________________________________________________________________________________________________

## We now focus on the amino acid composition

# We first define a helper function that takes an amino acid sequence (as a string) and returns a named vector of the fractions for each of the 20 standard amino acids.

# Naturally we first need to define the 20 standard amino acids
standard_aa <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", 
                 "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")

# then comes the function to calculate composition fractions
aa_composition <- function(seq) {
  # remember to ensure sequence is uppercase and split into individual amino acids
  aa <- unlist(strsplit(toupper(seq), split = ""))
  total <- length(aa)
  # we need a count for each amino acid present in the sequence
  counts <- table(aa)
  # Then we can calculate the fraction for each standard amino acid, making sure that the length of the proteins doesn't ruin the result
  comp <- setNames(rep(0, length(standard_aa)), standard_aa)
  comp[names(counts)] <- as.numeric(counts) / total
  return(comp)
}

# We can do a quick test of the function on a sample sequence
sample_seq <- "MKTIIALSYIFCLVFADYKDDDDK"
aa_composition(sample_seq)

# We get the following
# A          R          N          D          C 
# 0.08333333 0.00000000 0.00000000 0.20833333 0.04166667 
# Q          E          G          H          I 
# 0.00000000 0.00000000 0.00000000 0.00000000 0.12500000 
# L          K          M          F          P 
# 0.08333333 0.12500000 0.04166667 0.08333333 0.00000000 
# S          T          W          Y          V 
# 0.04166667 0.04166667 0.00000000 0.08333333 0.04166667

# Very good, now we can use it against the sequences we have from UniProt

# We calculate amino acid composition for each protein and store as a list-column
df_seq_composition <- df_seq %>%
  mutate(aa_comp = map(Sequence, aa_composition))

# Then unnest the list-column to long format:
df_aa <- df_seq_composition %>%
  select(Entry, match_status, aa_comp) %>%
  unnest_wider(aa_comp) %>% 
  pivot_longer(cols = all_of(standard_aa), names_to = "AA", values_to = "fraction")

# Finally view a sample of the resulting data
head(df_aa)

# Entry      match_status        X AA    fraction
# <chr>      <chr>           <dbl> <chr>    <dbl>
# 1 A0A3B3IU09 Autoimmune Only    NA A       0.0722
# 2 A0A3B3IU09 Autoimmune Only    NA R       0.0425
# 3 A0A3B3IU09 Autoimmune Only    NA N       0.0382
# 4 A0A3B3IU09 Autoimmune Only    NA D       0.0510
# 5 A0A3B3IU09 Autoimmune Only    NA C       0.0170
# 6 A0A3B3IU09 Autoimmune Only    NA Q       0.0361

# Good to see it going through, now we can summarize by group

aa_summary <- df_aa %>%
  group_by(match_status, AA) %>%
  summarise(avg_fraction = mean(fraction, na.rm = TRUE), .groups = "drop")

# View the summary
aa_summary

# match_status    AA    avg_fraction
# <chr>           <chr>        <dbl>
# 1 Autoimmune Only A           0.0712
# 2 Autoimmune Only C           0.0227
# 3 Autoimmune Only D           0.0443
# 4 Autoimmune Only E           0.0590
# 5 Autoimmune Only F           0.0438
# 6 Autoimmune Only G           0.0718
# 7 Autoimmune Only H           0.0253
# 8 Autoimmune Only I           0.0474
# 9 Autoimmune Only K           0.0539
# 10 Autoimmune Only L           0.103 

# To get a better look at it, let us visualize the information

frequency_first_plot <- ggplot(aa_summary, aes(x = AA, y = avg_fraction, fill = match_status)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
       x = "Amino Acid",
       y = "Average Fraction",
       fill = NULL) +
  theme_minimal() + 
  theme(
    legend.text     = element_text(size = 16),
    legend.title    = element_text(size = 16)
    )

ggsave(
  filename = "results/first_frequency_plot.png",
       plot = frequency_first_plot,
       width = 12,
       height = 8,
       dpi = 400
       )


# We here get a plot comparing the composition, but let us run a statistical test

results_amino_acid_p <- df_aa %>%
  group_by(AA) %>%
  wilcox_test(fraction ~ match_status) %>%
  adjust_pvalue(method = "BH") %>%   # adjust for multiple testing (Benjamini-Hockberg)
  arrange(p)

head(results_amino_acid_p, 10)

# AA    .y.      group1  group2    n1    n2 statistic       p
# <chr> <chr>    <chr>   <chr>  <int> <int>     <dbl>   <dbl>
# 1 C     fraction Autoim… Healt…    69    47     2166. 0.00219
# 2 L     fraction Autoim… Healt…    69    47     2155  0.00272
# 3 S     fraction Autoim… Healt…    69    47     2146. 0.00318
# 4 G     fraction Autoim… Healt…    69    47     1117  0.00459
# 5 W     fraction Autoim… Healt…    69    47     2114. 0.00567
# 6 K     fraction Autoim… Healt…    69    47     1162. 0.00992
# 7 F     fraction Autoim… Healt…    69    47     2036. 0.02   
# 8 H     fraction Autoim… Healt…    69    47     2022  0.0245 
# 9 E     fraction Autoim… Healt…    69    47     1367  0.153  
# 10 A     fraction Autoim… Healt…    69    47     1374. 0.164 

## We do see some results that could be meaningful, but after p.adj they quickly change due to taking into acount the multiple testing. There is therefore little significance

#___________________________________________________________________________________________________

# Now we move on to the gene enrichment analysis using the gprofiler2 package here in R
# The reason for not making use of the subcellular location, function and pathway columns in the UniProt metadata, is due to many missing values, and untidy data it contains
# Here using this package we can make use of the gene names from UniProt to run the enrichment which also provides us with p-values calculated hypergeometically

# First we extract unique gene symbols for each group.

genes_autoimmune_only <- full_data_subset %>%
  filter(match_status == "Autoimmune Only") %>%
  pull(`Gene Names`) %>% 
  unique()

genes_healthy_autoimmune <- full_data_subset %>%
  filter(match_status == "Healthy & Autoimmune") %>%
  pull(`Gene Names`) %>% 
  unique()

# Print number of unique genes for each group
cat("Number of unique genes in Autoimmune Only:", length(genes_autoimmune_only), "\n")
cat("Number of unique genes in Healthy & Autoimmune:", length(genes_healthy_autoimmune), "\n")

# Number of unique genes in Autoimmune Only: 69 
# Number of unique genes in Healthy & Autoimmune: 47 

# Now, we perform the GO Enrichment Analysis using gprofiler2's gost function, starting off with the autoimmune only group.

gostres_autoimmune_only <- gost(
  query = genes_autoimmune_only,
  organism = "hsapiens",       # human
  domain_scope = "annotated",  # use only annotated genes
  sources = c("GO:BP"),        # restrict to GO Biological Process
  correction_method = "fdr",   # multiple testing correction method
  user_threshold = 0.05        # p-value cutoff
)

# Then the healthy & autoimmune group

gostres_healthy_autoimmune <- gost(
  query = genes_healthy_autoimmune,
  organism = "hsapiens",
  domain_scope = "annotated",
  sources = c("GO:BP"),
  correction_method = "fdr",
  user_threshold = 0.05
)

# Now let us check the results

# We filter significant terms and take top 10
res_auto <- gostres_autoimmune_only$result %>%
  filter(p_value < 0.05) %>%
  arrange(p_value) %>%
  slice_head(n = 10)

# and Plot BP with ggplot2
gprofiler2_auto_BP <- ggplot(res_auto, aes(x = reorder(term_name, -log10(p_value)), y = -log10(p_value))) +
  geom_col(fill = "#F8766D") +
  coord_flip() +
  labs(
    x = "GO BP",
    y = NULL
  ) +
  theme_minimal()

# Repeat for Healthy

res_healthy <- gostres_healthy_autoimmune$result %>%
  filter(p_value < 0.05) %>%
  arrange(p_value) %>%
  slice_head(n = 10)

gprofiler2_health_BP <- ggplot(res_healthy, aes(x = reorder(term_name, -log10(p_value)), y = -log10(p_value))) +
  geom_col(fill = "#00BFC4") +
  coord_flip() +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_minimal()

# That was for Biological Processes, now we run for Cellular Components

# For the autoimmune only
gostres_cc_autoimmune_only <- gost(
  query = genes_autoimmune_only,
  organism = "hsapiens",         # Human
  sources = c("GO:CC"),           # Limit to Cellular Component terms
  correction_method = "fdr",      # Use FDR correction
  user_threshold = 0.05           # p-value cutoff
)

# For the Healthy & Autoimmune group:
gostres_cc_healthy_autoimmune <- gost(
  query = genes_healthy_autoimmune,
  organism = "hsapiens",
  sources = c("GO:CC"),
  correction_method = "fdr",
  user_threshold = 0.05
)

# The Visualization (They will all be included in one plot later)

# Autoimmune Only group:
df_cc_autoimmune <- gostres_cc_autoimmune_only$result %>% 
  filter(p_value < 0.05) %>% 
  arrange(p_value) %>% 
  slice_head(n = 10)

gprofiler2_auto_CC <- ggplot(df_cc_autoimmune, aes(x = reorder(term_name, -log10(p_value)), y = -log10(p_value))) +
  geom_col(fill = "#F8766D") +
  coord_flip() +
  labs(
    x = "GO CC",
    y = NULL
  ) +
  theme_minimal()

# Healthy & Autoimmune group:
df_cc_healthy <- gostres_cc_healthy_autoimmune$result %>% 
  filter(p_value < 0.05) %>% 
  arrange(p_value) %>% 
  slice_head(n = 10)

gprofiler2_health_CC <- ggplot(df_cc_healthy, aes(x = reorder(term_name, -log10(p_value)), y = -log10(p_value))) +
  geom_col(fill = "#00BFC4") +
  coord_flip() +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_minimal()

# Now for the molecular function

gostres_mf_autoimmune_only <- gost(
  query = genes_autoimmune_only,
  organism = "hsapiens",
  sources = c("GO:MF"),
  correction_method = "fdr",
  user_threshold = 0.05
)

gostres_mf_healthy_autoimmune <- gost(
  query = genes_healthy_autoimmune,
  organism = "hsapiens",
  sources = c("GO:MF"),
  correction_method = "fdr",
  user_threshold = 0.05
)

# Autoimmune Only group:

df_mf_autoimmune <- gostres_mf_autoimmune_only$result %>% 
  filter(p_value < 0.05) %>% 
  arrange(p_value) %>% 
  slice_head(n = 10)

gprofiler2_auto_MF <- ggplot(df_mf_autoimmune, aes(x = reorder(term_name, -log10(p_value)), y = -log10(p_value))) +
  geom_col(fill = "#F8766D") +
  coord_flip() +
  labs(
    x = "GO MF",
    y = "-log10(p-value)"
  ) +
  theme_minimal()

df_mf_healthy <- gostres_mf_healthy_autoimmune$result %>% 
  filter(p_value < 0.05) %>% 
  arrange(p_value) %>% 
  slice_head(n = 10)

gprofiler2_health_MF <- ggplot(df_mf_healthy, aes(x = reorder(term_name, -log10(p_value)), y = -log10(p_value))) +
  geom_col(fill = "#00BFC4") +
  coord_flip() +
  labs(
    x = NULL,
    y = "-log10(p-value)"
  ) +
  theme_minimal()

# We take the six figures and patchwork it into one plot
combined_go_gprofiler <- 
  (gprofiler2_auto_BP  | gprofiler2_health_BP)  /
  (gprofiler2_auto_CC  | gprofiler2_health_CC)  /
  (gprofiler2_auto_MF  | gprofiler2_health_MF)  +
  plot_layout(guides = "collect") &                      
  theme_minimal(base_size = 18) &              
  theme(
    legend.position = "none",                  
    plot.margin     = margin(5, 5, 5, 5)       
  )
# save to file
ggsave(
  filename = "results/combined_gprofiler_enrichment.png",
  plot     = combined_go_gprofiler,
  width    = 20,    # inches
  height   = 20,    # inches
  dpi      = 400
)

# Laslty we have Pathways

gostres_pathway_autoimmune <- gost(
  query = genes_autoimmune_only,
  organism = "hsapiens",         # Human
  sources = c("REAC"),           # Reactome pathways
  correction_method = "fdr",     # FDR correction for multiple testing
  user_threshold = 0.05          # p-value cutoff
)

# Check the result
df_pathway_autoimmune <- gostres_pathway_autoimmune$result %>%
  filter(p_value < 0.05) %>%  # Keep significant terms
  arrange(p_value) %>%
  slice_head(n = 10)        # Select top 10 enriched pathways (adjust as needed)

# We here reate a dot plot using ggplot2:
gprofiler2_auto_Path <- ggplot(df_pathway_autoimmune, aes(
  x = -log10(p_value),
  y = reorder(term_name, -log10(p_value)),
  size = intersection_size,      # number of genes that overlap the pathway
)) +
  geom_point() +
  scale_color_viridis_c() +
  labs(
    x = "-log10(p-value)",
    y = "Pathway",
    size = "Intersection size (Auto)",
  ) +
  theme_minimal()

gostres_pathway_healthy <- gost(
  query = genes_healthy_autoimmune,
  organism = "hsapiens",         # Human
  sources = c("REAC"),           # Reactome pathways
  correction_method = "fdr",     # FDR correction
  user_threshold = 0.05          # p-value cutoff
)

df_pathway_healthy <- gostres_pathway_healthy$result %>% 
  filter(p_value < 0.05) %>% 
  arrange(p_value) %>% 
  slice_head(n = 10)

gprofiler2_health_Path <- ggplot(df_pathway_healthy, aes(
  x = -log10(p_value),
  y = reorder(term_name, -log10(p_value)),
  size = intersection_size,
)) +
  geom_point() +
  scale_color_viridis_c() +
  labs(
    x = "-log10(p-value)",
    y = NULL,
    size = "Intersection Size (Healthy & Auto)",
  ) +
  theme_minimal()

# Now we wish to combine all of these into one plot, this is done using patchwork

# 1. Combine the eight plots with patchwork
combined_pathway_enrichment_gprofiler2 <- 
  (gprofiler2_auto_Path   | gprofiler2_health_Path) +  
  plot_layout(guides = "collect") theme_minimal(base_size = 28) &            
  theme(
    axis.title    = element_text(size = 24),  
    axis.text     = element_text(size = 24),  
    legend.title  = element_text(size = 24),  
    legend.text   = element_text(size = 24),  
    plot.margin   = margin(5, 5, 5, 5)        
  )

ggsave(
  filename = "results/combined_pathway_enrichment_gprofiler2.png",
  plot     = combined_pathway_enrichment_gprofiler2,
  width    = 16,    # inches
  height   = 8,    # inches
  dpi      = 400
)

#_______________________________________________________________________________________________

## Besides the UniProt metadata, we also took Fasta files of the two group, which allowed us to run them through the various systems of DTU Healthtech
# The first one was SignalP 6
# Here we downloaded two prediction summaries in txt format and ran then through to test is there was any difference between the two groups

sig_matched <- read_delim(
  "matched_prediction_results.txt",
  delim      = "\t",
  comment    = "#",
  col_names  = c("Entry","Description","prob_OTHER","prob_SP","CSinfo"),
  trim_ws    = TRUE
) %>%
  mutate(match_status = "Healthy & Autoimmune")

sig_unmatched <- read_delim(
  "unmatched_prediction_results.txt",
  delim      = "\t",
  comment    = "#",
  col_names  = c("Entry","Description","prob_OTHER","prob_SP","CSinfo"),
  trim_ws    = TRUE
) %>%
  mutate(match_status = "Autoimmune Only")

# We then chose to combine them
signalp <- bind_rows(sig_matched, sig_unmatched)

# It was then pivoted to a long format
signalp_long <- signalp %>%
  select(Entry, match_status, prob_OTHER, prob_SP) %>%
  pivot_longer(
    cols      = c(prob_OTHER, prob_SP),
    names_to  = "Compartment",
    values_to = "Probability"
  ) %>%
  mutate(
    Compartment = recode(
      Compartment,
      prob_OTHER = "Other",
      prob_SP    = "Signal peptide"
    )
  )

# From this we ran a per compartment wilcoxon test to see if there were any interesting results
pvals <- signalp_long %>%
  group_by(Compartment) %>%
  wilcox_test(Probability ~ match_status) %>%
  adjust_pvalue(method = "BH") %>%
  select(Compartment, p.adj) %>%
  mutate(p.label = paste0("p=", signif(p.adj, 2)))

# Afterwich we plottet the results

signalP_plot <- ggplot(signalp_long, aes(Compartment, Probability, color = match_status)) +
  geom_jitter(position = position_jitter(width = .2), alpha = .8, size = 1.5) +
  geom_boxplot(aes(fill = match_status), outlier.shape = NA, alpha = .5, width = .6) +
  geom_text(
    data = pvals,
    aes(x = Compartment, y = 1.05, label = p.label),
    inherit.aes = FALSE, size = 3
  ) +
  scale_y_continuous(expand = expansion(c(0, .15))) +
  # …and the rest of your styling…
  labs(x = NULL, y = "Probability", fill = NULL, color = NULL) +
  theme_minimal() +
  theme(
    legend.text = element_text(20),
    legend.title = element_text(24)
  )

ggsave(
  filename = "results/signalP_plot.png",
  plot     = signalP_plot,
  width    = 10,    # inches
  height   = 5,    # inches
  dpi      = 400
)

# The plot and p-values showed nothing significant, therefore we moved on to NetSurfP 3.0

#_________________________________________________________________________________________________

# Using the Fasta file once again, we this time made use of the service NetSurfP 3.0, which gave us two csv result files. 

ns_matched <- read_csv("matched_netsurfp_results.csv", show_col_types = FALSE) %>%
  rename_with(str_trim) %>%
  mutate(
    rsa         = as.numeric(rsa),
    disorder    = as.numeric(disorder),
    q3          = str_trim(q3),
    match_status = "Healthy & Autoimmune"
  )

ns_unmatched <- read_csv("unmatched_netsurfp_results.csv", show_col_types = FALSE) %>%
  rename_with(str_trim) %>%
  mutate(
    rsa         = as.numeric(rsa),
    disorder    = as.numeric(disorder),
    q3          = str_trim(q3),
    match_status = "Autoimmune Only"
  )

# We combined and cleaned the entry id's
ns_all <- bind_rows(ns_matched, ns_unmatched) %>%
  mutate(Entry = str_extract(id, "^>?(\\S+)"))

# In the beginning, we wished to take a look at all of them, and therefore did a per-protein summary
protein_summary <- ns_all %>%
  group_by(Entry, match_status) %>%
  summarise(
    mean_RSA      = mean(rsa,       na.rm = TRUE),
    prop_exposed  = mean(rsa > 0.5, na.rm = TRUE),
    mean_disorder = mean(disorder,  na.rm = TRUE),
    frac_H        = mean(q3 == "H", na.rm = TRUE),
    frac_E        = mean(q3 == "E", na.rm = TRUE),
    frac_C        = mean(q3 == "C", na.rm = TRUE),
    .groups = "drop"
  )

# This was then pivoted
long_summary <- protein_summary %>%
  pivot_longer(
    cols      = mean_RSA:frac_C,
    names_to  = "Metric",
    values_to = "Value"
  )

# We wished to only keep metrics with ≥2 observations in each group
valid_metrics <- long_summary %>%
  group_by(Metric, match_status) %>%
  summarise(n = sum(!is.na(Value)), .groups = "drop") %>%
  pivot_wider(names_from = match_status, values_from = n) %>%
  filter(`Autoimmune Only` >= 2, `Healthy & Autoimmune` >= 2) %>%
  pull(Metric)

long_summary <- filter(long_summary, Metric %in% valid_metrics)

# The we ran the Wilcoxon test
pvals <- long_summary %>%
  group_by(Metric) %>%
  wilcox_test(Value ~ match_status) %>%
  adjust_pvalue(method = "BH") %>%
  select(Metric, p.adj) %>%
  mutate(p.label = paste0("p=", signif(p.adj, 2)))

# We then plottet the results
primary_summary_netsurfp <- ggplot(long_summary, aes(match_status, Value, fill = match_status)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.3) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
  facet_wrap(~ Metric, scales = "free_y", ncol = 2) +
  geom_text(
    data        = pvals,
    aes(x = 2, y = Inf, label = p.label),
    inherit.aes = FALSE,
    vjust       = 1.5,
    size        = 5
  ) +
  labs(
    x     = "Group",
    y     = NULL,
    title = "NetSurfP 3: Structural summaries by group"
  ) +
  theme_minimal() +
  theme(
    legend.position   = "none",
    axis.text.x       = element_text(angle = 45, hjust = 1),
    strip.text        = element_text(face = "bold")
  )

ggsave("results/primary_summary_netsurfp.png", 
       primary_summary_netsurfp,
       width = 12, height = 10, dpi = 400)

# This did show us where there were some interesting information to go on, but there were also elements we missed out on when doing this.
# The first was, that this way was not a directly viable way for us to look at the disorder. 
# We therefore choose to look at two other ways to test the disorder

ns_matched   <- read_csv("matched_netsurfp_results.csv",   show_col_types = FALSE) %>%
  mutate(match_status = "Healthy & Autoimmune")
ns_unmatched <- read_csv("unmatched_netsurfp_results.csv", show_col_types = FALSE) %>%
  mutate(match_status = "Autoimmune Only")

ns_all <- bind_rows(ns_matched, ns_unmatched) %>%
  mutate(Entry = str_extract(id, "^\\S+"))

# First summarise each protein’s disorder runs into one row per protein
region_summary <- ns_all %>%
  group_by(Entry, match_status) %>%
  summarise(
    runs = list( rle(disorder > 0.5)$lengths[ rle(disorder > 0.5)$values ] ),
    max_run_length   = ifelse(length(runs[[1]])>0, max(runs[[1]]), 0),
    total_run_length = sum(runs[[1]]),
    .groups = "drop"
  )

# We naturally did some statistical tests
p_max   <- wilcox_test(region_summary, max_run_length   ~ match_status) %>% pull(p)
p_total <- wilcox_test(region_summary, total_run_length ~ match_status) %>% pull(p)

# Now the first thing we wished to take a look at was the longest disordered segment per protein
g1 <- ggplot(region_summary, aes(x = match_status, y = max_run_length, fill = match_status)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.8) +
  annotate(
    "text", x = 1.5, y = max(region_summary$max_run_length)*0.95,
    label = paste0("p = ", signif(p_max, 2)),
    size = 5
  ) +
  scale_fill_manual(values = c(
    "Autoimmune Only"        = "#F8766D",
    "Healthy & Autoimmune"   = "#00BFC4"
  )) +
  labs(
    x     = NULL,
    y     = NULL
  ) +
  theme_minimal(base_size = 18) +
  theme(legend.position = "none")

print(g1)


# This was the second thing, the total disorder across all runs per protein
g2 <- ggplot(region_summary, aes(x = match_status, y = total_run_length, fill = match_status)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.8) +
  annotate(
    "text", x = 1.5, y = max(region_summary$total_run_length)*0.95,
    label = paste0("p = ", signif(p_total, 2)),
    size = 5
  ) +
  scale_fill_manual(values = c(
    "Autoimmune Only"        = "#F8766D",
    "Healthy & Autoimmune"   = "#00BFC4"
  )) +
  labs(
    x     = NULL,
    y     = NULL
  ) +
  theme_minimal(base_size = 18) +
  theme(legend.position = "none")

print(g2)

combined_disorder_fig <- g1 + g2 + 
  plot_layout(ncol = 2) +                      
  plot_annotation(tag_levels = "A",
    theme = theme_minimal() & theme(plot.tag = element_text(face="bold",
                                                            size = 12)
                                    ))

ggsave(
  filename = "results/combined_disorder_fig.png",
  plot     = combined_disorder_fig,
  width    = 12,    # adjust as needed (in inches)
  height   = 6,
  dpi      = 400
)

# While the figure did show somewhat of a difference, nothing statistically important came from it. 
# We therefore moved on to focus on Alpha helices and beta sheets
# The reason for this being that there has to be at least four H's in a row for it to actually be a helix.
# In the same sense, there needs to be 3 E's in a row for it to be considered a sheet

ns_matched   <- read_csv("matched_netsurfp_results.csv")   %>% mutate(match_status = "Healthy & Auto")
ns_unmatched <- read_csv("unmatched_netsurfp_results.csv") %>% mutate(match_status = "Autoimmune Only")
ns_all       <- bind_rows(ns_matched, ns_unmatched)

# First we computed the raw count per protein
helix_sheet_counts <- ns_all %>%
  rename(Entry = id) %>%            # here is the per-residue NetSurfP ID
  select(Entry, match_status, q3) %>% # q3 is the secondary‐structure cell
  group_by(Entry, match_status) %>%
  summarize(
    Length   = n(),                 # the total residues will equal the protein length
    count_H4 = {                     # and we make sure it runs with H ≥ 4
      rle_h <- rle(q3 == "H")
      sum(rle_h$values & rle_h$lengths >= 4)
    },
    count_E3 = {                     # here it runs with E ≥ 3
      rle_e <- rle(q3 == "E")
      sum(rle_e$values & rle_e$lengths >= 3)
    },
    .groups="drop"
  ) %>%
  # We must remember to normalize to runs per 100 aa, so that all lengths are comparable:
  mutate(
    H4_per_100 = count_H4 / Length * 100,
    E3_per_100 = count_E3 / Length * 100
  )

# we do a quick sanity check:
helix_sheet_counts %>% 
  summarize(
    proteins = n_distinct(Entry),
    avg_len = mean(Length),
    avg_H4 = mean(count_H4),
    avg_H4_norm = mean(H4_per_100)
  ) %>% print()

# proteins avg_len avg_H4 avg_H4_norm
# <int>   <dbl>  <dbl>       <dbl>
#  116    652.   10.1        1.78

# With the standardized results we now run wilcoxon tests
tests <- helix_sheet_counts %>%
  pivot_longer(
    cols      = c(count_H4, H4_per_100, count_E3, E3_per_100),
    names_to  = "metric",
    values_to = "value"
  ) %>%
  group_by(metric) %>%
  summarise(
    p = wilcox.test(value ~ match_status)$p.value,
    .groups = "drop"
  )

print(tests)

# metric           p
# <chr>        <dbl>
# 1 E3_per_100 0.628  
# 2 H4_per_100 0.0787 
# 3 count_E3   0.699  
# 4 count_H4   0.00199

# We could here see that after normalizing the data, there was actually a significant p-value
# Yet this p-value would have not been usuable due to the diffence in amount of proteins and difference in length

# Now we finally come to the plotting part
p_H4_norm <- tests %>% filter(metric=="H4_per_100") %>% pull(p)
p_E3_norm <- tests %>% filter(metric=="E3_per_100") %>% pull(p)
group_cols <- c("Healthy & Auto" = "#00BFC4", "Autoimmune Only" = "#F8766D")

p1_norm <- ggplot(helix_sheet_counts, aes(
  x = match_status,
  y = H4_per_100,
  fill = match_status
)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.8) +
  scale_fill_manual(
    name   = NULL,
    values = group_cols
  ) +
  annotate(
    "text", x = 1.5, 
    y = max(helix_sheet_counts$H4_per_100, na.rm=TRUE) * 1.05,
    label = paste0("p = ", signif(p_H4_norm,3)),
    size  = 4
  ) +
  labs(
    x     = NULL,
    y     = NULL
  ) +
  theme_minimal(base_size = 18) +
  theme(
    legend.position = "none",
    plot.title      = element_text(face="bold", hjust=0.5)
  )

p2_norm <- ggplot(helix_sheet_counts, aes(
  x = match_status,
  y = E3_per_100,
  fill = match_status
)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.8) +
  scale_fill_manual(
    name   = NULL,
    values = group_cols
  ) +
  annotate(
    "text", x = 1.5,
    y = max(helix_sheet_counts$E3_per_100, na.rm=TRUE) * 1.05,
    label = paste0("p = ", signif(p_E3_norm,3)),
    size  = 4
  ) +
  labs(
    x     = NULL,
    y     = NULL
  ) +
  theme_minimal(base_size = 18) +
  theme(
    legend.position = "none",
    plot.title      = element_text(face="bold", hjust=0.5)
  )

combined_norm <- p1_norm + p2_norm + 
  plot_layout(ncol = 2, guides = "collect") & 
  theme(legend.position = "bottom")

final_fig <- combined_norm + 
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face="bold"))

print(final_fig)

ggsave(
  "results/helix_sheet_normalized.png",
  final_fig,
  width  = 16, height = 8, units = "in", dpi = 400)


# Now for the last part we wished to check from the NetSurfP data, we actually saw something interesting in regards to prop_exposed, along with mean_RSA
# Now, to explain:
# The mean RSA is the average relative surface accessibility across all its residues
# % Exposed or prop_exposed is the fraction of residues with RSA > 0.5, but this will be changed, since the recommended threshold is 0.25 according to NetSurfP

# We therefore call this part Surface Accessibility Analysis

matched   <- read_csv("matched_netsurfp_results.csv",   show_col_types = FALSE)
unmatched <- read_csv("unmatched_netsurfp_results.csv", show_col_types = FALSE)

# We compute the metrics
compute_metrics <- function(df) {
  df %>%
    group_by(id) %>%
    summarise(
      avg_rsa      = mean(rsa, na.rm = TRUE),
      frac_exposed = mean(rsa > 0.25, na.rm = TRUE),
      .groups      = "drop"
    )
}

metrics_matched   <- compute_metrics(matched)   %>% mutate(group = "Healthy & Autoimmune")
metrics_unmatched <- compute_metrics(unmatched) %>% mutate(group = "Autoimmune Only")
metrics_all       <- bind_rows(metrics_matched, metrics_unmatched)

# Run a couple of Wilcoxon tests
test_avg  <- wilcox.test(avg_rsa      ~ group, data = metrics_all)
test_frac <- wilcox.test(frac_exposed ~ group, data = metrics_all)

pv_avg  <- signif(test_avg$p.value, 3)
pv_frac <- signif(test_frac$p.value,3)

# Here we actually see some positive results
# pv_avg   0.00272
# pv_frac  0.012

# We set a limit, since it goes from 0-1
x_limits <- c(0, 1)
base_theme <- theme_minimal(base_size=14)


my_colors <- c(
  "Healthy & Autoimmune"   = "#00BFC4",  # ggplot’s default “Unmatched” color
  "Autoimmune Only" = "#F8766D"   # ggplot’s default “Matched” color
)

# We look at a histogram of the average rsa, with density scaled
p_avg_hist <- ggplot(metrics_all, aes(x = avg_rsa, fill = group)) +
  geom_histogram(aes(y = ..density..), 
                 position = "identity", 
                 alpha = 0.5, 
                 bins = 30) +
  scale_x_continuous(limits = x_limits) +
  scale_fill_manual(values = my_colors) +
  labs(
    subtitle = paste0("Wilcoxon p = ", pv_avg),
    x        = NULL,
    y        = "Density",
    fill     = NULL
  ) +
  theme_minimal(base_size = 18)

# And one for the frac_exposed
p_frac_hist <- ggplot(metrics_all, aes(x = frac_exposed, fill = group)) +
  geom_histogram(aes(y = ..density..), 
                 position = "identity", 
                 alpha = 0.5, 
                 bins = 30) +
  scale_x_continuous(limits = x_limits) +
  scale_fill_manual(values = my_colors) +
  labs(
    subtitle = paste0("Wilcoxon p = ", pv_frac),
    x        = NULL,
    y        = NULL,
    fill     = NULL
  ) +
  theme_minimal(base_size = 18)

combined_rsa_plots_hist <- (p_avg_hist | p_frac_hist) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave(
  filename = "results/combined_rsa_hist_plots.png",
  plot = combined_rsa_plots_hist,
  width = 12,
  height = 6,
  dpi = 400
)


# We choose to also make a smooth density plot of average RSA
p_avg <- ggplot(metrics_all, aes(x=avg_rsa, fill=group)) +
  geom_density(alpha=0.5, adjust=1.2) +
  scale_x_continuous(limits=x_limits) +
  labs(
    title    = "Average RSA per Protein",
    subtitle = paste0("Wilcoxon p = ", pv_avg),
    x        = "Mean RSA",
    y        = "Density",
    fill     = NULL
  ) +
  base_theme

# and a density of fraction of exposed residues
p_frac <- ggplot(metrics_all, aes(x=frac_exposed, fill=group)) +
  geom_density(alpha=0.5, adjust=1.2) +
  scale_x_continuous(limits=x_limits) +
  labs(
    title    = "Fraction of Exposed Residues (RSA>0.25)",
    subtitle = paste0("Wilcoxon p = ", pv_frac),
    x        = "Proportion Exposed",
    y        = NULL,
    fill     = NULL
  ) +
  base_theme

# We once again combine using patchwork
combined_rsa_plots <- (p_avg | p_frac) +
  plot_layout(guides = "collect") +
  plot_annotation(title = "RSA Metrics: Matched vs Unmatched") &
  theme(
    legend.position = "bottom"
  )

ggsave(
  filename = "results/combined_rsa_plots.png",
  plot     = combined_rsa_plots,
  width    = 12,   
  height   = 6,
  dpi      = 300
)

# We wish to see another version of it, with the colors switched around and one where the density hasn't been smoothed out.




# Here we actually saw some significant results, and this is the end of the NetSurfP testing

#___________________________________________________________________________________________


# The next part we moved onto once again required the fasta files, which we ran through DeepLoc-2.1. 
# This was for a better prediction of the localization. We downloaded the csv files, which contain the Winning location and vectors from 0-1
# We started off by loading the data and tagging them to their own groups

matched_dl <- read_csv("matched_proteins_deeploc_results.csv", show_col_types = FALSE) %>%
  mutate(Group = "Healthy & Autoimmune")  # labels the matched set

unmatched_dl <- read_csv("unmatched_proteins_deeploc_results.csv", show_col_types = FALSE) %>%
  mutate(Group = "Autoimmune Only")       # labels the unmatched set

# We then merged them
all_dl <- bind_rows(matched_dl, unmatched_dl)

# And then looked for the "Winner" location
loc_counts <- all_dl %>%
  count(Group, Localizations, name = "ProteinCount") %>%
  arrange(Group, desc(ProteinCount))

# A small print
cat("\n=== DeepLoc Localization Counts ===\n")
print(loc_counts)

# We then plottet them side by side in a bar plot
win_deeploc <- ggplot(loc_counts, aes(
  x    = Localizations,            # The compartment on x-axis
  y    = ProteinCount,             # The count on y
  fill = Group                     # And then it is colored by group
)) +
  geom_col(
    position = position_dodge(width = 0.8),  # We dodge bars side-by-side
    width    = 0.7
  ) +
  labs(
    title    = "DeepLoc-2.1 Predicted Localization by Group",
    subtitle = "Counts of proteins in each compartment",
    x        = NULL,
    y        = "Number of Proteins"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1),  # and lastly tilt labels
  ) + guides(fill = FALSE)

ggsave(
  filename = "results/win_deeploc.png",
  plot     = win_deeploc,
  width    = 12,   
  height   = 6,
  dpi      = 300
)

# We also did a run on signal type, which we wish to insert into the plot as it can also show us interesting values

# We break out the semicolon-separated signal column
signals_tbl <- all_dl %>%
  select(Protein_ID, Group, Signals) %>%
  separate_rows(Signals, sep = ";") %>%
  mutate(Signals = str_trim(Signals)) %>%
  distinct(Protein_ID, Group, Signals)

# Do a fisher test per signal with threshold
signal_results <- signals_tbl %>%
  mutate(present = TRUE) %>%
  pivot_wider(names_from = Signals, values_from = present, values_fill = FALSE) %>%
  pivot_longer(-c(Protein_ID, Group), names_to = "Signal", values_to = "Present") %>%
  group_by(Signal) %>%
  summarize(
    ct    = list(table(Group, Present)),
    p.adj = p.adjust(fisher.test(table(Group, Present))$p.value, method = "BH"),
    .groups = "drop"
  ) %>%
  filter(p.adj < 0.05)

# Then we build a plotting table of counts and percentage
viz_df <- signal_results %>%
  mutate(
    auto_yes    = map_int(ct, ~ .x["Autoimmune Only",    "TRUE"]),
    healthy_yes = map_int(ct, ~ .x["Healthy & Autoimmune","TRUE"]),
    total_auto    = map_int(ct, ~ sum(.x["Autoimmune Only",])),
    total_healthy = map_int(ct, ~ sum(.x["Healthy & Autoimmune",]))
  ) %>%
  select(Signal, auto_yes, healthy_yes, total_auto, total_healthy) %>%
  pivot_longer(c(auto_yes, healthy_yes), names_to = "Group", values_to = "n_yes") %>%
  mutate(
    Group = recode(Group,
                   auto_yes    = "Autoimmune Only",
                   healthy_yes = "Healthy & Autoimmune"),
    total = if_else(Group == "Autoimmune Only", total_auto, total_healthy),
    pct   = n_yes / total * 100,
    label = paste0(n_yes, "/", total)
  )

# Now we can plot it
win_signals <- ggplot(viz_df, aes(x = Signal, y = pct, fill = Group)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = label),
            position = position_dodge(width = 0.8),
            vjust = -0.2, size = 3) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  labs(
    x        = NULL,
    y        = "Percent of Proteins with Signal",
    fill     = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

ggsave(
  filename = "results/winning_signals_loc.png",
  plot = win_signals,
  width = 12,
  height = 6,
  dpi = 400
)

# We then combine the two figures
combined_win_signal_loc <- 
  win_deeploc / win_signals +               # stack them
  plot_layout(heights = c(1, 1), guides = "collect") +
  plot_annotation(
    title = "DeepLoc-2.1 Predictions: Subcellular Localizations & Signal Types"
  ) &
  theme(legend.position = "bottom")

ggsave(
  filename = "results/combined_win_signal_loc.png",
  plot     = combined_win_signal_loc,
  width    = 12,
  height   = 12,
  dpi      = 300
)

# Now, this already looked very interesting, and we could see clear differences, but we can expand it. 
# In the csv files, every location for every protein has a vector value, ranging from 0-1, with the highest one winning.
# But these values are still here, and can give us much information. 

# So let us reload the documents, and try this again
matched <- read_csv("matched_proteins_deeploc_results.csv", show_col_types=FALSE) %>%
  mutate(Group = "Healthy & Auto")
unmatched <- read_csv("unmatched_proteins_deeploc_results.csv", show_col_types=FALSE) %>%
  mutate(Group = "Autoimmune Only")

deeploc <- bind_rows(matched, unmatched)

# We remember to pivot
loc_cols <- c(
  "Cytoplasm", "Nucleus", "Extracellular", "Cell membrane",
  "Mitochondrion", "Endoplasmic reticulum",
  "Lysosome/Vacuole", "Golgi apparatus", "Peroxisome"
)

deeploc_long <- deeploc %>%
  pivot_longer(
    cols      = all_of(loc_cols),
    names_to  = "Subcellular localization",
    values_to = "Likelihood"
  ) %>%
  mutate(
    `Subcellular localization` = factor(`Subcellular localization`, levels = loc_cols)
  )

# We use Wilcoxon test per location
pvals <- deeploc_long %>%
  group_by(`Subcellular localization`) %>%
  summarise(
    p.value = wilcox.test(Likelihood ~ Group, data = cur_data_all())$p.value
  ) %>%
  ungroup() %>%
  mutate(
    signif = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE            ~ "ns"
    ),
    label = paste0("P = ", signif(p.value, 2), " ", signif)
  )

# Now we move on to the plotting part

p_final_deeploc <- ggplot(deeploc_long, aes(
  x = Likelihood,
  y = `Subcellular localization`,
  colour = Group,
  fill   = Group
)) +
  geom_jitter(position = position_jitter(height = 0.15), alpha = 0.6, size = 1.5) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.3) +
  scale_colour_manual(
    name   = NULL,
    values = c("Autoimmune Only" = "#F8766D", "Healthy & Auto" = "#00BFC4")
  ) +
  scale_fill_manual(
    name   = NULL,
    values = c("Autoimmune Only" = "#F8766D", "Healthy & Auto" = "#00BFC4")
  ) +
  labs(
    x     = "Predicted Likelihood",
    y     = NULL
  ) +
  theme_bw(base_size = 18) +
  theme(
    legend.position    = "right",
    panel.grid.major.y = element_blank()
  )


print(p_final_deeploc)
print(pvals)

ggsave("results/p_final_deeploc.png",
       plot = p_final_deeploc,
       width  = 12,   
       height = 8,
       units  = "in",
       dpi    = 400)

# We here have some interesting results
# `Subcellular localization`   p.value signif label          
# <fct>                          <dbl> <chr>  <chr>          
# 1 Cytoplasm                  0.202     ns     P = 0.2 ns     
# 2 Nucleus                    0.393     ns     P = 0.39 ns    
# 3 Extracellular              0.0248    *      P = 0.025 *    
# 4 Cell membrane              0.0000259 ***    P = 2.6e-05 ***
# 5 Mitochondrion              0.482     ns     P = 0.48 ns    
# 6 Endoplasmic reticulum      0.0122    *      P = 0.012 *    
# 7 Lysosome/Vacuole           0.0289    *      P = 0.029 *    
# 8 Golgi apparatus            0.00937   **     P = 0.0094 **  
# 9 Peroxisome                 0.207     ns     P = 0.21 ns 

## Now that we have worked with this, we will move on.

#______________________________________________________________________________________________


# In the next part we wished to once again do a gene enrichment analysis, this time using enrichR package
# Using this we also in accossiation everything with the diseases. First we loaded the files
# We loaded both the autoimmune and the supplementary table

supp <- read_csv("supplementary_table1.csv", show_col_types=FALSE)

supp_unique <- supp %>%
  filter(!is.na(Source_proteinID)) %>%
  mutate(ProteinID = sub("\\..*$", "", Source_proteinID)) %>%
  distinct(ProteinID)

write_csv(supp_unique, "unique_supplementary_proteins.csv")
message("Distinct proteins in supplementary table: ", nrow(supp_unique))



auto <- read_csv("autoimmune_iedb_pepetides.csv", show_col_types=FALSE)

auto_unique <- auto %>%
  filter(!is.na(`Epitope: Parent Protein Accession`)) %>%
  distinct(`Epitope: Parent Protein Accession`) %>%
  rename(ProteinID = `Epitope: Parent Protein Accession`)

write_csv(auto_unique, "unique_autoimmune_proteins.csv")
message("Distinct proteins in autoimmune IEDB peptides: ", nrow(auto_unique))

# Distinct proteins in supplementary table: 6810
# Distinct proteins in autoimmune IEDB peptides: 116

# Now we focus on from tibble to fasta, setting up a helper. This is quite the long part, that will be reused in further figures

tibble_to_fasta <- function(mytibble, filename='Output.fasta'){
  count <- 0
  for(x in t(mytibble)){
    if(count %% 2 == 1){
      y <- strsplit(x,'(?<=[A-Z]{80})', perl=TRUE)[[1]]
      cat(y, file=filename, sep="\n", append=TRUE)
    } else {
      cat(">", x, "\n", file=filename, sep="", append=TRUE)
    }
    count <- count + 1
  }
}

negative_selection <- function(positives, proteome, ratio, seed=10){
  set.seed(seed)
  proteome_clean <- anti_join(proteome, positives, by='Entry')
  pos_lengths   <- positives$Length
  negatives     <- tibble()
  for(i in seq_len(ratio)){
    for(l in pos_lengths){
      pool <- proteome_clean %>% filter(Length == l)
      if(nrow(pool)==0){
        # jiggle length until we find some
        delta <- 0
        while(nrow(pool)==0){
          delta <- sample(c(-1,0,1),1)
          pool  <- proteome_clean %>% filter(Length == l + delta)
        }
      }
      pick <- pool %>% slice_sample(n=1)
      negatives    <- bind_rows(negatives, pick)
      proteome_clean <- anti_join(proteome_clean, pick, by='Entry')
    }
  }
  negatives
}

aa_composition <- function(string){
  aa_list <- c('A','R','N','D','C','E','Q','G','H',
               'I','L','K','M','F','P','S','T','W','Y','V','X')
  tib <- tibble(sequence=string)
  for(aa in aa_list){
    tib <- tib %>% mutate(!!aa := str_count(sequence, aa))
  }
  tib
}

# We encode a color palette for each disease
distinct_20_blind <- c(
  "#d82772","#ba004c","#ff8e9a","#a20022","#560d00",
  "#880200","#ffa876","#803c00","#664c00","#e1c336",
  "#264d00","#00a847","#76eba5","#91a0ff","#024ebc",
  "#8f77bb","#a051c4","#fdb0ff","#e369d5","#a40071"
)
distinct_colors_21 <- c(
  "alopecia areata"             = "#fc87a3",
  "autoimmune hemolytic anemia" = "#fe0040",
  "autoimmune hepatitis"        = "#da4f00",
  "autoimmune thyroiditis"      = "#f2953f",
  "autoimmune uveitis"          = "#433c01",
  "autoimmune vasculitis"       = "#009030",
  "Behcet's disease"            = "#3bc266",
  "Goodpasture syndrome"        = "#00a177",
  "Graves disease"              = "#00b6fa",
  "multiple sclerosis"          = "#014a78",
  "myasthenia gravis"           = "#004ddd",
  "neuromyelitis optica"        = "#bf95ff",
  "pemphigus"                   = "#741bce",
  "primary biliary cholangitis" = "#502a6b",
  "psoriasis"                   = "#9d74a1",
  "rheumatoid arthritis"        = "#ff5eee",
  "Sjogren's syndrome"          = "#5e2752",
  "systemic lupus erythematosus"= "#ac008a",
  "type 1 diabetes mellitus"    = "#cc0031",
  "vitiligo"                    = "#fdb0ff",
  "Vogt-Koyanagi-Harada disease"= "#e1c336"
)

goplot <- function(data, title="", top=15, pvalue=0.05){
  data %>%
    separate(col=Overlap, into=c("Count","Total"), sep="/", remove=FALSE) %>%
    mutate(
      Count      = as.numeric(Count),
      Gene.Ratio = Count / as.numeric(Total)
    ) %>%
    filter(Adjusted.P.value < pvalue) %>%
    group_by(Data_Type) %>%
    slice_min(order_by=Adjusted.P.value, n=top, with_ties=FALSE) %>%
    ungroup() %>%
    ggplot(aes(
      x     = Gene.Ratio,
      y     = reorder(str_wrap(Term,30), -Adjusted.P.value),
      size  = Count,
      color = Adjusted.P.value
    )) +
    geom_point() +
    facet_wrap(~ Data_Type, scales="free_y", nrow=5) +
    scale_color_viridis(option="D") +
    scale_size(range=c(2,6)) +
    labs(
      title = title,
      x     = "Gene Ratio",
      y     = NULL,
      color = "Adj. P-value",
      size  = "Gene Count"
    ) +
    theme_bw() +
    theme(
      plot.title      = element_text(size=12, face="bold", hjust=0.5),
      strip.text      = element_text(size=10, face="bold"),
      axis.text.y     = element_text(size=6),
      legend.key.size = unit(0.4,"cm")
    )
}

# Here we read the IEDB peptides and clean them

iedb <- read_csv("autoimmune_iedb_pepetides.csv",
                 show_col_types=FALSE)

iedb <- iedb %>%
  rename(DiseaseName = `DiseaseName;;;;;;;;;;;;;`) %>%
  mutate(DiseaseName = str_remove(DiseaseName, ";+$"))

disease_gene <- iedb %>%
  select(DiseaseName, GeneSymbol) %>%
  distinct()

gene_list <- unique(disease_gene$GeneSymbol)
message("▶︎ Testing ", length(gene_list), " unique genes")

# ▶︎ Testing 116 unique genes

# We now run the enrichR

dbs <- c(
  "GO_Molecular_Function_2021",
  "GO_Cellular_Component_2021",
  "GO_Biological_Process_2021",
  "KEGG_2021_Human",
  "InterPro_Domains_2019"
)
enrich_res <- enrichr(gene_list, databases=dbs)
# Uploading data to Enrichr... Done.
# Querying GO_Molecular_Function_2021... Done.
# Querying GO_Cellular_Component_2021... Done.
# Querying GO_Biological_Process_2021... Done.
# Querying KEGG_2021_Human... Done.
# Querying InterPro_Domains_2019... Done.
# Parsing results... Done.

# We now combine and count diseases per term
pval_cut   <- 0.05
combi_data <- tibble()
combi_count<- tibble()

for(db in dbs){
  # count distinct diseases
  count_df <- enrich_res[[db]] %>%
    separate_rows(Genes, sep=";") %>%
    left_join(disease_gene,
              by=c("Genes"="GeneSymbol"),
              relationship="many-to-many") %>%
    filter(!is.na(DiseaseName)) %>%
    group_by(Term) %>%
    summarise(
      DiseaseCount = n_distinct(DiseaseName),
      DiseaseNames = paste(unique(DiseaseName), collapse=", ")
    ) %>%
    ungroup() %>%
    mutate(Data_Type=db)
  
  data_df <- enrich_res[[db]] %>%
    left_join(count_df, by="Term") %>%
    filter(DiseaseCount>1) %>%
    mutate(Data_Type=db)
  
  bar_df <- data_df %>%
    filter(Adjusted.P.value < pval_cut) %>%
    slice_min(order_by=Adjusted.P.value,
              n=10, with_ties=FALSE) %>%
    separate_rows(DiseaseNames, sep=", ") %>%
    select(Term, Data_Type, DiseaseNames, DiseaseCount)
  
  combi_data  <- bind_rows(combi_data,  data_df)
  combi_count <- bind_rows(combi_count, bar_df)
}
# We changed the names, making them more pleasing
label_map <- c(
  GO_Molecular_Function_2021 = "GO Molecular Function",
  GO_Cellular_Component_2021 = "GO Cellular Component",
  GO_Biological_Process_2021 = "GO Biological Process",
  KEGG_2021_Human            = "KEGG Human",
  InterPro_Domains_2019      = "InterPro Domains"
)
combi_data  <- combi_data  %>% mutate(Data_Type = recode(Data_Type, !!!label_map))
combi_count <- combi_count %>% mutate(Data_Type = recode(Data_Type, !!!label_map))

# Now we come to the plotting part
# bubble plot
bubble <- goplot(combi_data,
                 top=10,
                 pvalue=pval_cut)


# Bar plot, we compute per-term disease counts, then find the overall max:
term_counts <- combi_count %>%
  group_by(Data_Type, Term) %>%
  tally(name = "DiseaseCount") %>%
  ungroup()

max_d <- max(term_counts$DiseaseCount)

# then prepare a plotting version of combi_count with wrapped & reordered term labels
plot_df <- combi_count %>%
  group_by(Data_Type, Term) %>%
  mutate(DiseaseCount = n()) %>%        # each row is one disease, so n() gives count
  ungroup() %>%
  mutate(
    PlotTerm = str_wrap(Term, 30),      # wrap long terms
    PlotTerm = fct_reorder(PlotTerm, DiseaseCount)
  ) %>%
  mutate(Value = 1)                     # each row contributes one unit

# after we can draw the bar‐plot
bar <- ggplot(plot_df, aes(
  x    = Value,
  y    = PlotTerm,
  fill = DiseaseNames
)) +
  geom_col() +
  facet_wrap(~ Data_Type,
             scales = "free_y",
             nrow   = 5) +
  scale_fill_manual(
    name   = "Disease",
    values = distinct_colors_21,
    guide  = guide_legend(
      ncol      = 1,
      byrow     = TRUE,
      keywidth  = unit(1, "cm"),
      keyheight = unit(0.5, "cm")
    ),
    drop = FALSE
  ) +
  scale_x_continuous(
    name   = "Number of Diseases",
    limits = c(0, max_d),
    breaks = 0:max_d,
    expand = c(0,0)
  ) +
  labs(y = NULL) +
  theme_bw() +
  theme(
    strip.text      = element_text(size = 10),
    axis.text.y     = element_blank(),
    axis.ticks.y    = element_blank(),
    axis.title.x    = element_text(size = 10),
    legend.position = "right",
    legend.text     = element_text(size = 8),
    legend.title    = element_text(size = 10)
  )

# and combine and save
fig2 <- bubble + bar +
  plot_layout(guides="collect") +
  plot_annotation(tag_levels="A")

print(fig2)

dir.create("results", showWarnings=FALSE)
ggsave("results/enrichR_autoimmune.png",
       fig2,
       width=30, height=40, units="cm", dpi=400)

# Now we move on by splitting the group into two, the Matched and the Unmatched in regards to the healthy peptides
# We make use of the UniProt files, and then much of the process is repeated for the two new figures

matched_df   <- read_delim("protein_data/uniprotkb_matched_proteins.csv",
                           delim = ";", col_types = cols())
unmatched_df <- read_delim("protein_data/uniprotkb_unmatched_proteins.csv",
                           delim = ";", col_types = cols())

matched_entries   <- matched_df$Entry
unmatched_entries <- unmatched_df$Entry

# we Split IEDB into two subsets
iedb_matched   <- iedb %>% filter(`Epitope: Parent Protein Accession` %in% matched_entries)
iedb_unmatched <- iedb %>% filter(`Epitope: Parent Protein Accession` %in% unmatched_entries)

# then define a helper to run the full enrich→plot pipeline
run_enrichment <- function(iedb_subset, title_prefix){
  # build disease gene map and gene list
  dg   <- iedb_subset %>% select(DiseaseName, GeneSymbol) %>% distinct()
  gl   <- unique(dg$GeneSymbol)
  res  <- enrichr(gl, databases = dbs)
  
  # combine & count
  cd   <- tibble(); cc <- tibble()
  for(db in dbs){
    cnt <- res[[db]] %>%
      separate_rows(Genes, sep=";") %>%
      left_join(dg, by=c("Genes"="GeneSymbol"), relationship="many-to-many") %>%
      filter(!is.na(DiseaseName)) %>%
      group_by(Term) %>%
      summarise(
        DiseaseCount = n_distinct(DiseaseName),
        DiseaseNames = paste(unique(DiseaseName), collapse=", ")
      ) %>%
      ungroup() %>%
      mutate(Data_Type = db)
    
    dat <- res[[db]] %>%
      left_join(cnt, by="Term") %>%
      filter(DiseaseCount>1) %>%
      mutate(Data_Type = db)
    
    bar_tmp <- dat %>%
      filter(Adjusted.P.value < pval_cut) %>%
      slice_min(order_by=Adjusted.P.value, n=10, with_ties=FALSE) %>%
      separate_rows(DiseaseNames, sep=", ") %>%
      select(Term, Data_Type, DiseaseNames)
    
    cd <- bind_rows(cd,   dat)
    cc <- bind_rows(cc, bar_tmp)
  }
  
  
  cd <- cd %>% mutate(Data_Type = recode(Data_Type, !!!label_map))
  cc <- cc %>% mutate(Data_Type = recode(Data_Type, !!!label_map))
  
  # bubble plot
  bubble <- goplot(cd,
                   top    = 10,
                   pvalue = pval_cut)
  
  # bar plot
  maxd <- cc %>% group_by(Data_Type, Term) %>% tally(name="n") %>% pull(n) %>% max()
  pdf  <- cc %>%
    group_by(Data_Type, Term) %>%
    mutate(DiseaseCount = n()) %>%
    ungroup() %>%
    mutate(
      PlotTerm = str_wrap(Term, 30),
      PlotTerm = fct_reorder(PlotTerm, DiseaseCount),
      Value    = 1
    )
  bar <- ggplot(pdf, aes(x=Value, y=PlotTerm, fill=DiseaseNames)) +
    geom_col() +
    facet_wrap(~ Data_Type, scales="free_y", nrow=5) +
    scale_fill_manual(
      name   = "Disease",
      values = distinct_colors_21,
      guide  = guide_legend(ncol=1, byrow=TRUE,
                            keywidth=unit(1,"cm"),
                            keyheight=unit(0.5,"cm")),
      drop = FALSE
    ) +
    scale_x_continuous(
      name   = "Number of Diseases",
      limits = c(0, maxd),
      breaks = 0:maxd,
      expand = c(0,0)
    ) +
    labs(y=NULL) +
    theme_bw() +
    theme(
      strip.text      = element_text(size=10),
      axis.text.y     = element_blank(),
      axis.ticks.y    = element_blank(),
      axis.title.x    = element_text(size=10),
      legend.position = "right",
      legend.text     = element_text(size=8),
      legend.title    = element_text(size=10)
    )
  
  # we combine and return
  patchwork::wrap_plots(
    bubble + bar +
      plot_layout(guides="collect") +
      plot_annotation(tag_levels="A",
                      title = paste("Enrichment for", title_prefix)),
    ncol = 1
  )
}

# then run for matched & unmatched
fig_matched   <- run_enrichment(iedb_matched,   "Matched Proteins")
fig_unmatched <- run_enrichment(iedb_unmatched, "Unmatched Proteins")

# lastly print & save
print(fig_matched)
print(fig_unmatched)
ggsave("results/Matched_Enrichment.png",   fig_matched,   width=30, height=40, units="cm", dpi=400)
ggsave("results/Unmatched_Enrichment.png", fig_unmatched, width=30, height=40, units="cm", dpi=400)

# Now, using the same method, we tried to test the bigger dataset with healthy self peptides.
# For this we also downloaded the human proteosome from UniProt to help this analysis

supp <- read_csv("unique_supplementary_proteins.csv", show_col_types = FALSE)
genes_sup <- supp$ProteinID
message("▶︎ Testing ", length(genes_sup), " self-peptide proteins")

# ▶︎ Testing 6810 self-peptide proteins

iedb <- read_csv("autoimmune_iedb_pepetides.csv", show_col_types = FALSE) %>%
  # We fix the trailing semicolons in the DiseaseName column
  rename(DiseaseName = `DiseaseName;;;;;;;;;;;;;`) %>%
  mutate(DiseaseName = str_remove(DiseaseName, ";+$"))

# build a distinct disease↔gene table
disease_gene <- iedb %>%
  select(DiseaseName, GeneSymbol) %>%
  distinct()

# and list of unique autoimmune genes
auto_genes <- unique(disease_gene$GeneSymbol)
message("▶︎ Autoimmune genes: ", length(auto_genes))

# ▶︎ Autoimmune genes: 116

supp_unique <- read_csv("unique_supplementary_proteins.csv", show_col_types = FALSE)

# load UniProt human proteome (with GeneSymbol column)
uniprot_human <- read_tsv("rawdata/Uniprot/human_proteosome.tsv",
                          show_col_types = FALSE)

# map the ProteinIDs → GeneSymbol
healthy_genes <- uniprot_human %>%
  filter(Entry %in% supp_unique$ProteinID) %>%
  pull("Gene Names") %>%
  unique()
message("▶︎ Healthy self-peptide genes: ", length(healthy_genes))

# ▶︎ Healthy self-peptide genes: 5105
# As we can see not all is actually included

dbs <- c(
  "GO_Molecular_Function_2021",
  "GO_Cellular_Component_2021",
  "GO_Biological_Process_2021",
  "KEGG_2021_Human",
  "InterPro_Domains_2019"
)

res_auto    <- enrichr(auto_genes,    databases = dbs)
res_healthy <- enrichr(healthy_genes, databases = dbs)

label_map <- c(
  GO_Molecular_Function_2021 = "GO Molecular Function",
  GO_Cellular_Component_2021 = "GO Cellular Component",
  GO_Biological_Process_2021 = "GO Biological Process",
  KEGG_2021_Human            = "KEGG Human",
  InterPro_Domains_2019      = "InterPro Domains"
)

tidy_res <- function(res_list, grp){
  map_dfr(names(res_list), function(db){
    res_list[[db]] %>%
      mutate(
        Data_Type = recode(db, !!!label_map),
        Group     = grp
      ) %>%
      # we split the Overlap into numeric Count/Total
      separate(col = Overlap, into = c("Count","Total"), sep = "/", remove = FALSE) %>%
      mutate(
        Count     = as.integer(Count),
        GeneRatio = Count / as.numeric(Total)
      )
  })
}

comp_df <- bind_rows(
  tidy_res(res_healthy,  "Healthy Self"),
  tidy_res(res_auto,     "Autoimmune IEDB")
)

# We filter for significant hits
comp_sig <- comp_df %>%
  filter(Adjusted.P.value < 0.05)

# Now we go mapping the top10 healthy
healthy_top10 <- comp_df %>%                                       
  filter(Group == "Healthy Self", Adjusted.P.value < 0.05) %>%     
  group_by(Data_Type) %>%
  slice_min(order_by = Adjusted.P.value, n = 10, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    # we clean the facet labels
    Data_Type = recode(Data_Type,
                       GO_Molecular_Function_2021 = "GO Molecular Function",
                       GO_Cellular_Component_2021 = "GO Cellular Component",
                       GO_Biological_Process_2021 = "GO Biological Process",
                       KEGG_2021_Human            = "KEGG Human",
                       InterPro_Domains_2019      = "InterPro Domains"
    ),
    # Then we wrap each term to a single line
    Term      = str_wrap(Term, width = 50)
  )

# Finally for the plot
p_healthy <- ggplot(healthy_top10, aes(
  x     = reorder(Term, GeneRatio),
  y     = GeneRatio,
  size  = Count,
  color = Adjusted.P.value
)) +
  geom_point() +
  facet_wrap(~ Data_Type, ncol = 1, scales = "free_y") +
  scale_color_viridis(option = "D", name = "Adj. P-value") +
  scale_size_area(name = "Gene Count", max_size = 6) +
  coord_flip() +
  labs(
    x     = NULL,
    y     = "Gene Ratio"
  ) +
  theme_bw() +
  theme(
    plot.title      = element_text(face = "bold", hjust = 0.5),
    strip.text      = element_text(face = "bold"),
    axis.text.x     = element_text(angle = 45, hjust = 1),
    axis.text.y     = element_text(size = 8),
    legend.position = "right",
    plot.margin     = unit(c(1,1,1,1), "cm")
  )

# We display and save
print(p_healthy)

dir.create("results", showWarnings = FALSE)
ggsave(
  filename = "results/healthy_self_enrichment_top10.png",
  plot     = p_healthy,
  width    = 20,
  height   = 25,
  units    = "cm",
  dpi      = 300
)

# We saw that very very little significant evidence, so we will perhaps test this a bit further.

#____________________________________________________________________________________________________________

# Let us now run in the csv's and since we also have the tsv, we will begin to test aa composition, and net charge

self <- read_csv("unique_supplementary_proteins.csv", show_col_types = FALSE) %>%
  rename(ProteinID = ProteinID) %>%
  mutate(Group = "Healthy Self")

auto <- read_csv("unique_autoimmune_proteins.csv", show_col_types = FALSE) %>%
  rename(ProteinID = ProteinID) %>%
  mutate(Group = "Autoimmune IEDB")

all_ids <- bind_rows(self, auto)
message("▶︎ Healthy: ",   sum(all_ids$Group=="Healthy Self"),
        "   Autoimmune: ", sum(all_ids$Group=="Autoimmune IEDB"))

# ▶︎ Healthy: 6810   Autoimmune: 116

# We pull the information
prot <- read_tsv("rawdata/Uniprot/human_proteosome.tsv",
                 show_col_types = FALSE)

all_seqs <- prot %>%
  filter(Entry %in% all_ids$ProteinID) %>%
  select(Entry, Sequence) %>%
  left_join(all_ids, by = c("Entry" = "ProteinID"))

# Compute per-protein for the aa composition
aa_long <- all_seqs %>%
  mutate(aa_tab = map(Sequence, aa_composition)) %>%
  select(Entry, Group, aa_tab) %>%
  unnest(aa_tab) %>%
  select(-sequence)

# Let us summarise
aa_summary <- aa_long %>%
  pivot_longer(cols = -c(Entry, Group),
               names_to  = "AminoAcid",
               values_to = "Count") %>%
  group_by(Group, AminoAcid) %>%
  summarise(MeanCount = mean(Count), .groups="drop") %>%
  arrange(Group, desc(MeanCount))

# Let us try to plot it side by side and save it
p_compare_aa <- ggplot(aa_summary, aes(
  x     = MeanCount,
  y     = reorder(AminoAcid, MeanCount),
  fill  = Group
)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = c("Healthy Self"="#2c7fb8",
                               "Autoimmune IEDB"="#e41a1c")) +
  coord_flip() +
  labs(
    x     = NULL,
    y = "Mean Count per Protein",
    fill  = NULL
  ) +
  theme_minimal(base_size = 18) +
  theme(
    plot.title  = element_text(face="bold", hjust=0.5),
    axis.text.y = element_text(size=8),
    legend.position = "right"
  )

print(p_compare_aa)

ggsave(
  "results/aa_composition_self_vs_autoimmune.png",
  p_compare_aa,
  width  = 16,
  height = 12,
  dpi    = 400
)

# We saw little significance from the test, but will for a good test also check if we can see something from the net charge analysis

self_ids <- read_csv("unique_supplementary_proteins.csv", show_col_types=FALSE) %>%
  rename(ProteinID = ProteinID) %>%
  mutate(Group = "Healthy Self")

auto_ids <- read_csv("unique_autoimmune_proteins.csv", show_col_types=FALSE) %>%
  rename(ProteinID = ProteinID) %>%
  mutate(Group = "Autoimmune IEDB")

all_ids <- bind_rows(self_ids, auto_ids)

proteome <- read_tsv("rawdata/Uniprot/human_proteosome.tsv", show_col_types=FALSE)

seq_df <- proteome %>%
  filter(Entry %in% all_ids$ProteinID) %>%
  select(Entry, Sequence) %>%
  left_join(all_ids, by=c("Entry"="ProteinID"))

# This is the part where we compute the net charge at a pH of 7
charge_df <- seq_df %>%
  mutate(
    net_charge = map_dbl(Sequence, ~ charge(.x, pH=7.0))
  ) %>%
  select(Entry, Group, net_charge)

summary_stats <- charge_df %>%
  group_by(Group) %>%
  summarize(
    n           = n(),
    mean_charge = mean(net_charge),
    sd_charge   = sd(net_charge),
    .groups="drop"
  )

print(summary_stats)

# Group         n mean_charge sd_charge
# <chr>     <int>       <dbl>     <dbl>
# Autoimmu…   116       -5.00      38.1
# Healthy …  5109       -5.95      24.7

# We used a Welch two sample t test
ttest_res <- t.test(net_charge ~ Group, data=charge_df)
tidy(ttest_res)

# Even though we see little significance from the p value, let us still plot
p_charge <- ggplot(charge_df, aes(x=Group, y=net_charge, fill=Group)) +
  geom_violin(trim=TRUE, alpha=0.6) +
  geom_jitter(width=0.15, size=1, alpha=0.7) +
  scale_fill_manual(values = c("Healthy Self"="#377eb8", "Autoimmune IEDB"="#e41a1c")) +
  labs(
    x     = NULL,
    y     = "Net Charge"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    plot.title     = element_text(face="bold", hjust=0.5),
    legend.position= "none"
  )

print(p_charge)

ggsave("results/net_charge_comparison.png",
       p_charge,
       width  = 14, height = 10, dpi=400)

# Now we move on

#____________________________________________________________________________________________________


# We formerly had a problem with the gene enrichment analysis for the healthy dataset.
# This time we are going to make use of ClusterProfiler, org.Hs.eg.db, AnnotationDbi, UniProt.ws and ReactomePA
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(UniProt.ws)

autoimmune <- read_csv("autoimmune_iedb_pepetides.csv",
                       show_col_types = FALSE)
healthy    <- read_csv("supplementary_table1.csv",
                       show_col_types = FALSE)

# We extract the unique extensions

auto_ids   <- autoimmune  %>% pull(`Epitope: Parent Protein Accession`) %>% unique() %>% str_remove("\\.\\d+$")
healthy_ids<- healthy %>% pull(Source_proteinID) %>% unique() %>% str_remove("\\.\\d+$")

all_ids    <- unique(c(auto_ids, healthy_ids))

# We build a UniProt.ws object (human)
up <- UniProt.ws(taxId = 9606)
# And See what columns are actually available for the keytype “UniProtKB”:
columns(up) %>% head(100)

# We check which keytypes we can actually use
valid_kts <- keytypes(up)
print(valid_kts)

# [1] "Allergome"                    
# [2] "ArachnoServer"                
# [3] "Araport"                      
# [4] "BioCyc"                       
# [5] "BioGRID"                      
# [6] "BioMuta"                      
# [7] "CCDS"                         
# [8] "CGD"                          
# [9] "ChEMBL"   
# And so on

# Now we map to gene-symbol, EntrezID & sequence
up    <- UniProt.ws(taxId=9606)
cols  <- c("UniProtKB","gene_names","gene_primary")
batches <- split(all_ids, ceiling(seq_along(all_ids)/500))

mapping_up <- purrr::map_dfr(batches, function(batch) {
  tryCatch(
    AnnotationDbi::select(up,
                          keys    = batch,
                          columns = cols,
                          keytype = "UniProtKB"
    ),
    error = function(e){
      warning("batch failed, returning empty tibble: ", e$message)
      tibble()
    }
  )
})

# This part does take a bit of time, and we are pulling the data in batches of 500 due to the API
# Now let us map them to the group they are in

mapping_labeled <- mapping_up %>%
  # first, we keep only those rows whose original key was one of ours
  filter(From %in% c(auto_ids, healthy_ids)) %>%
  # now we tag by which list they came from
  mutate(
    match_status = case_when(
      From %in% auto_ids    ~ "Autoimmune",
      From %in% healthy_ids ~ "Healthy"
    )
  )

# We do a sanity‐check: 
mapping_labeled %>% count(match_status)

#  match_status    n
#    Autoimmune  116
#       Healthy 5148

# We extract the symbol list, and see that some of the many health proteins are removed

genes_auto    <- mapping_labeled %>%
  filter(match_status=="Autoimmune") %>%
  pull(Gene.Names..primary.) %>%
  unique()

genes_healthy <- mapping_labeled %>%
  filter(match_status=="Healthy") %>%
  pull(Gene.Names..primary.) %>%
  unique()

length(genes_auto)    
length(genes_healthy) 

# 116
# 5058

# Now we make a helper to run and plot an enrichGO
do_enrichGO <- function(gene.vector, ont, title) {
  res <- enrichGO(
    gene         = gene.vector,
    OrgDb        = org.Hs.eg.db,
    keyType      = "SYMBOL",
    ont          = ont,
    pAdjustMethod= "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05
  )
  dotplot(res, showCategory=10) +
    ggtitle(title) +
    scale_color_continuous(trans="reverse") +
    labs(x = "-log10(adj.p)", color="adj.p")
}

# We also make a helper for the reactome
# helper for Reactome
do_pathway <- function(gene.vector, title) {
  res <- enrichPathway(
    gene         = gene.vector,
    organism     = "human",
    pAdjustMethod= "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable     = TRUE
  )
  dotplot(res, showCategory=10) +
    ggtitle(title) +
    labs(x = "-log10(adj.p)", color="adj.p")
}

# We make sure R uses Cairo for bitmaps
options(bitmapType="cairo")

# Let us run this, first the symbol lists
genes_auto    <- mapping_labeled %>% filter(match_status=="Autoimmune") %>% pull(Gene.Names..primary.) %>% unique()
genes_healthy <- mapping_labeled %>% filter(match_status=="Healthy")    %>% pull(Gene.Names..primary.) %>% unique()

# We convert to Entrez (this is for Reactome only)
auto_e    <- bitr(genes_auto,    fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)$ENTREZID
healthy_e <- bitr(genes_healthy, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)$ENTREZID

# Now we have the helper for enrichGO (symbols)
do_enrichGO <- function(genes, ont, title){
  res <- enrichGO(gene = genes,
                  OrgDb = org.Hs.eg.db,
                  keyType = "SYMBOL",
                  ont = ont,
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.05)
  dotplot(res, showCategory=10) +
    ggtitle(title) +
    labs(x="-log10(adj.p)", color="adj.p") +
    scale_color_continuous(trans="reverse")
}

# And the helper for Reactome (Entrez)
do_pathway <- function(eids, title){
  res <- enrichPathway(gene = eids,
                       organism = "human",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05,
                       readable = TRUE)
  dotplot(res, showCategory=10) +
    ggtitle(title) +
    labs(x="-log10(adj.p)", color="adj.p")
}

# Now we run and plot
bp_auto    <- do_enrichGO(genes_auto,    "BP", "GO–BP: Autoimmune")
bp_healthy <- do_enrichGO(genes_healthy, "BP", "GO–BP: Healthy")

mf_auto    <- do_enrichGO(genes_auto,    "MF", "GO–MF: Autoimmune")
mf_healthy <- do_enrichGO(genes_healthy, "MF", "GO–MF: Healthy")

cc_auto    <- do_enrichGO(genes_auto,    "CC", "GO–CC: Autoimmune")
cc_healthy <- do_enrichGO(genes_healthy, "CC", "GO–CC: Healthy")

pw_auto    <- do_pathway(auto_e,    "Reactome: Autoimmune")
pw_healthy <- do_pathway(healthy_e, "Reactome: Healthy")

# We combine the figures, plot and save
combined_enrichment_auto <- 
  (bp_auto    | mf_auto)    /
  (cc_auto    | pw_auto) +
  plot_layout(
    guides   = "collect",
  ) +
  plot_annotation(
  ) &
  theme_minimal(base_size = 18) &        # bump text size across all panels
  theme(
    plot.title      = element_text(size = 20, face = "bold"),
    legend.position = "right"
  )

combined_enrichment_Healthy <- 
  (bp_healthy    | mf_healthy)    /
  (cc_healthy    | pw_healthy) +
  plot_layout(
    guides   = "collect",
  ) +
  plot_annotation(
  ) &
  theme_minimal(base_size = 18) &        # bump text size across all panels
  theme(
    plot.title      = element_text(size = 20, face = "bold"),
    legend.position = "right"
  )

ggsave("results/enrichment_autoimmune.png",
       combined_enrichment_auto, width = 20, height = 18, dpi = 300)

ggsave("results/enrichment_healthy.png",
       combined_enrichment_Healthy, width = 20, height = 18, dpi = 300)

# We have now also been able to complete a gene enrichment, that showcases the difference between the autoimmune and healthy
