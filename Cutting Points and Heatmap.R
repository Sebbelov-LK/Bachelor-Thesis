iedb_peptides <- read.csv("autoimmune_iedb_pepetides.csv", stringsAsFactors = FALSE)

######################## Extended IEDB Peptides ######################################

# Load required libraries
library(httr)

# Load data
iedb_peptides <- read.csv("autoimmune_iedb_pepetides.csv", stringsAsFactors = FALSE)

# Define UniProt sequence fetcher
fetch_uniprot_sequence <- function(uniprot_id) {
  url <- paste0("https://rest.uniprot.org/uniprotkb/", uniprot_id, ".fasta")
  response <- GET(url)
  
  if (status_code(response) == 200) {
    fasta_text <- content(response, as = "text", encoding = "UTF-8")
    lines <- strsplit(fasta_text, "\n")[[1]]
    seq <- paste(lines[-1], collapse = "")  # Skip header
    return(seq)
  } else {
    return(NA)
  }
}

# Define peptide extension function
extend_peptide <- function(peptide, full_seq) {
  match_pos <- regexpr(peptide, full_seq)[1]
  if (match_pos == -1) return(NA)
  
  start <- max(1, match_pos - 5)
  end <- min(nchar(full_seq), match_pos + nchar(peptide) + 4)
  substr(full_seq, start, end)
}

# Initialize vector for extended peptides
extended_peptides <- character(nrow(iedb_peptides))

# Loop through rows and extend each peptide
for (i in seq_len(nrow(iedb_peptides))) {
  peptide <- iedb_peptides[[1]][i]
  uniprot_id <- iedb_peptides[[4]][i]
  
  if (!is.na(peptide) && !is.na(uniprot_id) && !grepl("^NP_", uniprot_id)) {
    protein_seq <- fetch_uniprot_sequence(uniprot_id)
    if (!is.na(protein_seq)) {
      extended_peptides[i] <- extend_peptide(peptide, protein_seq)
    } else {
      extended_peptides[i] <- NA
    }
  } else {
    extended_peptides[i] <- NA
  }
}

# Add to dataframe
iedb_peptides$Extended_Peptide <- extended_peptides

# Preview result
head(iedb_peptides)

# Save to CSV
write.csv(iedb_peptides, "IEDB_Peptides_Extended.csv", row.names = FALSE)


################### ICElogo ###################
library(ggseqlogo)
# Load the extended data (from the previously saved CSV)
iedb_peptides <- read.csv("IEDB_Peptides_Extended.csv", stringsAsFactors = FALSE)

# Filter valid extended sequences (min 10 AAs)
valid_seqs <- na.omit(iedb_peptides$Extended_Peptide)
valid_seqs <- valid_seqs[nchar(valid_seqs) >= 10]

# First 10 AAs (N-terminal)
first_10 <- substr(valid_seqs, 1, 10)

# Last 10 AAs (C-terminal)
last_10 <- substr(valid_seqs, nchar(valid_seqs) - 9, nchar(valid_seqs))

# Create sequence logos
logo_first <- ggseqlogo(first_10, method = "prob") +
  ggtitle("N-terminal (First 10 AAs) - Extended Autoimmune Peptides")

logo_last <- ggseqlogo(last_10, method = "prob") +
  ggtitle("C-terminal (Last 10 AAs) - Extended Autoimmune Peptides")

# Save to files
ggsave("ICElogo_N_terminal.png", plot = logo_first, width = 8, height = 4, dpi = 300)
ggsave("ICElogo_C_terminal.png", plot = logo_last, width = 8, height = 4, dpi = 300)


####################### Save N terminal and C terminal ###################
# Ensure valid sequences (non-NA and long enough)
valid_seqs <- na.omit(iedb_peptides$Extended_Peptide)
valid_seqs <- valid_seqs[nchar(valid_seqs) >= 10]

# Extract first and last 10 amino acids
first_10 <- substr(valid_seqs, 1, 10)
last_10 <- substr(valid_seqs, nchar(valid_seqs) - 9, nchar(valid_seqs))

# Write to .txt files
writeLines(first_10, "First10_Extended_Autoimmune.txt")
writeLines(last_10, "Last10_Extended_Autoimmune.txt")


########################## Extended self-peptides #######################
# Load required file
supp_table <- read.csv("supplementary_table1 (1).csv", stringsAsFactors = FALSE)

# Define a regular expression to match any invalid amino acid
invalid_aa_pattern <- "[BJOUXZ]"

# Filter peptides: >5 AAs, no invalid AAs, and UniProt IDs not starting with "NP_"
supp_filtered <- supp_table[
  !grepl("^NP_", supp_table[[3]]) &
    !is.na(supp_table[[3]]) &
    supp_table[[3]] != "" &
    supp_table[[2]] != "" &
    !grepl(invalid_aa_pattern, supp_table[[2]]) &   # Exclude if peptide contains invalid AAs
    nchar(supp_table[[2]]) > 5,
]

# Load autoimmune peptides and determine target count
autoimmune_peptides <- read.csv("IEDB_Peptides_Extended.csv", stringsAsFactors = FALSE)
valid_autoimmune <- na.omit(autoimmune_peptides$Extended_Peptide)
valid_autoimmune <- valid_autoimmune[nchar(valid_autoimmune) >= 10]
autoimmune_count <- length(valid_autoimmune)

# Shuffle to randomize sampling
set.seed(42)
supp_filtered <- supp_filtered[sample(nrow(supp_filtered)), ]

# Fetch and extend self peptides (ensure uniqueness, valid length, and no invalid AAs)
extended_self <- character()
seen_seqs <- character()

for (i in seq_len(nrow(supp_filtered))) {
  if (length(extended_self) >= autoimmune_count) break
  
  peptide <- supp_filtered[[2]][i]
  uniprot_id <- supp_filtered[[3]][i]
  
  seq <- fetch_uniprot_sequence(uniprot_id)
  if (!is.na(seq)) {
    extended <- extend_peptide(peptide, seq)
    
    if (!is.na(extended) &&
        nchar(extended) >= 10 &&
        !grepl(invalid_aa_pattern, extended) &&    # Exclude if extended seq contains invalid AAs
        !(extended %in% seen_seqs)) {
      
      extended_self <- c(extended_self, extended)
      seen_seqs <- c(seen_seqs, extended)
    }
  }
}

# Extract N- and C-terminal 10-mers
self_first10 <- substr(extended_self, 1, 10)
self_last10  <- substr(extended_self, nchar(extended_self) - 9, nchar(extended_self))

# Save results
writeLines(extended_self, "Extended_Self_Peptides.txt")
writeLines(self_first10,   "First10_Self_Peptides.txt")
writeLines(self_last10,    "Last10_Self_Peptides.txt")
#################### Self-peptide ICElogo ##################

# Load the ggseqlogo package
library(ggseqlogo)

# Load N- and C-terminal 10-mers from self-peptides
self_first10 <- readLines("First10_Self_Peptides.txt")
self_last10  <- readLines("Last10_Self_Peptides.txt")

# Plot first 10 AAs (N-terminal)
self_logo_first <- ggseqlogo(self_first10, method = "prob") +
  ggtitle("ICElogo-style Plot: First 10 AAs of Self-Peptides")

# Plot last 10 AAs (C-terminal)
self_logo_last <- ggseqlogo(self_last10, method = "prob") +
  ggtitle("ICElogo-style Plot: Last 10 AAs of Self-Peptides")

# Save the logo plots
ggsave("Self_First10_ICElogo.png", plot = self_logo_first, width = 8, height = 4, dpi = 300)
ggsave("Self_Last10_ICElogo.png",  plot = self_logo_last, width = 8, height = 4, dpi = 300)


############################### Task 3 ########################################

library(httr)
library(ggplot2)

# Load autoimmune peptide data
iedb <- read.csv("autoimmune_iedb_pepetides.csv", stringsAsFactors = FALSE)

# Function: Fetch UniProt protein sequence
fetch_uniprot_sequence <- function(uniprot_id) {
  url <- paste0("https://rest.uniprot.org/uniprotkb/", uniprot_id, ".fasta")
  response <- GET(url)
  if (status_code(response) == 200) {
    fasta_text <- content(response, as = "text", encoding = "UTF-8")
    lines <- strsplit(fasta_text, "\n")[[1]]
    paste(lines[-1], collapse = "")
  } else {
    NA
  }
}

# Loop to collect normalized positions
normalized_positions <- c()

for (i in seq_len(nrow(iedb))) {
  peptide <- iedb[[1]][i]
  uniprot_id <- iedb[[4]][i]
  
  if (!is.na(peptide) && !is.na(uniprot_id) && !grepl("^NP_", uniprot_id)) {
    full_seq <- fetch_uniprot_sequence(uniprot_id)
    if (!is.na(full_seq)) {
      match_pos <- regexpr(peptide, full_seq)[1]
      prot_len <- nchar(full_seq)
      if (match_pos != -1 && prot_len > 0) {
        norm_pos <- match_pos / prot_len
        normalized_positions <- c(normalized_positions, norm_pos)
      }
    }
  }
}

# Bin positions manually
binwidth <- 0.05
bins <- cut(normalized_positions, breaks = seq(0, 1, by = binwidth), include.lowest = TRUE)
bin_levels <- levels(bins)
freqs <- table(bins)

# Midpoints for x-axis
midpoints <- seq(binwidth / 2, 1 - binwidth / 2, by = binwidth)

# Assemble data frame for plotting
plot_data <- data.frame(
  Bin = bin_levels,
  Midpoint = midpoints,
  Frequency = as.numeric(freqs)
)

# Plot with red-to-green fill gradient
p <- ggplot(plot_data, aes(x = Midpoint, y = Frequency, fill = Frequency)) +
  geom_col(color = "black") +
  scale_fill_gradient(low = "red", high = "green") +
  labs(
    title = "Autoimmune Peptide Frequency by Normalized Protein Position",
    x = "Normalized Protein Length",
    y = "Peptide Count"
  ) +
  theme_minimal()

# Save plot
ggsave("Autoimmune_Peptide_HeatHistogram.png", plot = p, width = 8, height = 4, dpi = 300)


############# Heatmap for self-peptides(the sample from last) ##################
# Map peptide positions in normalized space
normalized_positions <- c()
seen_seqs <- character()

for (i in seq_len(nrow(supp_filtered))) {
  if (length(normalized_positions) >= autoimmune_count) break
  
  peptide <- supp_filtered[[2]][i]
  uniprot_id <- supp_filtered[[3]][i]
  
  seq <- fetch_uniprot_sequence(uniprot_id)
  
  if (!is.na(seq) && !peptide %in% seen_seqs) {
    match_pos <- regexpr(peptide, seq)[1]
    
    if (match_pos != -1 && nchar(seq) > 0) {
      norm_pos <- match_pos / nchar(seq)
      normalized_positions <- c(normalized_positions, norm_pos)
      seen_seqs <- c(seen_seqs, peptide)
    }
  }
}

# Create histogram-style heatmap plot
df <- data.frame(Normalized_Position = normalized_positions)

heat_plot <- ggplot(df, aes(x = Normalized_Position)) +
  geom_histogram(bins = 50, aes(fill = ..count..)) +
  scale_fill_gradient(low = "red", high = "green") +
  labs(title = "Peptide Frequency by Normalized Protein Position (Self-Peptides)",
       x = "Normalized Position (0=start, 1=end)",
       y = "Frequency") +
  theme_minimal()

# Save
ggsave("Self_Peptides_Position_Heatmap.png", plot = heat_plot, width = 8, height = 4, dpi = 300)


##################### Save seed so i have the same ########################
# Extract the sampled data (all columns)
seed_sample <- supp_filtered[1:autoimmune_count, ]

# Save to CSV
write.csv(seed_sample, "Seed_Sample_Peptides.csv", row.names = FALSE)

