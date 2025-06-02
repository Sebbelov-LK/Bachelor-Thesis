autoimmune_data <- read.csv("autoimmune_iedb_pepetides.csv")
supplementary_data <- read.csv("supplementary_table1 (1).csv")
inhouse_data <- read.csv("Inhouse_donor_hlas.csv")


# Find sekvenserne og vær sikker på at der ikke er nogen gentagelser
AI_seq <- autoimmune_data[[1]]
AI_seq <- unique(AI_seq)
self_seq <- supplementary_data[[2]]
self_seq <- unique(self_seq)


#Test hvor mange procent af peptiderne der kan blive fundet i naturlige peptider
AI_counts <- sapply(AI_seq, function(x) sum(self_seq == x))
non_zero_counts <- AI_counts[AI_counts > 0]
total_sum <- sum(non_zero_counts)
total_elements_AI_seq <- length(AI_seq)
percentage <- (total_sum / total_elements_AI_seq) * 100

#Jeg finder alle sekvenser der kan bliver ni-merer eller mere
AI_long <- AI_seq[nchar(AI_seq) >= 9]
self_long <- self_seq[nchar(self_seq) >= 9]

#Lav 9-merer
extract_9mers <- function(sequence) {
  n <- nchar(sequence)
  if (n < 9) return(character(0))  
  
  sapply(1:(n - 8), function(i) substr(sequence, i, i + 8))
}

expanded_AI_long <- unique(unlist(sapply(AI_long, extract_9mers)))
expanded_self_long <- unique(unlist(sapply(self_long, extract_9mers)))

#Opdatér original
AI_seq <- unique(expanded_AI_long)
self_seq <- unique(expanded_self_long)

#Test procent ligesom før
AI_counts <- sapply(AI_seq, function(x) sum(self_seq == x))
non_zero_counts <- AI_counts[AI_counts > 0]
total_sum <- sum(non_zero_counts)
total_elements_AI_seq <- length(AI_seq)
percentage <- (total_sum / total_elements_AI_seq) * 100
print(percentage)

#Put matchesne i en liste
matched_sequences <- intersect(AI_seq, self_seq)
writeLines(matched_sequences, "matched_sequences.txt")



########################## Forsøg på at sortere matchene sekvenser til donorer ##################
# Load supplementary_data
donors <- unique(supplementary_data[[1]])  # Extract unique donor names

# Function to extract 9-mers
extract_9mers <- function(sequence) {
  n <- nchar(sequence)
  if (n < 9) return(character(0))  # Skip short sequences
  sapply(1:(n - 8), function(i) substr(sequence, i, i + 8))
}

# Initialize an empty list to store donor-specific 9-mers
donor_9mers <- list()

# Loop through each donor and extract 9-mers
for (donor in donors) {
  # Filter sequences for the current donor
  donor_sequences <- supplementary_data[supplementary_data[[1]] == donor, 2]  # Assuming sequences are in column 2
  
  # Generate 9-mers for all sequences of the donor
  donor_9mers[[donor]] <- unique(unlist(sapply(donor_sequences, extract_9mers)))
}

# Print a summary of results
for (donor in names(donor_9mers)) {
  cat(paste(donor, "has", length(donor_9mers[[donor]]), "unique 9-mers\n"))
}

########################## Prøver at plotte distribution ##################
match_counts <- numeric(length(matched_sequences))

# Loop through matched_sequences and count occurrences across all donors
for (i in seq_along(matched_sequences)) {
  sequence <- matched_sequences[i]
  count <- sum(sapply(donor_9mers, function(donor_list) sequence %in% donor_list))
  match_counts[i] <- count
}

# Convert to a data frame for plotting
distribution_df <- as.data.frame(table(match_counts))
colnames(distribution_df) <- c("Times_Found", "Sequence_Count")

# Convert to numeric
distribution_df$Times_Found <- as.numeric(as.character(distribution_df$Times_Found))
distribution_df$Sequence_Count <- as.numeric(distribution_df$Sequence_Count)

# Create the plot
plot <- ggplot(distribution_df, aes(x = Times_Found, y = Sequence_Count)) +
  geom_bar(stat = "identity", fill = "steelblue", color = "black") +
  labs(
    title = "Distribution of Matched Sequence Occurrences Across Donors",
    x = "Number of Donors Found In (X-axis)",
    y = "Number of Sequences (Y-axis)"
  ) +
  theme_minimal()

# Save the plot using ggsave
ggsave("matched_sequence_donor_distribution.png", plot = plot, width = 8, height = 6, dpi = 300)

########################## Prøver Donor1 9-mers og dets HLA'er #########################
# Get unique matched 9-mers from Donor1
donor1_matched_9mers <- intersect(donor_9mers[["Donor1"]], matched_sequences)

# Extract full peptides from supplementary_data where Donor1's 9-mers appear
donor1_full_peptides <- unique(
  supplementary_data[supplementary_data[[1]] == "Donor1" & grepl(paste(donor1_matched_9mers, collapse = "|"), supplementary_data[[2]]), 2]
)

# Save to a text file (no repetitions)
writeLines(donor1_full_peptides, "Donor1_matched_full_peptides.txt")

######################## Laver HLA FIL ######################
# Ensure the first column is 'Donor' and all others are HLAs
colnames(inhouse_data)[1] <- "Donor"

# Combine all HLA columns into a single string per row, separated by spaces
inhouse_data$HLAs <- apply(inhouse_data[, -1, drop = FALSE], 1, function(row) {
  hla_list <- unique(na.omit(row))  # Remove NA values and ensure uniqueness
  return(paste(hla_list, collapse = " "))  # Join HLAs with spaces
})

# Create formatted output with "Donor HLAs" ensuring correct spacing
formatted_lines <- paste(inhouse_data$Donor, inhouse_data$HLAs, sep = " ")

# Ensure HLA entries are unique again after replacements
formatted_lines <- sapply(strsplit(formatted_lines, " "), function(x) {
  paste(unique(x), collapse = " ")  # Remove duplicate HLAs
})

# Convert back to character vector
formatted_lines <- as.character(formatted_lines)

# Replace formatting: spaces -> commas, semicolons -> spaces
formatted_lines <- gsub(" ", ",", formatted_lines)
formatted_lines <- gsub(";", " ", formatted_lines)

# Remove trailing commas (if any line ends with one)
formatted_lines <- gsub(",$", "", formatted_lines)

# Save to a .txt file
writeLines(formatted_lines, "Formatted_Donor_HLAs.txt")


####### Backtracker core peptides or sætter dem sammen med donorer ###################
matched_peptides_donors <- list()

# Loop through each matched sequence
for (seq in matched_sequences) {
  # Find the peptides that contain this matched sequence
  matching_rows <- supplementary_data[grepl(seq, supplementary_data[[2]]), ]
  
  # If matches found, store the full peptide and donor name
  if (nrow(matching_rows) > 0) {
    results <- paste(matching_rows[[2]], matching_rows[[1]])  # "Peptide Donor"
    matched_peptides_donors <- c(matched_peptides_donors, results)  # Append to list
  }
}

matched_peptides_donors <- as.character(matched_peptides_donors)

# Save to a text file
writeLines(matched_peptides_donors, "Matched_Peptides_Donors.txt")



