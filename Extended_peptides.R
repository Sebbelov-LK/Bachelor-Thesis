# Load the dataset
df <- read.csv("Matched_Core_Peptides.csv", stringsAsFactors = FALSE)

# Initialize the new column
df$Cleaned_Self_Peptide <- NA

# Get unique autoimmune peptides
unique_auto <- unique(df$Autoimmune_Peptide)

# Function to get longest coherent region covered by matched peptides
get_designated_region <- function(auto_peptide, matched_peptides) {
  positions <- lapply(matched_peptides, function(pep) {
    start <- regexpr(pep, auto_peptide)[1]
    if (start == -1) return(NULL)
    end <- start + nchar(pep) - 1
    return(c(start, end))
  })
  
  # Remove any NULLs (non-matching peptides)
  positions <- Filter(Negate(is.null), positions)
  if (length(positions) == 0) return(NA)
  
  # Get global min and max positions
  min_start <- min(sapply(positions, `[`, 1))
  max_end   <- max(sapply(positions, `[`, 2))
  
  # Extract substring from autoimmune peptide
  substr(auto_peptide, min_start, max_end)
}

# Loop over each autoimmune peptide group
for (auto_pep in unique_auto) {
  matched_rows <- df[df$Autoimmune_Peptide == auto_pep, ]
  matched_peps <- matched_rows$Matched_Peptide
  extended <- get_designated_region(auto_pep, matched_peps)
  
  df$Cleaned_Self_Peptide[df$Autoimmune_Peptide == auto_pep] <- extended
}

# Save to file
write.csv(df, "Matched_Core_Peptides_Cleaned.csv", row.names = FALSE)