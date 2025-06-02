################## Name of longest peptide ###############################
# Extract unique autoimmune peptides (column 1)
unique_peptides <- unique(data[[1]])

# Find the longest one
longest_peptide <- unique_peptides[which.max(nchar(unique_peptides))]

# Print it
cat("🧬 The longest peptide is:\n", longest_peptide, "\n")# Load your cleaned data
data <- read.csv("Matched_Core_Peptides_Cleaned.csv", stringsAsFactors = FALSE)

# Extract unique autoimmune peptides (column 1)
unique_peptides <- unique(data[[1]])

# Find the longest one
longest_peptide <- unique_peptides[which.max(nchar(unique_peptides))]

# Print it
cat("🧬 The longest peptide is:\n", longest_peptide, "\n")
cat("📏 Length:", nchar(longest_peptide), "amino acids\n")

######### Avg length between self-peptide, autoimmune and combined #############
# Calculate average length for each group
avg_autoimmune <- mean(nchar(unique(data[[1]])))        # Autoimmune peptides (column 1)
avg_self        <- mean(nchar(unique(data[[2]])))        # Original self-peptides (column 2)
avg_cleaned     <- mean(nchar(unique(data$Cleaned_Self_Peptide)))  # Cleaned peptides

# Print results
cat("🧪 Average Peptide Lengths:\n")
cat(" - Autoimmune peptides: ", round(avg_autoimmune, 2), "aa\n")
cat(" - Matched self-peptides: ", round(avg_self, 2), "aa\n")
cat(" - Cleaned peptides (longest form): ", round(avg_cleaned, 2), "aa\n")

############# P-value on difference between cleaned and auto ###################
# Load the cleaned data
data <- read.csv("Matched_Core_Peptides_Cleaned.csv", stringsAsFactors = FALSE)

# Get unique sequences from relevant columns
unique_autoimmune <- unique(data[[1]])
unique_cleaned    <- unique(data[[7]])

# Get lengths
auto_lengths <- nchar(unique_autoimmune)
cleaned_lengths <- nchar(unique_cleaned)

# Perform t-test
t_result <- t.test(auto_lengths, cleaned_lengths)

# Show result
print(t_result)



############# Histograms ###############
# Get lengths
auto_lengths   <- nchar(unique(data[[1]]))
self_lengths   <- nchar(unique(data[[2]]))
cleaned_lengths <- nchar(unique(data$Cleaned_Self_Peptide))

# Create a combined dataframe
length_df <- data.frame(
  Length = c(auto_lengths, self_lengths, cleaned_lengths),
  Group = c(
    rep("Autoimmune", length(auto_lengths)),
    rep("Matched Self", length(self_lengths)),
    rep("Cleaned Self", length(cleaned_lengths))
  )
)

# Load ggplot2
library(ggplot2)

# Plot the histogram
plot <- ggplot(length_df, aes(x = Length, fill = Group)) +
  geom_histogram(binwidth = 1, position = "identity", alpha = 0.6, color = "black") +
  theme_minimal() +
  scale_fill_manual(values = c("Autoimmune" = "darkorange", "Matched Self" = "steelblue", "Cleaned Self" = "darkgreen")) +
  labs(
    title = "Peptide Length Distributions",
    x = "Peptide Length (amino acids)",
    y = "Count",
    fill = "Peptide Group"
  )

# Save the histogram
ggsave("Peptide_Length_Distributions.png", plot, width = 10, height = 6, dpi = 300)




######################### 3 different ######################
# Get lengths
auto_lengths    <- nchar(unique(data[[1]]))
self_lengths    <- nchar(unique(data[[2]]))
cleaned_lengths <- nchar(unique(data$Cleaned_Self_Peptide))

# Create a combined dataframe with group labels
length_df <- data.frame(
  Length = c(auto_lengths, self_lengths, cleaned_lengths),
  Group = factor(c(
    rep("Autoimmune", length(auto_lengths)),
    rep("Matched Self", length(self_lengths)),
    rep("Extended Self", length(cleaned_lengths))
  ), levels = c("Autoimmune", "Matched Self", "Extended Self"))
)

# Load ggplot2
library(ggplot2)

# Create faceted histograms
plot <- ggplot(length_df, aes(x = Length)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "black") +
  facet_wrap(~ Group, scales = "free_y", ncol = 1) +
  theme_minimal() +
  labs(
    title = "Length Distributions for Different Peptide Groups",
    x = "Peptide Length (amino acids)",
    y = "Count"
  )

# Save the plot
ggsave("Peptide_Length_Distributions_Faceted.png", plot, width = 8, height = 9, dpi = 300)

