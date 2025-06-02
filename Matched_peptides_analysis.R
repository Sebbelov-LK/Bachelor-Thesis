###############################################################################
####################### Overextensions Treeplots ##############################
###############################################################################

# Define amino acid biochemical properties
aa_properties <- data.frame(
  AA = c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"),
  Property = c(
    "Hydrophobic", "Positive", "Polar", "Negative", "Polar",
    "Polar", "Negative", "Small", "Positive", "Hydrophobic",
    "Hydrophobic", "Positive", "Hydrophobic", "Hydrophobic",
    "Cyclic", "Polar", "Polar", "Hydrophobic", "Polar", "Hydrophobic"
  ),
  stringsAsFactors = FALSE
)

# Function to summarize amino acid counts by biochemical property
summarize_by_property <- function(sequence_vector) {
  aa_vector <- unlist(strsplit(paste(sequence_vector, collapse = ""), split = ""))
  tbl <- table(aa_vector)
  df <- data.frame(AA = names(tbl), Count = as.vector(tbl))
  df <- merge(df, aa_properties, by = "AA")
  df %>% group_by(Property) %>% summarise(Count = sum(Count), .groups = "drop")
}

# Function to plot treemap
if (!requireNamespace("treemapify", quietly = TRUE)) install.packages("treemapify")
library(treemapify)
library(ggplot2)

plot_treemap <- function(df, title, filename) {
  p <- ggplot(df, aes(area = Count, fill = Property, label = paste0(Property, "\n", Count))) +
    geom_treemap() +
    geom_treemap_text(colour = "white", grow = TRUE, reflow = TRUE) +
    labs(title = title) +
    theme(legend.position = "none")
  ggsave(filename, plot = p, width = 8, height = 6, dpi = 300)
}

###############################################################################
###################### Load and Filter the Dataset ############################
###############################################################################

# Load the dataset
df <- read.csv("Peptide_Overlap_Decomposition.csv", stringsAsFactors = FALSE)

# Apply filtering to remove pairs with both Auto & Self prefixes or suffixes
filtered_df <- df %>%
  filter(!((Auto_Prefix != "" & Self_Prefix != "") |
             (Auto_Suffix != "" & Self_Suffix != "")))

# Show how many were kept and filtered out
remaining_rows <- nrow(filtered_df)
filtered_out_rows <- nrow(df) - remaining_rows
cat("Remaining peptide pairs after filtering:", remaining_rows, "\n")
cat("Number of peptide pairs filtered out:", filtered_out_rows, "\n")

###############################################################################
######################## Treemap Visualization ###############################
###############################################################################

# Autoimmune extensions from filtered data
auto_ext_sequences <- c(filtered_df$Auto_Prefix, filtered_df$Auto_Suffix)
auto_ext_summary <- summarize_by_property(auto_ext_sequences)
plot_treemap(auto_ext_summary, "Amino Acid Properties - Autoimmune Extensions", "Treemap_Autoimmune_Extensions.png")

# Self extensions from filtered data
if ("Self_Prefix" %in% names(filtered_df) && "Self_Suffix" %in% names(filtered_df)) {
  self_ext_sequences <- c(filtered_df$Self_Prefix, filtered_df$Self_Suffix)
  self_ext_summary <- summarize_by_property(self_ext_sequences)
  plot_treemap(self_ext_summary, "Amino Acid Properties - Self Extensions", "Treemap_Self_Extensions.png")
}

###############################################################################
####################### Center Distance Analysis ##############################
###############################################################################

filtered_df <- filtered_df %>%
  mutate(
    Auto_Center = floor(nchar(Autoimmune_Peptide) / 2),
    Self_Center = floor(nchar(Self_Peptide) / 2),
    Distance = Auto_Center - Self_Center
  )

upstream_count <- sum(filtered_df$Distance < 0)
downstream_count <- sum(filtered_df$Distance > 0)
equal_count <- sum(filtered_df$Distance == 0)

avg_abs_distance <- mean(abs(filtered_df$Distance))
avg_upstream_distance <- mean(abs(filtered_df$Distance[filtered_df$Distance < 0]))
avg_downstream_distance <- mean(abs(filtered_df$Distance[filtered_df$Distance > 0]))

cat("Autoimmune centers upstream of self center:", upstream_count, "\n")
cat("Autoimmune centers downstream of self center:", downstream_count, "\n")
cat("Autoimmune centers equal to self center:", equal_count, "\n")
cat("Average absolute distance from autoimmune to self center:", avg_abs_distance, "\n")
cat("Average upstream distance:", avg_upstream_distance, "\n")
cat("Average downstream distance:", avg_downstream_distance, "\n")