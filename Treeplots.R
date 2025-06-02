# ------------------------------------------------------------
# Load Libraries
# ------------------------------------------------------------
install.packages("treemapify")  # Only needed once
library(httr)
library(ggplot2)
library(dplyr)
library(treemapify)

###############################################################################
######################### Autoimmune Peptides #################################
###############################################################################
# ------------------------------------------------------------
# Amino Acid Property Definitions
# ------------------------------------------------------------
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

# ------------------------------------------------------------
# Helper Functions
# ------------------------------------------------------------
fetch_uniprot_sequence <- function(uniprot_id) {
  url <- paste0("https://rest.uniprot.org/uniprotkb/", uniprot_id, ".fasta")
  response <- GET(url)
  if (status_code(response) == 200) {
    fasta_text <- content(response, as = "text", encoding = "UTF-8")
    lines <- strsplit(fasta_text, "\n")[[1]]
    seq <- paste(lines[-1], collapse = "")
    return(seq)
  } else {
    return(NA)
  }
}

get_cleavage_pair <- function(peptide, full_seq, direction = "upstream") {
  match_pos <- regexpr(peptide, full_seq)[1]
  pep_len <- nchar(peptide)
  if (match_pos == -1) return(NA)
  
  if (direction == "upstream" && match_pos > 1) {
    return(substr(full_seq, match_pos - 1, match_pos))  # -1, +1 of N-term
  }
  
  if (direction == "downstream") {
    end_pos <- match_pos + pep_len - 1
    if (end_pos < nchar(full_seq)) {
      return(substr(full_seq, end_pos, end_pos + 1))  # -1, +1 of C-term
    }
  }
  return(NA)
}

get_overrepresented_list <- function(auto_vec, self_vec, alpha = 0.05) {
  auto_table <- table(auto_vec)
  self_table <- table(self_vec)
  shared <- intersect(names(auto_table), names(self_table))
  result_list <- list()
  
  for (d in shared) {
    auto_count <- auto_table[[d]]
    self_count <- self_table[[d]]
    matrix <- matrix(c(auto_count, length(auto_vec) - auto_count,
                       self_count, length(self_vec) - self_count),
                     nrow = 2, byrow = TRUE)
    test <- fisher.test(matrix, alternative = "greater")
    if (test$p.value < alpha) {
      result_list[[d]] <- auto_count
    }
  }
  
  df <- data.frame(Dipeptide = names(result_list), AutoCount = unlist(result_list), stringsAsFactors = FALSE)
  df <- df %>% arrange(Dipeptide)
  df$Index <- seq_len(nrow(df))
  return(df)
}

summarize_by_property <- function(aa_vector) {
  tbl <- table(aa_vector)
  df <- data.frame(AA = names(tbl), Count = as.vector(tbl))
  df <- merge(df, aa_properties, by = "AA")
  df %>% group_by(Property) %>% summarise(Count = sum(Count), .groups = "drop")
}

plot_treemap <- function(df, title, filename) {
  p <- ggplot(df, aes(area = Count, fill = Property, label = paste0(Property, "\n", Count))) +
    geom_treemap() +
    geom_treemap_text(colour = "white", grow = TRUE, reflow = TRUE) +
    labs(title = title) +
    theme(legend.position = "none")
  ggsave(filename, plot = p, width = 8, height = 6, dpi = 300)
}

# ------------------------------------------------------------
# Load Data
# ------------------------------------------------------------
auto_data <- read.csv("autoimmune_iedb_pepetides.csv", stringsAsFactors = FALSE)
self_data <- read.csv("Seed_Sample_Peptides.csv", stringsAsFactors = FALSE)

# ------------------------------------------------------------
# Extract Cleavage Pairs
# ------------------------------------------------------------
upstream_auto <- character()
downstream_auto <- character()
for (i in seq_len(nrow(auto_data))) {
  pep <- auto_data[[1]][i]
  uid <- auto_data[[4]][i]
  if (!is.na(pep) && !is.na(uid) && !grepl("^NP_", uid)) {
    seq <- fetch_uniprot_sequence(uid)
    if (!is.na(seq)) {
      up <- get_cleavage_pair(pep, seq, "upstream")
      down <- get_cleavage_pair(pep, seq, "downstream")
      if (!is.na(up)) upstream_auto <- c(upstream_auto, up)
      if (!is.na(down)) downstream_auto <- c(downstream_auto, down)
    }
  }
}

upstream_self <- character()
downstream_self <- character()
for (i in seq_len(nrow(self_data))) {
  pep <- self_data[[2]][i]
  uid <- self_data[[3]][i]
  if (!is.na(pep) && !is.na(uid)) {
    uid_clean <- sub("\\..*", "", uid)
    seq <- fetch_uniprot_sequence(uid_clean)
    if (!is.na(seq)) {
      up <- get_cleavage_pair(pep, seq, "upstream")
      down <- get_cleavage_pair(pep, seq, "downstream")
      if (!is.na(up)) upstream_self <- c(upstream_self, up)
      if (!is.na(down)) downstream_self <- c(downstream_self, down)
    }
  }
}

# ------------------------------------------------------------
# Identify Exclusive and Overrepresented Dipeptides
# ------------------------------------------------------------
exclusive_upstream <- setdiff(unique(upstream_auto), unique(upstream_self))
exclusive_downstream <- setdiff(unique(downstream_auto), unique(downstream_self))

sig_upstream_df <- get_overrepresented_list(upstream_auto, upstream_self)
sig_downstream_df <- get_overrepresented_list(downstream_auto, downstream_self)

combined_upstream <- unique(c(sig_upstream_df$Dipeptide, exclusive_upstream))
combined_downstream <- unique(c(sig_downstream_df$Dipeptide, exclusive_downstream))

# ------------------------------------------------------------
# Extract AA Positions from Cleavage Pairs
# ------------------------------------------------------------
aa_minus1_nterm <- substr(combined_upstream, 1, 1)  # residue before peptide starts
aa_plus1_nterm  <- substr(combined_upstream, 2, 2)  # first residue in peptide
aa_minus1_cterm <- substr(combined_downstream, 1, 1)  # last residue in peptide
aa_plus1_cterm  <- substr(combined_downstream, 2, 2)  # residue after peptide

# ------------------------------------------------------------
# Summarize by Biochemical Properties
# ------------------------------------------------------------
df_minus1_nterm_summary <- summarize_by_property(aa_minus1_nterm)
df_plus1_nterm_summary  <- summarize_by_property(aa_plus1_nterm)
df_minus1_cterm_summary <- summarize_by_property(aa_minus1_cterm)
df_plus1_cterm_summary  <- summarize_by_property(aa_plus1_cterm)

# ------------------------------------------------------------
# Generate Treemaps
# ------------------------------------------------------------
plot_treemap(df_minus1_nterm_summary, "Cleavage Site (-1): Before N-terminal (Autoimmune-specific)", "Cleavage_Minus1_Nterm_Treemap.png")
plot_treemap(df_plus1_nterm_summary,  "Cleavage Site (+1): First AA in Peptide (Autoimmune-specific)", "Cleavage_Plus1_Nterm_Treemap.png")
plot_treemap(df_minus1_cterm_summary, "Cleavage Site (-1): Last AA in Peptide (Autoimmune-specific)", "Cleavage_Minus1_Cterm_Treemap.png")
plot_treemap(df_plus1_cterm_summary,  "Cleavage Site (+1): After C-terminal (Autoimmune-specific)", "Cleavage_Plus1_Cterm_Treemap.png")

###############################################################################
############################ Self-Peptides ####################################
###############################################################################
# Invert the logic to detect self-specific dipeptides
exclusive_upstream_self <- setdiff(unique(upstream_self), unique(upstream_auto))
exclusive_downstream_self <- setdiff(unique(downstream_self), unique(downstream_auto))

# Use inverted Fisher's test for overrepresentation in self
get_self_overrepresented <- function(self_vec, auto_vec, alpha = 0.05) {
  self_table <- table(self_vec)
  auto_table <- table(auto_vec)
  shared <- intersect(names(self_table), names(auto_table))
  result_list <- list()
  
  for (d in shared) {
    self_count <- self_table[[d]]
    auto_count <- auto_table[[d]]
    matrix <- matrix(c(self_count, length(self_vec) - self_count,
                       auto_count, length(auto_vec) - auto_count),
                     nrow = 2, byrow = TRUE)
    test <- fisher.test(matrix, alternative = "greater")
    if (test$p.value < alpha) {
      result_list[[d]] <- self_count
    }
  }
  
  df <- data.frame(Dipeptide = names(result_list), SelfCount = unlist(result_list), stringsAsFactors = FALSE)
  df <- df %>% arrange(Dipeptide)
  df$Index <- seq_len(nrow(df))
  return(df)
}

sig_upstream_self_df <- get_self_overrepresented(upstream_self, upstream_auto)
sig_downstream_self_df <- get_self_overrepresented(downstream_self, downstream_auto)

combined_upstream_self <- unique(c(sig_upstream_self_df$Dipeptide, exclusive_upstream_self))
combined_downstream_self <- unique(c(sig_downstream_self_df$Dipeptide, exclusive_downstream_self))

# ------------------------------------------------------------
# Extract AA positions from cleavage dipeptides
# ------------------------------------------------------------
aa_minus1_nterm_self <- substr(combined_upstream_self, 1, 1)
aa_plus1_nterm_self  <- substr(combined_upstream_self, 2, 2)
aa_minus1_cterm_self <- substr(combined_downstream_self, 1, 1)
aa_plus1_cterm_self  <- substr(combined_downstream_self, 2, 2)

# ------------------------------------------------------------
# Summarize by property
# ------------------------------------------------------------
df_minus1_nterm_self_summary <- summarize_by_property(aa_minus1_nterm_self)
df_plus1_nterm_self_summary  <- summarize_by_property(aa_plus1_nterm_self)
df_minus1_cterm_self_summary <- summarize_by_property(aa_minus1_cterm_self)
df_plus1_cterm_self_summary  <- summarize_by_property(aa_plus1_cterm_self)

# ------------------------------------------------------------
# Plot and save treemaps
# ------------------------------------------------------------
plot_treemap(df_minus1_nterm_self_summary, "Cleavage Site (-1): Before N-terminal (Self-specific)", "Cleavage_Minus1_Nterm_Treemap_Self.png")
plot_treemap(df_plus1_nterm_self_summary,  "Cleavage Site (+1): First AA in Peptide (Self-specific)", "Cleavage_Plus1_Nterm_Treemap_Self.png")
plot_treemap(df_minus1_cterm_self_summary, "Cleavage Site (-1): Last AA in Peptide (Self-specific)", "Cleavage_Minus1_Cterm_Treemap_Self.png")
plot_treemap(df_plus1_cterm_self_summary,  "Cleavage Site (+1): After C-terminal (Self-specific)", "Cleavage_Plus1_Cterm_Treemap_Self.png")



###############################################################################
###############################################################################
######################### W/O BACKGROUND BIASES ###############################
###############################################################################
###############################################################################

# ------------------------------------------------------------
# 1. Calculate Background Amino Acid Frequencies
# ------------------------------------------------------------

# Function to calculate background amino acid frequencies
calculate_background_freq <- function(peptide_vector) {
  all_aa <- unlist(strsplit(peptide_vector, split = ""))
  total_count <- length(all_aa)
  freq_table <- table(all_aa) / total_count  # Fractional frequency
  return(freq_table)
}

# Calculate for autoimmune and self datasets
background_auto_freq <- calculate_background_freq(auto_data[[1]])   # Autoimmune peptides
background_self_freq <- calculate_background_freq(self_data[[2]])   # Self peptides

# ------------------------------------------------------------
# 2. Normalize Cleavage Site Counts by Background
# ------------------------------------------------------------

# Updated summarize_by_property function
summarize_normalized_by_property <- function(aa_vector, background_freq) {
  aa_table <- table(aa_vector)
  
  # Observed frequencies
  observed_freq <- aa_table / sum(aa_table)
  
  # Merge into data frame
  df <- data.frame(AA = names(observed_freq), ObservedFreq = as.numeric(observed_freq))
  
  # Add background frequency
  df$BackgroundFreq <- background_freq[df$AA]
  
  # Handle cases where AA was not seen in background (rare but possible)
  df$BackgroundFreq[is.na(df$BackgroundFreq)] <- min(df$BackgroundFreq, na.rm = TRUE)
  
  # Calculate enrichment ratio
  df$Enrichment <- df$ObservedFreq / df$BackgroundFreq
  
  # Map properties
  df <- merge(df, aa_properties, by = "AA")
  
  # Average enrichment by property
  summary_df <- df %>%
    group_by(Property) %>%
    summarise(AverageEnrichment = mean(Enrichment), .groups = "drop")
  
  return(summary_df)
}

# ------------------------------------------------------------
# 3. Generate Normalized Summaries
# ------------------------------------------------------------

# For autoimmune-specific cleavage points (already defined: combined_upstream, combined_downstream)

# Extract positions from autoimmune-specific cleavages
aa_minus1_nterm_auto <- substr(combined_upstream, 1, 1)
aa_plus1_nterm_auto  <- substr(combined_upstream, 2, 2)
aa_minus1_cterm_auto <- substr(combined_downstream, 1, 1)
aa_plus1_cterm_auto  <- substr(combined_downstream, 2, 2)

# Summarize with normalization
df_minus1_nterm_auto_summary <- summarize_normalized_by_property(aa_minus1_nterm_auto, background_auto_freq)
df_plus1_nterm_auto_summary  <- summarize_normalized_by_property(aa_plus1_nterm_auto, background_auto_freq)
df_minus1_cterm_auto_summary <- summarize_normalized_by_property(aa_minus1_cterm_auto, background_auto_freq)
df_plus1_cterm_auto_summary  <- summarize_normalized_by_property(aa_plus1_cterm_auto, background_auto_freq)

# ------------------------------------------------------------
# 4. Plot Normalized Treemaps
# ------------------------------------------------------------

plot_normalized_treemap <- function(df, title, filename) {
  p <- ggplot(df, aes(area = AverageEnrichment, fill = Property, label = paste0(Property, "\n", round(AverageEnrichment, 2)))) +
    geom_treemap() +
    geom_treemap_text(colour = "white", grow = TRUE, reflow = TRUE) +
    labs(title = title) +
    theme(legend.position = "none")
  
  ggsave(filename, plot = p, width = 8, height = 6, dpi = 300)
}

# Autoimmune normalized treemaps
plot_normalized_treemap(df_minus1_nterm_auto_summary, "Autoimmune: -1 N-terminal (Bias Corrected)", "Norm_Cleavage_Minus1_Nterm_Auto.png")
plot_normalized_treemap(df_plus1_nterm_auto_summary,  "Autoimmune: +1 N-terminal (Bias Corrected)", "Norm_Cleavage_Plus1_Nterm_Auto.png")
plot_normalized_treemap(df_minus1_cterm_auto_summary, "Autoimmune: -1 C-terminal (Bias Corrected)", "Norm_Cleavage_Minus1_Cterm_Auto.png")
plot_normalized_treemap(df_plus1_cterm_auto_summary,  "Autoimmune: +1 C-terminal (Bias Corrected)", "Norm_Cleavage_Plus1_Cterm_Auto.png")

# ------------------------------------------------------------
# 5. Repeat for Self-specific Cleavage Points
# ------------------------------------------------------------

# Self-specific combined dipeptides (already defined: combined_upstream_self, combined_downstream_self)

# Extract positions
aa_minus1_nterm_self <- substr(combined_upstream_self, 1, 1)
aa_plus1_nterm_self  <- substr(combined_upstream_self, 2, 2)
aa_minus1_cterm_self <- substr(combined_downstream_self, 1, 1)
aa_plus1_cterm_self  <- substr(combined_downstream_self, 2, 2)

# Summarize with normalization
df_minus1_nterm_self_summary <- summarize_normalized_by_property(aa_minus1_nterm_self, background_self_freq)
df_plus1_nterm_self_summary  <- summarize_normalized_by_property(aa_plus1_nterm_self, background_self_freq)
df_minus1_cterm_self_summary <- summarize_normalized_by_property(aa_minus1_cterm_self, background_self_freq)
df_plus1_cterm_self_summary  <- summarize_normalized_by_property(aa_plus1_cterm_self, background_self_freq)

# Self normalized treemaps
plot_normalized_treemap(df_minus1_nterm_self_summary, "Self: -1 N-terminal (Bias Corrected)", "Norm_Cleavage_Minus1_Nterm_Self.png")
plot_normalized_treemap(df_plus1_nterm_self_summary,  "Self: +1 N-terminal (Bias Corrected)", "Norm_Cleavage_Plus1_Nterm_Self.png")
plot_normalized_treemap(df_minus1_cterm_self_summary, "Self: -1 C-terminal (Bias Corrected)", "Norm_Cleavage_Minus1_Cterm_Self.png")
plot_normalized_treemap(df_plus1_cterm_self_summary,  "Self: +1 C-terminal (Bias Corrected)", "Norm_Cleavage_Plus1_Cterm_Self.png")



###############################################################################
###############################################################################
######################## Comparison of Averages ###############################
###############################################################################
###############################################################################

# ------------------------------------------------------------
# Create Pie Chart for Background Amino Acid Frequencies
# ------------------------------------------------------------

# Convert background frequencies into data frames
df_auto_bg <- data.frame(
  AA = names(background_auto_freq),
  Frequency = as.numeric(background_auto_freq)
)

df_self_bg <- data.frame(
  AA = names(background_self_freq),
  Frequency = as.numeric(background_self_freq)
)

# Function to plot pie chart
plot_piechart <- function(df, title, filename) {
  p <- ggplot(df, aes(x = "", y = Frequency, fill = AA)) +
    geom_col(width = 1) +
    coord_polar(theta = "y") +
    labs(title = title, x = NULL, y = NULL) +
    theme_void() +
    theme(legend.position = "right")
  
  ggsave(filename, plot = p, width = 8, height = 6, dpi = 300)
}

# Generate and save pie charts
plot_piechart(df_auto_bg, "Amino Acid Composition - Autoimmune Peptides", "AminoAcid_Pie_Autoimmune.png")
plot_piechart(df_self_bg, "Amino Acid Composition - Self Peptides", "AminoAcid_Pie_Self.png")

# ------------------------------------------------------------
# Create Bar Plot for Comparison of Background Amino Acid Frequencies
# ------------------------------------------------------------

# Combine for side-by-side barplot
df_auto_bg$Group <- "Autoimmune"
df_self_bg$Group <- "Self"
df_combined <- rbind(df_auto_bg, df_self_bg)

# Plot side-by-side barplot
p_bar <- ggplot(df_combined, aes(x = AA, y = Frequency, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Amino Acid Composition Comparison", x = "Amino Acid", y = "Frequency") +
  theme_minimal()

ggsave("AminoAcid_Barplot_Comparison.png", plot = p_bar, width = 10, height = 6, dpi = 300)


###############################################################################
###############################################################################
######################## Statistical Test Avg  ################################
###############################################################################
###############################################################################

# ------------------------------------------------------------
# 6. Amino Acid Specific Statistical Tests (Fisher's Exact Test per AA)
# ------------------------------------------------------------

# Again, calculate total counts
total_aa_auto <- sum(nchar(auto_data[[1]]))  # Total amino acids autoimmune
total_aa_self <- sum(nchar(self_data[[2]]))  # Total amino acids self

# Raw counts
aa_auto_counts <- background_auto_freq * total_aa_auto
aa_self_counts <- background_self_freq * total_aa_self

# Align amino acids
common_aas <- intersect(names(aa_auto_counts), names(aa_self_counts))
aa_auto_counts <- aa_auto_counts[common_aas]
aa_self_counts <- aa_self_counts[common_aas]

# Initialize results list
aa_test_results <- data.frame(
  AminoAcid = character(),
  PValue = numeric(),
  Autoimmune_Fraction = numeric(),
  Self_Fraction = numeric(),
  stringsAsFactors = FALSE
)

# Run Fisher's Exact Test for each amino acid
for (aa in common_aas) {
  # Build 2x2 table
  matrix_2x2 <- matrix(c(
    aa_auto_counts[aa], total_aa_auto - aa_auto_counts[aa],
    aa_self_counts[aa], total_aa_self - aa_self_counts[aa]
  ), nrow = 2, byrow = TRUE)
  
  fisher_res <- fisher.test(matrix_2x2)
  
  aa_test_results <- rbind(aa_test_results, data.frame(
    AminoAcid = aa,
    PValue = fisher_res$p.value,
    Autoimmune_Fraction = aa_auto_counts[aa] / total_aa_auto,
    Self_Fraction = aa_self_counts[aa] / total_aa_self,
    stringsAsFactors = FALSE
  ))
}

# Optional: Adjust p-values for multiple testing (Benjamini-Hochberg FDR correction)
aa_test_results$AdjustedPValue <- p.adjust(aa_test_results$PValue, method = "BH")

# Sort by adjusted p-value
aa_test_results <- aa_test_results %>% arrange(AdjustedPValue)

# View top significant amino acids
print(aa_test_results)

# Save to file
write.csv(aa_test_results, "Fisher_AA_Composition_Differences.csv", row.names = FALSE)

# ------------------------------------------------------------
# Volcano Plot: Amino Acid Differences (Autoimmune vs Self)
# ------------------------------------------------------------

# Calculate log2 fold-change
aa_test_results$Log2FoldChange <- log2(aa_test_results$Autoimmune_Fraction / aa_test_results$Self_Fraction)

# Calculate -log10 adjusted p-value
aa_test_results$NegLog10Pval <- -log10(aa_test_results$AdjustedPValue)

# Classify points for coloring
aa_test_results$Category <- "Not Significant"
aa_test_results$Category[aa_test_results$AdjustedPValue < 0.05 & aa_test_results$Log2FoldChange > 0] <- "Enriched in Autoimmune"
aa_test_results$Category[aa_test_results$AdjustedPValue < 0.05 & aa_test_results$Log2FoldChange < 0] <- "Depleted in Autoimmune"

# Plot
p_volcano <- ggplot(aa_test_results, aes(x = Log2FoldChange, y = NegLog10Pval, label = AminoAcid)) +
  geom_point(aes(color = Category), size = 3) +
  scale_color_manual(values = c(
    "Enriched in Autoimmune" = "red",
    "Depleted in Autoimmune" = "blue",
    "Not Significant" = "grey"
  )) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(title = "Volcano Plot: Amino Acid Differences (Autoimmune vs Self)",
       x = "Log2 Fold Change (Autoimmune vs Self)",
       y = "-Log10 Adjusted P-value") +
  theme_minimal()

# Save the volcano plot
ggsave("VolcanoPlot_AA_Comparison.png", plot = p_volcano, width = 10, height = 7, dpi = 300)


###############################################################################
###############################################################################
########################## Combine Backgrounds ################################
###############################################################################
###############################################################################

# ------------------------------------------------------------
# 7. New: Combined Background Frequencies
# ------------------------------------------------------------

# Calculate combined counts
total_aa_auto <- sum(nchar(auto_data[[1]]))
total_aa_self <- sum(nchar(self_data[[2]]))

# Raw counts
aa_auto_counts <- background_auto_freq * total_aa_auto
aa_self_counts <- background_self_freq * total_aa_self

# Combine counts
combined_counts <- aa_auto_counts + aa_self_counts
combined_total <- total_aa_auto + total_aa_self

# Combined background frequency
background_combined_freq <- combined_counts / combined_total

# ------------------------------------------------------------
# 8. New: Re-normalize Cleavage Sites Against Combined Background
# ------------------------------------------------------------

# Autoimmune-specific cleavage sites (already extracted earlier)
aa_minus1_nterm_auto <- substr(combined_upstream, 1, 1)
aa_plus1_nterm_auto  <- substr(combined_upstream, 2, 2)
aa_minus1_cterm_auto <- substr(combined_downstream, 1, 1)
aa_plus1_cterm_auto  <- substr(combined_downstream, 2, 2)

# Self-specific cleavage sites
aa_minus1_nterm_self <- substr(combined_upstream_self, 1, 1)
aa_plus1_nterm_self  <- substr(combined_upstream_self, 2, 2)
aa_minus1_cterm_self <- substr(combined_downstream_self, 1, 1)
aa_plus1_cterm_self  <- substr(combined_downstream_self, 2, 2)

# Summarize autoimmune vs self with combined background
df_minus1_nterm_auto_combined <- summarize_normalized_by_property(aa_minus1_nterm_auto, background_combined_freq)
df_plus1_nterm_auto_combined  <- summarize_normalized_by_property(aa_plus1_nterm_auto, background_combined_freq)
df_minus1_cterm_auto_combined <- summarize_normalized_by_property(aa_minus1_cterm_auto, background_combined_freq)
df_plus1_cterm_auto_combined  <- summarize_normalized_by_property(aa_plus1_cterm_auto, background_combined_freq)

df_minus1_nterm_self_combined <- summarize_normalized_by_property(aa_minus1_nterm_self, background_combined_freq)
df_plus1_nterm_self_combined  <- summarize_normalized_by_property(aa_plus1_nterm_self, background_combined_freq)
df_minus1_cterm_self_combined <- summarize_normalized_by_property(aa_minus1_cterm_self, background_combined_freq)
df_plus1_cterm_self_combined  <- summarize_normalized_by_property(aa_plus1_cterm_self, background_combined_freq)

# ------------------------------------------------------------
# 9. New: Treemaps normalized against Combined Background
# ------------------------------------------------------------

# Autoimmune
plot_normalized_treemap(df_minus1_nterm_auto_combined, "Autoimmune: -1 N-terminal (Combined Background)", "Norm_Combined_Cleavage_Minus1_Nterm_Auto.png")
plot_normalized_treemap(df_plus1_nterm_auto_combined,  "Autoimmune: +1 N-terminal (Combined Background)", "Norm_Combined_Cleavage_Plus1_Nterm_Auto.png")
plot_normalized_treemap(df_minus1_cterm_auto_combined, "Autoimmune: -1 C-terminal (Combined Background)", "Norm_Combined_Cleavage_Minus1_Cterm_Auto.png")
plot_normalized_treemap(df_plus1_cterm_auto_combined,  "Autoimmune: +1 C-terminal (Combined Background)", "Norm_Combined_Cleavage_Plus1_Cterm_Auto.png")

# Self
plot_normalized_treemap(df_minus1_nterm_self_combined, "Self: -1 N-terminal (Combined Background)", "Norm_Combined_Cleavage_Minus1_Nterm_Self.png")
plot_normalized_treemap(df_plus1_nterm_self_combined,  "Self: +1 N-terminal (Combined Background)", "Norm_Combined_Cleavage_Plus1_Nterm_Self.png")
plot_normalized_treemap(df_minus1_cterm_self_combined, "Self: -1 C-terminal (Combined Background)", "Norm_Combined_Cleavage_Minus1_Cterm_Self.png")
plot_normalized_treemap(df_plus1_cterm_self_combined,  "Self: +1 C-terminal (Combined Background)", "Norm_Combined_Cleavage_Plus1_Cterm_Self.png")
















###############################################################################
###############################################################################
############################ New Proposition ##################################
###############################################################################
###############################################################################

# ------------------------------------------------------------
# Load Libraries
# ------------------------------------------------------------
library(httr)
library(ggplot2)
library(dplyr)
library(treemapify)

# ------------------------------------------------------------
# Amino Acid Property Definitions
# ------------------------------------------------------------
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

# ------------------------------------------------------------
# Helper Functions
# ------------------------------------------------------------
fetch_uniprot_sequence <- function(uniprot_id) {
  url <- paste0("https://rest.uniprot.org/uniprotkb/", uniprot_id, ".fasta")
  response <- GET(url)
  if (status_code(response) == 200) {
    fasta_text <- content(response, as = "text", encoding = "UTF-8")
    lines <- strsplit(fasta_text, "\n")[[1]]
    seq <- paste(lines[-1], collapse = "")
    return(seq)
  } else {
    return(NA)
  }
}

remove_peptides <- function(protein_seq, peptides) {
  for (pep in peptides) {
    protein_seq <- gsub(pep, "", protein_seq, fixed = TRUE)
  }
  return(protein_seq)
}

calculate_background_freq <- function(sequences) {
  all_aa <- unlist(strsplit(sequences, split = ""))
  table(all_aa) / length(all_aa)
}

summarize_normalized_by_property <- function(aa_vector, background_freq) {
  aa_table <- table(aa_vector)
  observed_freq <- aa_table / sum(aa_table)
  df <- data.frame(AA = names(observed_freq), ObservedFreq = as.numeric(observed_freq))
  df$BackgroundFreq <- background_freq[df$AA]
  df$BackgroundFreq[is.na(df$BackgroundFreq)] <- min(df$BackgroundFreq, na.rm = TRUE)
  df$Enrichment <- df$ObservedFreq / df$BackgroundFreq
  df <- merge(df, aa_properties, by = "AA")
  df %>%
    group_by(Property) %>%
    summarise(AverageEnrichment = mean(Enrichment), .groups = "drop")
}

plot_normalized_treemap <- function(df, title, filename) {
  p <- ggplot(df, aes(area = AverageEnrichment, fill = Property, label = paste0(Property, "\n", round(AverageEnrichment, 2)))) +
    geom_treemap() +
    geom_treemap_text(colour = "white", grow = TRUE, reflow = TRUE) +
    labs(title = title) +
    theme(legend.position = "none")
  ggsave(filename, plot = p, width = 8, height = 6, dpi = 300)
}

# ------------------------------------------------------------
# Load Data
# ------------------------------------------------------------
auto_data <- read.csv("autoimmune_iedb_pepetides.csv", stringsAsFactors = FALSE)
self_data <- read.csv("Seed_Sample_Peptides.csv", stringsAsFactors = FALSE)

# ------------------------------------------------------------
# Process Each Group
# ------------------------------------------------------------
process_group <- function(data, peptide_col, uid_col) {
  peptides <- character()
  peptides_minus <- character()
  peptide_positions <- list()
  proteins <- list()
  proteins_cleaned <- list()
  
  for (i in seq_len(nrow(data))) {
    pep <- data[[peptide_col]][i]
    uid <- data[[uid_col]][i]
    uid_clean <- sub("\\..*", "", uid)
    if (!is.na(pep) && !is.na(uid_clean) && !grepl("^NP_", uid_clean)) {
      seq <- fetch_uniprot_sequence(uid_clean)
      if (!is.na(seq)) {
        pos <- regexpr(pep, seq, fixed = TRUE)[1]
        if (pos != -1) {
          peptides <- c(peptides, pep)
          proteins[[length(proteins)+1]] <- seq
          peptide_positions[[length(peptide_positions)+1]] <- list(seq=seq, pos=pos, len=nchar(pep))
          proteins_cleaned[[length(proteins_cleaned)+1]] <- gsub(pep, "", seq, fixed = TRUE)
        }
      }
    }
  }
  
  # Compute AA frequencies
  peptide_freq <- calculate_background_freq(peptides)
  protein_minus_freq <- calculate_background_freq(unlist(proteins_cleaned))
  
  # Extract cleavage AAs
  aa_minus1_nterm <- aa_plus1_nterm <- aa_minus1_cterm <- aa_plus1_cterm <- character()
  
  for (info in peptide_positions) {
    seq <- info$seq
    pos <- info$pos
    len <- info$len
    end_pos <- pos + len - 1
    
    if (pos > 1) aa_minus1_nterm <- c(aa_minus1_nterm, substr(seq, pos - 1, pos - 1))
    aa_plus1_nterm <- c(aa_plus1_nterm, substr(seq, pos, pos))
    aa_minus1_cterm <- c(aa_minus1_cterm, substr(seq, end_pos, end_pos))
    if (end_pos < nchar(seq)) aa_plus1_cterm <- c(aa_plus1_cterm, substr(seq, end_pos + 1, end_pos + 1))
  }
  
  return(list(
    peptide_freq = peptide_freq,
    protein_freq = protein_minus_freq,
    aa = list(
      minus1_nterm = aa_minus1_nterm,
      plus1_nterm = aa_plus1_nterm,
      minus1_cterm = aa_minus1_cterm,
      plus1_cterm = aa_plus1_cterm
    )
  ))
}

auto_group <- process_group(auto_data, 1, 4)
self_group <- process_group(self_data, 2, 3)

# ------------------------------------------------------------
# Generate and Save Treemaps
# ------------------------------------------------------------
generate_all_treemaps <- function(group, group_name) {
  # -1 N-term and +1 C-term use protein background
  df1 <- summarize_normalized_by_property(group$aa$minus1_nterm, group$protein_freq)
  df2 <- summarize_normalized_by_property(group$aa$plus1_cterm,  group$protein_freq)
  # +1 N-term and -1 C-term use peptide background
  df3 <- summarize_normalized_by_property(group$aa$plus1_nterm,  group$peptide_freq)
  df4 <- summarize_normalized_by_property(group$aa$minus1_cterm, group$peptide_freq)
  
  plot_normalized_treemap(df1, paste0(group_name, ": -1 N-term"), paste0(group_name, "_Minus1_Nterm.png"))
  plot_normalized_treemap(df2, paste0(group_name, ": +1 C-term"), paste0(group_name, "_Plus1_Cterm.png"))
  plot_normalized_treemap(df3, paste0(group_name, ": +1 N-term"), paste0(group_name, "_Plus1_Nterm.png"))
  plot_normalized_treemap(df4, paste0(group_name, ": -1 C-term"), paste0(group_name, "_Minus1_Cterm.png"))
}

generate_all_treemaps(auto_group, "Autoimmune")
generate_all_treemaps(self_group, "Self")


###############################################################################
###############################################################################
########################## -5s and 5's extended ###############################
###############################################################################
###############################################################################



# ------------------------------------------------------------
# Load Libraries
# ------------------------------------------------------------
library(httr)
library(ggplot2)
library(dplyr)
library(treemapify)

# ------------------------------------------------------------
# Amino Acid Properties
# ------------------------------------------------------------
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

# ------------------------------------------------------------
# Helper Functions
# ------------------------------------------------------------
fetch_uniprot_sequence <- function(uniprot_id) {
  url <- paste0("https://rest.uniprot.org/uniprotkb/", uniprot_id, ".fasta")
  response <- GET(url)
  if (status_code(response) == 200) {
    fasta_text <- content(response, as = "text", encoding = "UTF-8")
    lines <- strsplit(fasta_text, "\n")[[1]]
    seq <- paste(lines[-1], collapse = "")
    return(seq)
  } else {
    return(NA)
  }
}

calculate_background_freq <- function(sequences) {
  all_aa <- unlist(strsplit(sequences, split = ""))
  table(all_aa) / length(all_aa)
}

summarize_normalized_by_property <- function(aa_vector, background_freq) {
  aa_table <- table(aa_vector)
  observed_freq <- aa_table / sum(aa_table)
  df <- data.frame(AA = names(observed_freq), ObservedFreq = as.numeric(observed_freq))
  df$BackgroundFreq <- background_freq[df$AA]
  df$BackgroundFreq[is.na(df$BackgroundFreq)] <- min(df$BackgroundFreq, na.rm = TRUE)
  df$Enrichment <- df$ObservedFreq / df$BackgroundFreq
  df <- merge(df, aa_properties, by = "AA")
  df %>%
    group_by(Property) %>%
    summarise(AverageEnrichment = mean(Enrichment), .groups = "drop")
}

plot_normalized_treemap <- function(df, title, filename) {
  p <- ggplot(df, aes(area = AverageEnrichment, fill = Property, label = paste0(Property, "\n", round(AverageEnrichment, 2)))) +
    geom_treemap() +
    geom_treemap_text(colour = "white", grow = TRUE, reflow = TRUE) +
    labs(title = title) +
    theme(legend.position = "none")
  ggsave(filename, plot = p, width = 8, height = 6, dpi = 300)
}

# ------------------------------------------------------------
# Extended Cleavage Site Window Logic
# ------------------------------------------------------------
extract_extended_context <- function(seq, pos, len, direction = "N") {
  if (direction == "N") {
    start <- max(1, pos - 5)
    end <- pos - 1
  } else if (direction == "N+5") {
    start <- pos
    end <- min(nchar(seq), pos + 4)
  } else if (direction == "C-5") {
    start <- max(1, pos + len - 5)
    end <- pos + len - 1
  } else if (direction == "C") {
    start <- pos + len
    end <- min(nchar(seq), pos + len + 4)
  } else {
    return(NA)
  }
  if (start <= end && end <= nchar(seq)) {
    return(strsplit(substr(seq, start, end), "")[[1]])
  } else {
    return(NA)
  }
}

process_group_extended <- function(data, peptide_col, uid_col) {
  peptides <- character()
  proteins <- list()
  proteins_cleaned <- list()
  aa_minus5_to_minus1_nterm <- list()
  aa_plus1_to_plus5_nterm <- list()
  aa_minus5_to_minus1_cterm <- list()
  aa_plus1_to_plus5_cterm <- list()
  for (i in seq_len(nrow(data))) {
    pep <- data[[peptide_col]][i]
    uid <- data[[uid_col]][i]
    uid_clean <- sub("\\..*", "", uid)
    if (!is.na(pep) && !is.na(uid_clean) && !grepl("^NP_", uid_clean)) {
      seq <- fetch_uniprot_sequence(uid_clean)
      if (!is.na(seq)) {
        pos <- regexpr(pep, seq, fixed = TRUE)[1]
        if (pos != -1) {
          peptides <- c(peptides, pep)
          proteins[[length(proteins)+1]] <- seq
          proteins_cleaned[[length(proteins_cleaned)+1]] <- gsub(pep, "", seq, fixed = TRUE)
          aa_minus5_to_minus1_nterm[[length(aa_minus5_to_minus1_nterm)+1]] <- extract_extended_context(seq, pos, nchar(pep), "N")
          aa_plus1_to_plus5_nterm[[length(aa_plus1_to_plus5_nterm)+1]] <- extract_extended_context(seq, pos, nchar(pep), "N+5")
          aa_minus5_to_minus1_cterm[[length(aa_minus5_to_minus1_cterm)+1]] <- extract_extended_context(seq, pos, nchar(pep), "C-5")
          aa_plus1_to_plus5_cterm[[length(aa_plus1_to_plus5_cterm)+1]] <- extract_extended_context(seq, pos, nchar(pep), "C")
        }
      }
    }
  }
  context <- list(
    minus5_to_minus1_nterm = unlist(aa_minus5_to_minus1_nterm),
    plus1_to_plus5_nterm = unlist(aa_plus1_to_plus5_nterm),
    minus5_to_minus1_cterm = unlist(aa_minus5_to_minus1_cterm),
    plus1_to_plus5_cterm = unlist(aa_plus1_to_plus5_cterm)
  )
  return(list(
    peptide_freq = calculate_background_freq(peptides),
    protein_freq = calculate_background_freq(unlist(proteins_cleaned)),
    aa = context
  ))
}

# ------------------------------------------------------------
# Load Data and Run
# ------------------------------------------------------------
auto_data <- read.csv("autoimmune_iedb_pepetides.csv", stringsAsFactors = FALSE)
self_data <- read.csv("Seed_Sample_Peptides.csv", stringsAsFactors = FALSE)
auto_group <- process_group_extended(auto_data, 1, 4)
self_group <- process_group_extended(self_data, 2, 3)

# ------------------------------------------------------------
# Plot Treemaps for Autoimmune
# ------------------------------------------------------------
plot_normalized_treemap(summarize_normalized_by_property(auto_group$aa$minus5_to_minus1_nterm, auto_group$protein_freq),
                        "Autoimmune: -5 to -1 N-term", "Auto_Ext_Minus5_Nterm.png")
plot_normalized_treemap(summarize_normalized_by_property(auto_group$aa$plus1_to_plus5_nterm, auto_group$peptide_freq),
                        "Autoimmune: +1 to +5 N-term", "Auto_Ext_Plus5_Nterm.png")
plot_normalized_treemap(summarize_normalized_by_property(auto_group$aa$minus5_to_minus1_cterm, auto_group$peptide_freq),
                        "Autoimmune: -5 to -1 C-term", "Auto_Ext_Minus5_Cterm.png")
plot_normalized_treemap(summarize_normalized_by_property(auto_group$aa$plus1_to_plus5_cterm, auto_group$protein_freq),
                        "Autoimmune: +1 to +5 C-term", "Auto_Ext_Plus5_Cterm.png")

# ------------------------------------------------------------
# Plot Treemaps for Self
# ------------------------------------------------------------
plot_normalized_treemap(summarize_normalized_by_property(self_group$aa$minus5_to_minus1_nterm, self_group$protein_freq),
                        "Self: -5 to -1 N-term", "Self_Ext_Minus5_Nterm.png")
plot_normalized_treemap(summarize_normalized_by_property(self_group$aa$plus1_to_plus5_nterm, self_group$peptide_freq),
                        "Self: +1 to +5 N-term", "Self_Ext_Plus5_Nterm.png")
plot_normalized_treemap(summarize_normalized_by_property(self_group$aa$minus5_to_minus1_cterm, self_group$peptide_freq),
                        "Self: -5 to -1 C-term", "Self_Ext_Minus5_Cterm.png")
plot_normalized_treemap(summarize_normalized_by_property(self_group$aa$plus1_to_plus5_cterm, self_group$protein_freq),
                        "Self: +1 to +5 C-term", "Self_Ext_Plus5_Cterm.png")

