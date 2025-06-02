########### Differences between merged self-peptide and autoimmune #############
# Load data
data <- read.csv("Matched_Core_Peptides_Cleaned.csv", stringsAsFactors = FALSE)

# Use only the columns we need
pairs <- unique(data[, c(1, 7)])
colnames(pairs) <- c("Autoimmune_Peptide", "Self_Peptide")

# Function to find longest common substring
find_longest_common_substring <- function(a, b) {
  max_len <- 0
  best_match <- ""
  
  for (i in 1:nchar(a)) {
    for (j in i:nchar(a)) {
      substr_a <- substr(a, i, j)
      if (nchar(substr_a) >= 9 && grepl(substr_a, b) && nchar(substr_a) > max_len) {
        max_len <- nchar(substr_a)
        best_match <- substr_a
      }
    }
  }
  return(best_match)
}

# Initialize output
output <- data.frame(
  Pair_ID = integer(),
  Autoimmune_Peptide = character(),
  Self_Peptide = character(),
  Shared_Overlap = character(),
  Auto_Prefix = character(),
  Auto_Suffix = character(),
  Self_Prefix = character(),
  Self_Suffix = character(),
  stringsAsFactors = FALSE
)

# Loop over peptide pairs
for (i in 1:nrow(pairs)) {
  auto <- pairs$Autoimmune_Peptide[i]
  self <- pairs$Self_Peptide[i]
  
  shared <- find_longest_common_substring(auto, self)
  
  if (nchar(shared) == 0) next
  
  auto_start <- regexpr(shared, auto)[1]
  self_start <- regexpr(shared, self)[1]
  
  auto_prefix <- if (auto_start > 1) substr(auto, 1, auto_start - 1) else ""
  auto_suffix <- substr(auto, auto_start + nchar(shared), nchar(auto))
  
  self_prefix <- if (self_start > 1) substr(self, 1, self_start - 1) else ""
  self_suffix <- substr(self, self_start + nchar(shared), nchar(self))
  
  output <- rbind(output, data.frame(
    Pair_ID = i,
    Autoimmune_Peptide = auto,
    Self_Peptide = self,
    Shared_Overlap = shared,
    Auto_Prefix = auto_prefix,
    Auto_Suffix = auto_suffix,
    Self_Prefix = self_prefix,
    Self_Suffix = self_suffix,
    stringsAsFactors = FALSE
  ))
}

# Save to CSV
write.csv(output, "Peptide_Overlap_Decomposition.csv", row.names = FALSE)

###################### Align autoimmune and self-peptides ####################
# Load the decomposed result from the previous script
decomp <- read.csv("Peptide_Overlap_Decomposition.csv", stringsAsFactors = FALSE)

# Open PDF
pdf("All_Aligned_Peptides_FIXED.pdf", width = 10, height = 2)

for (i in 1:nrow(decomp)) {
  auto_prefix <- unlist(strsplit(decomp$Auto_Prefix[i], ""))   # red
  shared      <- unlist(strsplit(decomp$Shared_Overlap[i], ""))# yellow
  auto_suffix <- unlist(strsplit(decomp$Auto_Suffix[i], ""))   # red
  
  self_prefix <- unlist(strsplit(decomp$Self_Prefix[i], ""))   # green
  self_suffix <- unlist(strsplit(decomp$Self_Suffix[i], ""))   # green
  
  # Combine visual sequence: self_prefix + auto_prefix + shared + auto_suffix + self_suffix
  text <- c(self_prefix, auto_prefix, shared, auto_suffix, self_suffix)
  color <- c(
    rep("darkgreen", length(self_prefix)),
    rep("red",       length(auto_prefix)),
    rep("gold",      length(shared)),
    rep("red",       length(auto_suffix)),
    rep("darkgreen", length(self_suffix))
  )
  
  # Create data frame
  df <- data.frame(
    Position = 1:length(text),
    AminoAcid = text,
    Color = color
  )
  
  # Mini legend
  legend_text <- data.frame(
    Label = c("Autoimmune-only", "Self-only", "Shared"),
    Color = c("red", "darkgreen", "gold"),
    X = c(2, 15, 30),
    Y = rep(-1.5, 3)
  )
  
  # Plot
  p <- ggplot(df, aes(x = Position, y = 0, label = AminoAcid, color = Color)) +
    geom_text(size = 6, fontface = "bold") +
    geom_text(data = legend_text, aes(x = X, y = Y, label = Label, color = Color),
              inherit.aes = FALSE, size = 4) +
    scale_color_identity() +
    theme_void() +
    labs(title = paste("Peptide Comparison – Pair", decomp$Pair_ID[i])) +
    theme(plot.title = element_text(hjust = 0.5))
  
  print(p)
}

dev.off()
