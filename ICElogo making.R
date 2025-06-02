############################################################################
############### FUNCTION TO FORMAT & SAVE FREQUENCIES IN SEQ2LOGO FORMAT ###############
write_seq2logo_freq <- function(peptides, output_filename) {
  all_aa <- paste(peptides, collapse = "")
  aa_vector <- unlist(strsplit(all_aa, ""))
  
  # Ensure we include all 20 canonical amino acids in the output
  canonical_aa <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
                    "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
  
  freq_table <- table(factor(aa_vector, levels = canonical_aa))
  freq_normalized <- freq_table / sum(freq_table)
  
  # Write to file
  writeLines(paste(canonical_aa, collapse = "\t"), con = output_filename)
  write.table(t(freq_normalized), file = output_filename, append = TRUE,
              col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
}

####################### BACKGROUND FREQUENCY FILES #######################

# Self peptides - First10
peptides_self_first <- readLines("First10_Self_Peptides.txt")
write_seq2logo_freq(peptides_self_first, "Self_Background_Frequencies_First10.txt")

# Self peptides - Last10
peptides_self_last <- readLines("Last10_Self_Peptides.txt")
write_seq2logo_freq(peptides_self_last, "Self_Background_Frequencies_Last10.txt")

# Autoimmune peptides - First10
peptides_auto_first <- readLines("First10_Extended_Autoimmune.txt")
write_seq2logo_freq(peptides_auto_first, "Autoimmune_Background_Frequencies_First10.txt")

# Autoimmune peptides - Last10
peptides_auto_last <- readLines("Last10_Extended_Autoimmune.txt")
write_seq2logo_freq(peptides_auto_last, "Autoimmune_Background_Frequencies_Last10.txt")


################################################################################
######################## PSSM.txt into dataframe ###############################
# Function to load and clean a PSSM matrix from file
load_pssm <- function(file_path) {
  lines <- readLines(file_path)
  aa_labels <- unlist(strsplit(lines[1], "\\s+"))
  data_lines <- lines[-1]
  split_data <- strsplit(data_lines, "\\s+")
  pssm_values <- lapply(split_data, function(x) as.numeric(x[3:length(x)]))
  pssm_matrix <- do.call(rbind, pssm_values)
  colnames(pssm_matrix) <- aa_labels
  rownames(pssm_matrix) <- sapply(split_data, function(x) x[1])
  return(pssm_matrix)
}

# Function to save a PSSM matrix in Seq2Logo format
save_pssm <- function(matrix, output_file) {
  cat(colnames(matrix), file = output_file, sep = "\t")
  cat("\n", file = output_file, append = TRUE)
  for (i in 1:nrow(matrix)) {
    line <- paste(rownames(matrix)[i], "L", paste(sprintf("%.3f", matrix[i, ]), collapse = "\t"))
    cat(line, file = output_file, append = TRUE, sep = "\n")
  }
}

# First 10 comparison
pssm_ai_first <- load_pssm("First 10 AI.txt")
pssm_self_first <- load_pssm("First 10 Self.txt")
stopifnot(identical(dim(pssm_ai_first), dim(pssm_self_first)))
stopifnot(all(colnames(pssm_ai_first) == colnames(pssm_self_first)))
pssm_diff_first <- pssm_ai_first - pssm_self_first
save_pssm(pssm_diff_first, "PSSM_Difference_First10.txt")

# Last 10 comparison
pssm_ai_last <- load_pssm("Last 10 AI.txt")
pssm_self_last <- load_pssm("Last 10 Self.txt")
stopifnot(identical(dim(pssm_ai_last), dim(pssm_self_last)))
stopifnot(all(colnames(pssm_ai_last) == colnames(pssm_self_last)))
pssm_diff_last <- pssm_ai_last - pssm_self_last
save_pssm(pssm_diff_last, "PSSM_Difference_Last10.txt")

###################### Export Difference matrices ###########################

###################### Autoimmune minus Self
# Export First10 difference matrix
pos_labels_first <- paste(1:nrow(pssm_diff_first), rownames(pssm_diff_first))
export_first <- cbind(Position = pos_labels_first, round(pssm_diff_first, 3))
write.table(export_first, file = "PSSM_Difference_First10.txt",
            sep = "\t", quote = FALSE, row.names = FALSE,
            col.names = c("", colnames(pssm_diff_first)))

# Export Last10 difference matrix
pos_labels_last <- paste(1:nrow(pssm_diff_last), rownames(pssm_diff_last))
export_last <- cbind(Position = pos_labels_last, round(pssm_diff_last, 3))
write.table(export_last, file = "PSSM_Difference_Last10.txt",
            sep = "\t", quote = FALSE, row.names = FALSE,
            col.names = c("", colnames(pssm_diff_last)))




# Define export function in Seq2Logo format
save_pssm_seq2logo <- function(pssm_matrix, file_name) {
  # Get amino acid headers
  aa_headers <- colnames(pssm_matrix)
  
  # Write header (just the amino acid letters, tab-separated)
  header_line <- paste(aa_headers, collapse = "\t")
  writeLines(header_line, file_name)
  
  # Create reference AA column (default to "L" or update with real data if available)
  reference_aa <- rep("L", nrow(pssm_matrix))
  
  # Combine into export format
  export_matrix <- cbind(Position = seq_len(nrow(pssm_matrix)), AA = reference_aa, round(pssm_matrix, 3))
  
  # Write matrix rows
  write.table(export_matrix, file = file_name, append = TRUE, sep = "\t",
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}

# Now apply the export to your matrices

# First10 difference matrix
save_pssm_seq2logo(pssm_diff_first, "PSSM_Difference_First10.txt")

# Last10 difference matrix
save_pssm_seq2logo(pssm_diff_last, "PSSM_Difference_Last10.txt")

###################### Self minus autoimmune
# Reverse subtraction: Self - Autoimmune

# First 10 positions
pssm_diff_first_rev <- pssm_self_first - pssm_ai_first
pos_labels_first_rev <- paste(1:nrow(pssm_diff_first_rev), rownames(pssm_diff_first_rev))
export_first_rev <- cbind(Position = pos_labels_first_rev, round(pssm_diff_first_rev, 3))
write.table(export_first_rev, file = "PSSM_Difference_First10_Rev.txt",
            sep = "\t", quote = FALSE, row.names = FALSE,
            col.names = c("", colnames(pssm_diff_first_rev)))

# Last 10 positions
pssm_diff_last_rev <- pssm_self_last - pssm_ai_last
pos_labels_last_rev <- paste(1:nrow(pssm_diff_last_rev), rownames(pssm_diff_last_rev))
export_last_rev <- cbind(Position = pos_labels_last_rev, round(pssm_diff_last_rev, 3))
write.table(export_last_rev, file = "PSSM_Difference_Last10_Rev.txt",
            sep = "\t", quote = FALSE, row.names = FALSE,
            col.names = c("", colnames(pssm_diff_last_rev)))