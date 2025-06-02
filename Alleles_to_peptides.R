allele_output <- read.table("output.csv", stringsAsFactors = FALSE, sep = "\t", header = TRUE)
allele_output <- unique(allele_output)
supplementary_data <- read.csv("supplementary_table1 (1).csv")

####################### I will have to do it in 4 loops for processing purposes ###################

# Rename columns correctly (Run this only once)
colnames(supplementary_data)[1:3] <- c("Donor", "Peptide", "Gene_ID")
colnames(allele_output)[2:4] <- c("Peptide", "Binding_Core", "HLA_Allele")

# Get the list of unique donors
unique_donors <- unique(supplementary_data$Donor)


################# 1-10 ################
batch_1_10 <- data.frame(HLA_Allele = character(), Peptide = character(), Gene_ID = character(), Donor = character(), stringsAsFactors = FALSE)

filtered_supplementary <- supplementary_data[supplementary_data$Donor %in% unique_donors[1:10], ]

for (i in 1:nrow(allele_output)) {
  peptide_match <- allele_output$Peptide[i]
  matching_rows <- filtered_supplementary[filtered_supplementary$Peptide == peptide_match, ]
  
  if (nrow(matching_rows) > 0) {
    for (j in 1:nrow(matching_rows)) {
      batch_1_10 <- rbind(batch_1_10, data.frame(
        HLA_Allele = allele_output$HLA_Allele[i],
        Peptide = matching_rows$Peptide[j],
        Gene_ID = matching_rows$Gene_ID[j],
        Donor = matching_rows$Donor[j]
      ))
    }
  }
}



################# 11-20 ################
batch_11_20 <- data.frame(HLA_Allele = character(), Peptide = character(), Gene_ID = character(), Donor = character(), stringsAsFactors = FALSE)

filtered_supplementary <- supplementary_data[supplementary_data$Donor %in% unique_donors[11:20], ]

for (i in 1:nrow(allele_output)) {
  peptide_match <- allele_output$Peptide[i]
  matching_rows <- filtered_supplementary[filtered_supplementary$Peptide == peptide_match, ]
  
  if (nrow(matching_rows) > 0) {
    for (j in 1:nrow(matching_rows)) {
      batch_11_20 <- rbind(batch_11_20, data.frame(
        HLA_Allele = allele_output$HLA_Allele[i],
        Peptide = matching_rows$Peptide[j],
        Gene_ID = matching_rows$Gene_ID[j],
        Donor = matching_rows$Donor[j]
      ))
    }
  }
}



################# 21-30 ################
batch_21_30 <- data.frame(HLA_Allele = character(), Peptide = character(), Gene_ID = character(), Donor = character(), stringsAsFactors = FALSE)

filtered_supplementary <- supplementary_data[supplementary_data$Donor %in% unique_donors[21:30], ]

for (i in 1:nrow(allele_output)) {
  peptide_match <- allele_output$Peptide[i]
  matching_rows <- filtered_supplementary[filtered_supplementary$Peptide == peptide_match, ]
  
  if (nrow(matching_rows) > 0) {
    for (j in 1:nrow(matching_rows)) {
      batch_21_30 <- rbind(batch_21_30, data.frame(
        HLA_Allele = allele_output$HLA_Allele[i],
        Peptide = matching_rows$Peptide[j],
        Gene_ID = matching_rows$Gene_ID[j],
        Donor = matching_rows$Donor[j]
      ))
    }
  }
}

################# 31-40 ################
batch_31_40 <- data.frame(HLA_Allele = character(), Peptide = character(), Gene_ID = character(), Donor = character(), stringsAsFactors = FALSE)

filtered_supplementary <- supplementary_data[supplementary_data$Donor %in% unique_donors[31:40], ]

for (i in 1:nrow(allele_output)) {
  peptide_match <- allele_output$Peptide[i]
  matching_rows <- filtered_supplementary[filtered_supplementary$Peptide == peptide_match, ]
  
  if (nrow(matching_rows) > 0) {
    for (j in 1:nrow(matching_rows)) {
      batch_31_40 <- rbind(batch_31_40, data.frame(
        HLA_Allele = allele_output$HLA_Allele[i],
        Peptide = matching_rows$Peptide[j],
        Gene_ID = matching_rows$Gene_ID[j],
        Donor = matching_rows$Donor[j]
      ))
    }
  }
}


################ lægger det nu ind i kæmpe fil ###############
# Combine all four batches into one dataset
combined_data <- rbind(batch_1_10, batch_11_20, batch_21_30, batch_31_40)
combined_data <- unique(combined_data)

#Sortér efter alleller og så gener:
sorted_data <- combined_data[order(combined_data$HLA_Allele, combined_data$Gene_ID), ]

#Fjern Donorer
sorted_data <- sorted_data[, !(colnames(sorted_data) %in% "Donor")]
sorted_data <- unique(sorted_data)

# Save as a CSV file




################ Jeg tilføjer lige core_peptides ######################### 
# Ensure correct column names based on structure
colnames(allele_output)[2:3] <- c("Peptide", "Core_Peptide")

# Add a new Core_Peptide column while preserving existing data
sorted_data$Core_Peptide <- NA  # Initialize new column

# Scan each peptide in sorted_data and extract the core peptide from allele_output
for (i in 1:nrow(sorted_data)) {
  peptide_match <- sorted_data$Peptide[i]  # Get peptide from sorted_data
  
  # Find corresponding core peptide in allele_output
  match_row <- allele_output[allele_output$Peptide == peptide_match, "Core_Peptide"]
  
  # If a match is found, update Core_Peptide column
  if (length(match_row) > 0) {
    sorted_data$Core_Peptide[i] <- match_row[1]  # Take first match if multiple exist
  }
}

# Og skifter også lige søjler så peptider står ved siden af cores
sorted_data <- sorted_data[, c(1, 2, 4, 3)]

#Og nu gemmer jeg så
write.csv(sorted_data, "Peptides_classified.csv", row.names = FALSE)


