# Load required libraries
library(dplyr)
library(ggplot2)

# Step 1: Load data
my_data <- read.csv("autoimmune_iedb_pepetides.csv")

# Step 2: Filter to include only HLA class II interactions
filtered_data <- my_data

# Step 3: Sanitize column names for easier reference
colnames(filtered_data) <- make.names(colnames(filtered_data))

# Step 4: Rename key columns for easier use
filtered_data <- filtered_data %>%
  rename(
    Disease = DiseaseName,
    Protein = Epitope..Parent.Protein.Accession
  )

# Step 5: Group by Disease and Protein, and count occurrences
protein_disease_df <- filtered_data %>%
  group_by(Disease, Protein) %>%
  summarise(Occurrences = n(), .groups = "drop")

# Step 6: Create the ggplot object
p <- ggplot(protein_disease_df, aes(x = reorder(Protein, Occurrences), y = Occurrences, fill = Occurrences)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  facet_wrap(~ Disease, scales = "free_y", ncol = 3) +  
  labs(title = "Protein Occurrences by Disease", x = "Proteins", y = "Occurrences") +
  scale_fill_gradient(low = "red", high = "green") +  
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12, face = "bold", margin = margin(b = 8, t = 8)),  
    strip.background = element_rect(fill = "white", color = "black"),  
    axis.text.y = element_text(size = 10, margin = margin(r = 8)),  
    axis.text.x = element_text(angle = 45, hjust = 1)  
  )

# Step 7: Save the plot as PNG
ggsave(filename = "Final_diseases.png", plot = p, width = 22, height = 36, dpi = 300)
