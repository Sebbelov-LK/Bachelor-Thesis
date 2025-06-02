# Bachelor-Thesis
This repository, contains the coding required to replicate the result section of the bachelor thesis named Comparative Analysis of HLA Binding Patterns: Tolerogenic Self-Peptides vs. Autoimmune Disease-Associated Epitopes, created in collaboration with Herbert Grunkin. 

The first file contains the full protein analysis R script, spanning half of the result section up until the subsubsection Distribution of proteins over diseases. The second half of the result section is divided up into multiple files following the narrative: 

1. Disease by protein
  a. Overall Diseases 2.R
2. Cutting position heatmaps
  a. Cutting Points and Heatmaps.R
3. Input til deconvolution software
  a. Data_connections.R
4. Peptide lengths and histograms
  a. Alleles_to_peptides.R (Processing deconvolutions)
  b. Extended_peptides.R
  c. Peptide_lengths.R
5. Treemaps over cleavage pairs
  a. Treeplots.R
6. Prefix & Suffix analysis
  a. Autoimmune vs Self-peptide.R comparisons (appendix file)
  b. Matched_peptides_analysis.R (treemap & center distance)
7. Differential Logo Making
   a. ICElogo making.R (files to insert into seq2logo)
  b. BICElogo.R

