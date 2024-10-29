# Visualizing Blocking Primer Overlap Using ggmsa

# libraries
library(ggmsa)
library(Biostrings)
library(tidyverse)

# import MSA data
nt_sequences <- readDNAStringSet("./FishOnly_12S_align_JK_072821.fas")

# trim to first 20 sequences
trimmed_seqs <- nt_sequences[1:30]

# and save as new fasta file for ggmsa
writeXStringSet(trimmed_seqs, "./trimmed_sequences_20.fas")

# custom color palette
my_custom <- data.frame(names = c("A","G","T","C"), 
                        colors = c("#4371f0","#43f082","#f0b643","#f04943"))

# visualize
ggmsa("./trimmed_sequences_20.fas", 1, 75, 
      custom_color = my_custom, 
      font = "DroidSansMono", 
      char_width = 0.5, 
      seq_name = TRUE )




# seeing how long without primers
no_primer <- readDNAStringSet("./FishOnly_12S_align_JK_noprimers_091922.fas")
top_four <- no_primer[1:4]
pairwiseAlignment()

calc_differences <- function(seq1, seq2) {
  sum(as.character(seq1) != as.character(seq2))
}

# Create a matrix to store pairwise differences
combs <- combn(1:4, 2)
diff_matrix <- matrix(NA, ncol = 3, nrow = ncol(combs))
colnames(diff_matrix) <- c("Seq1", "Seq2", "Differences")

# Loop through the pairwise combinations
for (i in 1:ncol(combs)) {
  seq1 <- top_four[[combs[1, i]]]
  seq2 <- top_four[[combs[2, i]]]
  diff_matrix[i, ] <- c(combs[1, i], combs[2, i], calc_differences(seq1, seq2))
}

diff_matrix




# Comparing differences in bp between Coregoninae
library(ape)
# trim to Coregonus sequences
whitefish_sequences <- nt_sequences[grep("Coregonus|Prosopium|Stenodus", names(nt_sequences))]
# convert to DNAbin object
whitefish_dnabin <- as.DNAbin(whitefish_sequences)
# calculate distance matrix n
distance_matrix <- dist.dna(whitefish_dnabin, model = "N")  
range(distance_matrix)







