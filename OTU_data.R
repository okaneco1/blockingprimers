library(tidyverse)
library(readxl)
library(ggthemes)


options(scipen = 999)  # removing scientific notation

matrix <- read_excel("./12S_SL_BP_Micro.mothur.community.matrix.xls")
head(matrix)

# change name of first column to indicate the sample
colnames(matrix)[1] <- "Sample"

sums <- colSums(matrix[2:24])
sums <- as.data.frame(sums)

sums$species <- rownames(sums)
rownames(sums) <- NULL  # remove row names from the data frame
sums[,c(1,2)] <- sums[,c(2,1)] # reorder
colnames(sums) <- c("species", "readcounts") 

sums$species <- factor(sums$species, levels = sums$species[order(-sums$readcounts)]) #reorder based on read count

total_reads_barchart <- ggplot(data = sums, aes(x = species, y = readcounts)) +
  geom_bar(stat = "identity", fill = "#292a40") +
  geom_text(aes(label = readcounts), vjust = -0.5, size = 3, fontface = "bold") +
  theme_few()+
  labs(title = "Total Sequence Read Counts per OTU", x = "Operational Taxonomic Unit (OTU)", y = "Total Read Count") +
  theme(axis.text.x = element_text(angle = 70, vjust = 1, hjust=1),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold"))

total_reads_barchart

# pie chart for the total amount of reads
total_reads_piechart <- ggplot(data = sums, aes(x = "", y = readcounts, fill = species)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=58) +
  theme_void() +
  scale_fill_brewer(palette="Accent")+
  labs(fill = "Species")

total_reads_piechart 

# Final figure was created using these graphs and further visual manipulation 
# in Adobe Photoshop


#---------- some basic stats

nrow(sums) # 23 OTU

total <- sum(sums$readcounts) # 2090328 total reads
total4 <- sum(sums$readcounts[1:4]) # 2056721 total reads in first four
round((total4/total)*100, digits = 2) # 96.76% in first four














