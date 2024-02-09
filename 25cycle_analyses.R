# Comparison of 25-Cycle Sequences

# libraries
library(tidyverse)
library(readxl)

# data set up

matrix <- read_excel("./12S_SL_BP_Micro.mothur.community.matrix.xls")
# change name of first column to indicate the sample
colnames(matrix)[1] <- "Sample"
matrix1 <- matrix[,1:5]

matrix1$block <- ifelse(grepl("NoBlock", matrix1$Sample), "No", "Yes")
matrix1$cycles <- ifelse(grepl("25", matrix1$Sample), "25", "40")
matrix1 <- matrix1 %>%
  mutate(specific_blocker = case_when(
    grepl("40_PMBlock_A_C3_H_", Sample) ~ "40_A_C3_H",
    grepl("40_PMBlock_B_C3_H_", Sample) ~ "40_B_C3_H",
    grepl("40_PMBlock_A_dT_H_", Sample) ~ "40_A_dT_H",
    grepl("40_PMBlock_B_dT_H_", Sample) ~ "40_B_dT_H",
    grepl("40_PMBlock_A_C3", Sample) ~ "40_A_C3",
    grepl("40_PMBlock_B_C3", Sample) ~ "40_B_C3",
    grepl("40_PMBlock_A_dT", Sample) ~ "40_A_dT",
    grepl("40_PMBlock_B_dT", Sample) ~ "40_B_dT",
    grepl("25_PMBlock_B_C3_H_", Sample) ~ "25_B_C3_H",
    grepl("25_PMBlock_B_C3", Sample) ~ "25_B_C3",
    TRUE ~ "No_Blocker" # default case if none of the above conditions are met
  ))

# elongate for bar charts
matrix_long <- matrix1 %>%
  pivot_longer(cols = c(Salvelinus_namaycush, Catostomus_commersonii, Petromyzontidae_unclassified, Salmonidae_unclassified), 
               names_to = "Species", values_to = "Reads")
# simplify sample and blocker names
matrix_long$Sample <- gsub("\\.12S$", "", matrix_long$Sample)
matrix_long <-   mutate(matrix_long, specific_blocker = str_replace_all(specific_blocker, c("^40_" = "", "^25_" = "")))
# add digestive sample column
matrix_long <- matrix_long %>%
  mutate(Digestive_Sample = case_when(
    grepl("CA14", Sample) ~ "CA14",
    grepl("HP3", Sample) ~ "HP3",
    grepl("HP15", Sample) ~ "HP15",
    grepl("HP5", Sample) ~ "HP5",
    grepl("M1", Sample) ~ "M1",
    grepl("M4", Sample) ~ "M4",
    grepl("M5", Sample) ~ "M5",
    TRUE ~ "NTC" # default case
  ))

# some colors
custom_colors <- c("#50b299","#87a7e3", "#7d5388","#01161e")
custom_colors2 <- c("#7ac4ae", "#1f4a3d")
#-------- some basic stats n plots
#colors
custom_colors <- c("#50b299","#87a7e3", "#7d5388","#01161e")
# all samples
ggplot(matrix_long, aes(x = Sample, y = Reads, fill = Species)) +
  geom_col() +
  scale_fill_manual(values = custom_colors) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size =7))

# filter to just compare 25-cycle samples with respective 40-cycle samples
matrix_long_25comparison <- matrix_long %>%
  filter(specific_blocker %in% c("No_Blocker", "B_C3", "B_C3_H"))
# plot as side by side barplot
ggplot(matrix_long_25comparison, aes(x = specific_blocker, y = Reads, fill = cycles)) +
  geom_bar(stat="identity", position=position_dodge(width = 0.75)) + # Make sure to set the dodge width
  facet_wrap(~Species) + # create a separate panel for each species
  scale_fill_manual(values = custom_colors2) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7)) +
  labs(x = "Sample", 
       y = "Number of Sequence Reads", 
       title = "Sequence Reads per Cycle Number")















