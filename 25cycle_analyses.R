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
    grepl("PMBlock_A_C3_H_", Sample) ~ "Blocker1",
    grepl("PMBlock_B_C3_H_", Sample) ~ "Blocker2",
    grepl("PMBlock_A_dT_H_", Sample) ~ "Blocker3",
    grepl("PMBlock_B_dT_H_", Sample) ~ "Blocker4",
    grepl("PMBlock_A_C3", Sample) ~ "Blocker5",
    grepl("PMBlock_B_C3", Sample) ~ "Blocker6",
    grepl("PMBlock_A_dT", Sample) ~ "Blocker7",
    grepl("PMBlock_B_dT", Sample) ~ "Blocker8",
    TRUE ~ "No_Blocker" # default case if none of the above conditions are met
  ))

# elongate for bar charts
matrix_long <- matrix1 %>%
  pivot_longer(cols = c(Salvelinus_namaycush, Catostomus_commersonii, Petromyzontidae_unclassified, Salmonidae_unclassified), 
               names_to = "Species", values_to = "Reads")
# simplify sample and blocker names
matrix_long$Sample <- gsub("\\.12S$", "", matrix_long$Sample)
# matrix_long <- mutate(matrix_long, specific_blocker = str_replace_all(specific_blocker, c("^40_" = "", "^25_" = "")))
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


# ----------------------------------------------------------
# making a plot to directly compare 25 vs 40 cycle samples
# ----------------------------------------------------------

# filter to just compare 25-cycle samples with respective 40-cycle samples
matrix_long_25comparison <- matrix_long %>%
  filter(specific_blocker %in% c("No_Blocker", "Blocker2", "Blocker6"))

# create new data frame with combination of sample and blocker
matrix_long_25comparison <- matrix_long_25comparison %>%
  mutate(sample_blocker = paste(specific_blocker, Digestive_Sample, sep = "_"))

# NOTE: some of the 25 cycle columns (e.g. Blocker2_CA14 for 25 cycles) because these
# samples were removed during bioinformatic filtering, likely for low read coutnt.

# filter out NTC and M5 and only sea lamprey/lake trout
filtered_matrix_1 <- matrix_long_25comparison %>%
  filter(!str_detect(sample_blocker, "M5")) %>%
  filter(!str_detect(sample_blocker, "NTC")) %>%
  filter(Species %in% c("Petromyzontidae_unclassified", "Salvelinus_namaycush"))

# creating a factor interaction to adjust the colors 
filtered_matrix_2 <- filtered_matrix_1 %>%
  mutate(sample_group = interaction(specific_blocker, cycles, sep = "_"))

# custom color palette where:
#--- light colors = 25 cycles
#--- dark colors = 40 cycles

custom_colors2 <- c("Blocker2_25" = "#add8e6", 
                   "Blocker2_40" = "#324ca8", 
                   "Blocker6_25" = "#ffcc99", 
                   "Blocker6_40" = "#ff8c00", 
                   "No_Blocker_25" = "#90ee90", 
                   "No_Blocker_40" = "#0a8029") 

# reverse the factor levels
filtered_matrix_2 <- filtered_matrix_2 %>%
  mutate(sample_blocker = fct_relevel(sample_blocker, rev(levels(sample_blocker))))

# plot
ggplot(filtered_matrix_2, aes(x = sample_blocker, y = Reads, fill = sample_group)) +
  geom_bar(stat="identity", position=position_dodge(width = 0.75, preserve = "single")) + # preserve as some data is missing and this preserves column width
  facet_wrap(~Species) +
  scale_fill_manual(values = custom_colors2) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 70, vjust = 0.95, hjust=1, size=9),
        panel.spacing = unit(1, "lines")) +
  labs(x = "Sample",
       y = "Number of Sequence Reads",
       title = "Sequence Reads per Cycle Number",
       fill = "Blocker + Cycle Number") 


# This only looks at sea lamprey and lake trout -- similar results were shown for the single
# white sucker sample (M5)

# The no blocking primer samples (green) act as a baseline. The goal is to then look
# for differences between 40 and 25 cycle samples in their amplification of both sea
# lamprey and lake trout. Sea lamprey differences were minimal, while lake trout
# amplification was consistently higher in 40 cycle samples. Even detecting lake trout
# in sammples were the 25 cycle sample did not.
















