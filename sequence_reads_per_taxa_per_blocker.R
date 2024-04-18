# Sequence Read Counts per Taxa per Blocking Primer
# Conor O'Kane

library(tidyverse)
library(readxl)
library(ggthemes)


options(scipen = 999)  # removing scientific notation

matrix <- read_excel("./12S_SL_BP_Micro.mothur.community.matrix.xls")



#-------------------------------------------------------------------------------
# Data Organization

# change name of first column to indicate the sample
colnames(matrix)[1] <- "Sample"
matrix$Sample <- gsub("\\.12S$", "", matrix$Sample) #simplify sample names

# lot of extra columns, only need first four species
# these are lake trout, white sucker, sea lamprey, and Salmonidae 
matrix1 <- matrix[,1:5]

# filter out 25 cycles
matrix1 <- matrix1 %>%
  filter(grepl("40_", Sample))

# which sample?
matrix1 <- matrix1 %>%
  mutate(sample = case_when(
    grepl("CA14", Sample) ~ "CA14",
    grepl("HP3", Sample) ~ "HP3",
    grepl("HP15", Sample) ~ "HP15",
    grepl("HP5", Sample) ~ "HP5",
    grepl("M1", Sample) ~ "M1",
    grepl("M4", Sample) ~ "M4",
    grepl("M5", Sample) ~ "M5",
    TRUE ~ "NTC" # default case
  ))

# specific Sample column
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

# reorder to sort by blocker, and to start with no blocker
matrix1 <- matrix1[order(matrix1$specific_blocker), ]
matrix1 <- matrix1 %>%
  mutate(specific_blocker = factor(specific_blocker, levels = c(
    "No_Blocker", "Blocker1", "Blocker2", "Blocker3", "Blocker4",
    "Blocker5", "Blocker6", "Blocker7", "Blocker8"
  )))


# pivot longer to prepare for stacked bar chart
matrix1_long <- matrix1 %>%
  pivot_longer(cols = c(Salvelinus_namaycush, 
                        Catostomus_commersonii, 
                        Petromyzontidae_unclassified, 
                        Salmonidae_unclassified), 
               names_to = "Species", values_to = "Reads")

# add column to indicate sample names used 
matrix1_long <- matrix1_long %>%
  mutate(Digestive_Sample = case_when(
    grepl("CA14", Sample) ~ "CA14",
    grepl("HP3", Sample) ~ "HP3",
    grepl("HP15", Sample) ~ "HP15",
    grepl("HP5", Sample) ~ "HP5",
    grepl("M1", Sample) ~ "M1",
    grepl("M4", Sample) ~ "M4",
    grepl("M5", Sample) ~ "M5",
    TRUE ~ "NTC" # default case (non-template control)
  ))

matrix1_long$Species <- factor(matrix1_long$Species, levels=c("Salmonidae_unclassified",
                                                              "Salvelinus_namaycush",
                                                              "Catostomus_commersonii",
                                                              "Petromyzontidae_unclassified"))
# filter out homo sapiens and NTC samples
matrix1_long <- filter(matrix1_long, Species != "Homo_sapiens")
# matrix1_long <- filter(matrix1_long, !grepl("NTC", matrix1_long$Sample))




#-------------------------------------------------------------------------------
# Sequence reads per taxa per blocker plot

# adjust colors
custom_colors <- c("#6bc9b2", "#40aac2","#87a7e3","#01161e")

#change order of factor
matrix1_long$Digestive_Sample <- factor(matrix1_long$Digestive_Sample, levels = c("M1", "M4", "M5", "HP3", "HP5", "HP15", "CA14", "NTC"))

#plot
matrix1_long_plot <- ggplot(matrix1_long, aes(x = specific_blocker, y = Reads, fill = Species)) +
  geom_col() +
  scale_fill_manual(values = custom_colors, 
                    labels = c("Salmonid (unclassified)",  "Lake Trout", "White Sucker", "Sea Lamprey"))+
  theme_few()+
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1, size = 15),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("Sequence Read Counts per Taxa per Blocking Primer") +
  labs(x = "Blocking Primer", y = "Number of Sequence Reads") +
  facet_wrap(~ Digestive_Sample, ncol = 2)+  # toggle this on/off for wrapped/unwrapped plot
  theme(text = element_text(family = "Helvetica", face = "bold"),
        strip.text = element_text(size = 14),
        axis.text = element_text(family = "Helvetica", face = "plain"),
        axis.text.x = element_text(size = 9),
        legend.key.size = unit(1.5, "lines"),  
        legend.text = element_text(size = 12),  
        legend.title = element_text(size = 14, face = "bold"))

matrix1_long_plot

#-------------------------------------------------------------------------------
# Simplified plot for poster presentation

# only looking to plot M4, M5, HP15, and CA14
matrix1_simple <- matrix1_long %>%
  filter(Digestive_Sample %in% c("M4", "M5", "HP15", "CA14"))

#drop unused levels
matrix1_simple$Digestive_Sample <- droplevels(matrix1_simple$Digestive_Sample)

# change facet labels for clarity
levels(matrix1_simple$Digestive_Sample) <- c("M4 (90:10 SL:LT)",
                                             "M5 (90:10 SL:WS)",
                                             "HP15 (Wild Juvenile)",
                                             "CA14 (Wild Adult)")

# plot as before
matrix1_simple_plot <- ggplot(matrix1_simple, aes(x = specific_blocker, y = Reads, fill = Species)) +
  geom_col() +
  scale_fill_manual(values = custom_colors, 
                    labels = c("Salmonid",  "Lake Trout", "White Sucker", "Sea Lamprey"))+
  theme_few()+
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1, size = 15),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("Sequence Read Counts per Taxa per Blocking Primer") +
  labs(x = "Blocking Primer", y = "Number of Sequence Reads") +
  facet_wrap(~ Digestive_Sample, ncol = 2)+  # toggle this on/off for wrapped/unwrapped plot
  theme(text = element_text(family = "Helvetica", face = "bold"),
        strip.text = element_text(size = 14),
        axis.text = element_text(family = "Helvetica", face = "plain"),
        axis.text.x = element_text(size = 9),
        legend.key.size = unit(1.25, "lines"),  
        legend.text = element_text(size = 11),  
        #legend.title = element_text(size = 14, face = "bold"),
        legend.title = element_blank(),
        legend.position = "bottom"  
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

matrix1_simple_plot



#-------------------------------------------------------------------------------
# Statistical comparisons

head(matrix1_long)

# average sea lamprey reads in mock communities with blocking primer
lamprey_blocker_reads <- matrix1_long %>%
  filter(Species == "Petromyzontidae_unclassified") %>%
  filter(specific_blocker != "No_Blocker") %>%
  filter(sample %in% c("M1", "M4", "M5"))

round(mean(lamprey_blocker_reads$Reads), digits = 2) # 1.03
round(sd(lamprey_blocker_reads$Reads), digits = 2) # 1.38



