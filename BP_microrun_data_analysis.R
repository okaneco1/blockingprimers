library(tidyverse)
library(readxl)
library(ggthemes)


options(scipen = 999)  # removing scientific notation

matrix <- read_excel("./12S_SL_BP_Micro.mothur.community.matrix.xls")
head(matrix)

# change name of first column to indicate the sample
colnames(matrix)[1] <- "Sample"

# lot of extra columns, only need first 5
# these are lake trout, white sucker, sea lamprey, Salmonidae, and human 
matrix1 <- matrix[,1:6]
head(matrix1)


###############################################################
# Adding Categories/Organizing the data
###############################################################

# is there a block or no block?
matrix1$block <- ifelse(grepl("NoBlock", matrix1$Sample), "No", "Yes")

# cycle number
matrix1$cycles <- ifelse(grepl("25", matrix1$Sample), "25", "40")

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

# block A or block B?
matrix1$A_B <- ifelse(grepl("_A_", matrix1$Sample), "A", 
                      ifelse(grepl("_B_", matrix1$Sample), "B","No_Blocker"))

# which end modification?
matrix1$end_mod <- ifelse(grepl("_C3_", matrix1$Sample), "C3", 
                          ifelse(grepl("_dT_", matrix1$Sample), "dT","No_Blocker"))

# HPLC purified?
matrix1$purification <- ifelse(grepl("_H_", matrix1$Sample), "H", 
                               ifelse(grepl("NoBlock", matrix1$Sample), "No_Blocker","non_H"))

# specific Sample column
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

#####################################################################
# Making Graphs
#####################################################################

library(ggrepel) # some of the text was overlapping, so this fixes that 
# comparing average number of sea lamprey sequences in bp samples vs non-bp samples

# for sea lamprey
ggplot(data=matrix1,aes(x=specific_blocker, y=Petromyzontidae_unclassified)) +
  geom_boxplot()

# for lake trout      
ggplot(data=matrix1,aes(x=specific_blocker, y=Salvelinus_namaycush)) +
  geom_boxplot()
# but these don't account for the samples that don't include lake trout, so...
# select only lake trout samples by removing M5 and NTC
matrix1 %>%
  filter(!grepl("M5", Sample) & !grepl("NTC", Sample)) %>%
  ggplot(aes(specific_blocker, Salvelinus_namaycush, fill = cycles))+
  geom_boxplot(alpha=0.8, outlier.shape = NA)+
  # bit messy with jitter or text, so can turn on/off
  geom_jitter(width = 0.2, height = 0, color = "black", alpha = 0.5) +
  geom_text_repel(data = subset(matrix1, sample %in% c("CA14","HP3")),
    aes(label = sample, alpha=0.8), vjust = 2, color = "black", size = 3, fontface = "bold") +
  scale_fill_manual(values=c("darkgray","#64CCC5")) +
  labs(x="Blocking Primer", y="Lake Trout Reads") +
  theme_classic()

# looking at this graph and the data a bit more closely, it looks like CA14 may have
# some outliers, so let's remove those as well and see what it looks like


matrix1 %>%
  filter(!grepl("M5", Sample) & !grepl("NTC", Sample) & !grepl("CA14", Sample)) %>%
  ggplot(aes(specific_blocker, Salvelinus_namaycush, fill = cycles))+
  geom_boxplot(alpha=0.8, outlier.shape = NA)+
  geom_jitter(width = 0.2, height = 0, color = "black", alpha = 0.5) +
  geom_text_repel(aes(label = sample, alpha=0.8), vjust = 2, color = "black", size = 3, fontface = "bold") +
  scale_fill_manual(values=c("darkgray","#64CCC5")) +
  labs(x="Blocking Primer", y="Lake Trout Reads") +
  theme_classic()

# -- Some Notes --
# Those points at the bottom for the 25 cycle samples are HP3 and HP15. When comparing
# to the 40 cycle samples, the one bottom outlier is HP3 and HP15 shoots up in with the 
# rest of the samples. Much more effective at the higher cycle ranges.


# let's try for white sucker
matrix1 %>%
  filter(grepl("M5", Sample) & !grepl("NTC", Sample) & !grepl("CA14", Sample)) %>%
  ggplot(aes(x = specific_blocker, y = Catostomus_commersonii, fill = cycles))+
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6) +
  scale_fill_manual(values=c("darkgray","#64CCC5")) +
  labs(x="Blocking Primer", y="White Sucker Reads") +
  theme_classic()

# and finally for sea lamprey
matrix1 %>%
  filter(!grepl("NTC", Sample)) %>%
  ggplot(aes(specific_blocker, Petromyzontidae_unclassified, fill = cycles))+
  geom_boxplot(alpha=0.8, outlier.shape = NA) +
  geom_jitter(width = 0.2, height = 0, color = "black", alpha = 0.5) +
  #geom_text_repel(aes(label = sample, alpha=0.8), vjust = 2, color = "black", size = 3, fontface = "bold") +
  scale_fill_manual(values=c("darkgray","#64CCC5")) +
  labs(x="Blocking Primer", y="Sea Lamprey Reads") +
  theme_classic()

## Additional Plots?:
# ---- Add in "salmonidae" with the lake trout?
# ---- Can compare specific variables to each other (e.g. A vs B or H vs non-H, for paper)
# ---- What are the proportion changes that the Samples make?
# --------- Compare the percentages by which the Samples reduced sea lamprey/lake trout
# ---- How the fraction of reads changes across categories and Samples




##### random stats ##########

matrix1 %>% 
  filter(specific_blocker == "No_Blocker" & cycles == "25") %>%
  summarize(mean = mean(Petromyzontidae_unclassified))  # 12211

matrix1 %>% 
  filter(specific_blocker == "No_Blocker" & cycles == "40") %>%
  summarize(mean = mean(Petromyzontidae_unclassified))  # 19387
  

# comparing the two cycle groups for number of lake trout reads

group25 <- subset(matrix1, cycles == "25" & specific_blocker != "No_Blocker" & sample != "NTC" & sample != "M5")
group40 <- subset(matrix1, cycles == "40" & specific_blocker != "No_Blocker" & sample != "NTC" & sample != "M5")
groupNB <- subset(matrix1, cycles == "40" & specific_blocker == "No_Blocker" & sample != "NTC" & sample != "M5")

summarize(groupNB, mean = mean(groupNB$Salvelinus_namaycush))
summarize(group40, mean = mean(group40$Salvelinus_namaycush))


result <- wilcox.test(group25$Salvelinus_namaycush, group40$Salvelinus_namaycush)

group25 %>%
  filter(specific_blocker == "25_B_C3") %>%
  summarize(mean = mean(Salvelinus_namaycush))

group40 %>%
  filter(specific_blocker == "40_B_C3") %>%
  summarize(mean = mean(Salvelinus_namaycush))

wilcox.test(group40$Salvelinus_namaycush, groupNB$Salvelinus_namaycush)

#just B-C3 and B-C3-H

group40_BC3 <- group40 %>%
  filter(A_B == "B" & end_mod == "C3")
wilcox.test(group25$Salvelinus_namaycush, group40_BC3$Salvelinus_namaycush)


# looking at salmonid reads per blocking primer

# make a new column with total reads 
matrix1$total_reads <- rowSums(matrix1[3:6])

# combine lake trout and Salmonidae for total salmonids (probably all lake trout)
matrix1$total_salmonids <- matrix1$Salmonidae_unclassified + matrix1$Salvelinus_namaycush

matrix1 %>%
  filter(!grepl("M5", Sample) & !grepl("NTC", Sample)) %>%
  ggplot(aes(specific_blocker, total_salmonids, fill = cycles))+
  geom_boxplot(alpha=0.8, outlier.shape = NA)+
  geom_jitter(width = 0.2, height = 0, color = "black", alpha = 0.5) +
  geom_text_repel(data = subset(matrix1, sample %in% c("CA14","HP3")),
                  aes(label = sample, alpha=0.8), vjust = 2, color = "black", size = 3, fontface = "bold") +
  scale_fill_manual(values=c("darkgray","#64CCC5")) +
  labs(x="Blocking Primer", y="Lake Trout Reads") +
  theme_classic()

# reorder and create stacked bar column

reordered <- matrix1[order(matrix1$specific_blocker), ]

reordered_long <- reordered[ ,c(1:6,13,14)] %>%
  pivot_longer(cols = c(Salvelinus_namaycush, Catostomus_commersonii, Petromyzontidae_unclassified, Salmonidae_unclassified, Homo_sapiens), 
               names_to = "Species", values_to = "Reads")

reordered_long <- reordered_long %>%
  group_by(specific_blocker) %>%
  mutate(Proportion = Reads / sum(Reads))

reordered_long <- reordered_long %>%
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

# taking ".12S" out of the sample names
reordered_long$Sample <- gsub("\\.12S$", "", reordered_long$Sample)

#colors
custom_colors <- c("#50b299", "#40aac2","#87a7e3", "#7d5388","#01161e")

# all samples
ggplot(reordered_long, aes(x = Sample, y = Reads, fill = Species)) +
  geom_col() +
  scale_fill_manual(values = custom_colors) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size =7))

# grouped by blocker
ggplot(reordered_long, aes(x = specific_blocker, y = Reads, fill = Species)) +
  geom_col() +
  scale_fill_manual(values = custom_colors) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size =7))

# CA14 with different blockers
ggplot(subset(reordered_long, grepl("CA14", Sample)), aes(x = specific_blocker, y = Reads, fill = Species)) +
  geom_col() +
  scale_fill_manual(values = custom_colors) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size =7))+
  labs(x = "Blocking Primer", y = "Number of Sequence Reads", title = "Number of Sequence Reads per Blocking Primer of CA14")

# HP15 with different blockers
ggplot(subset(reordered_long, grepl("HP15", Sample)), aes(x = specific_blocker, y = Reads, fill = Species)) +
  geom_col() +
  scale_fill_manual(values = custom_colors) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 11),
        plot.title = element_text(face = "bold", family = "Helvetica", size = 14)) +
  labs(x = "Blocking Primer", y = "Number of Sequence Reads")

unique(matrix1$sample)


library(ggthemes)


ggplot(reordered_long, aes(x = specific_blocker, y = Reads, fill = Species)) +
  geom_col() +
  scale_fill_manual(values = custom_colors, 
                    labels = c("White Sucker", "Homo Sapiens", "Sea Lamprey", "Salmonid (unclassified)", "Lake Trout"))+
  theme_few()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 11)) +
  labs(x = "Blocking Primer", y = "Number of Sequence Reads") +
  facet_wrap(~ Digestive_Sample, ncol = 3)+
  theme(text = element_text(family = "Helvetica", face = "bold"),
        strip.text = element_text(size = 14),
        axis.text = element_text(family = "Helvetica", face = "plain"),
        axis.text.x = element_text(size = 9),
        legend.key.size = unit(1.5, "lines"),  # Adjust legend key size
        legend.text = element_text(size = 12),  # Adjust legend text size
        legend.title = element_text(size = 14, face = "bold"))

# reorder of the samples
reordered_long$Digestive_Sample <- factor(reordered_long$Digestive_Sample, 
                                          levels = c("M1", "M4", "M5", "HP3", "HP5", "HP15", "CA14", "NTC"))

#---------- reorganizing the data frame

# breaking up the NB sections into 25-cycle and 40-cycle
for (i in 1:nrow(reordered_long)) {
  if (grepl("25_NoBlock", reordered_long$Sample[i])) {
    reordered_long$specific_blocker[i] <- "No_Blocker_25"
  }
}

# new dataframe that does not include the 25-cycle PCR
seq_df_40 <- reordered_long %>%
  filter(!grepl("25", Sample))

#--- new graph with new df
seq_df_40$Species <- factor(seq_df_40$Species, levels=c("Salmonidae_unclassified",
                                                        "Salvelinus_namaycush",
                                                        "Catostomus_commersonii",
                                                        "Homo_sapiens",
                                                        "Petromyzontidae_unclassified"))

ggplot(seq_df_40, aes(x = specific_blocker, y = Reads, fill = Species)) +
  geom_col() +
  scale_fill_manual(values = custom_colors, 
                    labels = c("Salmonid (unclassified)",  "Lake Trout", "White Sucker", "Homo Sapiens", "Sea Lamprey"))+
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
        legend.key.size = unit(1.5, "lines"),  # Adjust legend key size
        legend.text = element_text(size = 12),  # Adjust legend text size
        legend.title = element_text(size = 14, face = "bold"))


# filter out homo sapiens
seq_df_40_nohs <- filter(seq_df_40, Species != "Homo_sapiens")
# adjust colors
custom_colors2 <- c("#6bc9b2", "#40aac2","#87a7e3","#01161e")

#plot again
no_hs_seq_plot <- ggplot(seq_df_40_nohs, aes(x = specific_blocker, y = Reads, fill = Species)) +
  geom_col() +
  scale_fill_manual(values = custom_colors2, 
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
        legend.key.size = unit(1.5, "lines"),  # Adjust legend key size
        legend.text = element_text(size = 12),  # Adjust legend text size
        legend.title = element_text(size = 14, face = "bold"))

no_hs_seq_plot







#---------- more basic stats on blocking primer numbers
head(seq_df_40_nohs)

m1_nb_total <- seq_df_40_nohs %>%
  filter(specific_blocker )





















