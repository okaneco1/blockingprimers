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
matrix1$sample <- ifelse(grepl("CA14", matrix1$Sample), "CA14",
                         ifelse(grepl("HP3", matrix1$Sample), "HP3",
                                ifelse(grepl("HP15", matrix1$Sample), "HP15",
                                       ifelse(grepl("HP5", matrix1$Sample), "HP5",
                                              ifelse(grepl("M1", matrix1$Sample), "M1",
                                                     ifelse(grepl("M4", matrix1$Sample), "M4",
                                                            ifelse(grepl("M5", matrix1$Sample), "M5", "NTC")))))))

# block A or block B?
matrix1$A_B <- ifelse(grepl("_A_", matrix1$Sample), "A", 
                      ifelse(grepl("_B_", matrix1$Sample), "B","NB"))

# which end modification?
matrix1$end_mod <- ifelse(grepl("_C3_", matrix1$Sample), "C3", 
                          ifelse(grepl("_dT_", matrix1$Sample), "dT","NB"))

# HPLC purified?
matrix1$purification <- ifelse(grepl("_H_", matrix1$Sample), "H", 
                               ifelse(grepl("NoBlock", matrix1$Sample), "NB","non_H"))

# specific Sample column
matrix1$specific_blocker <- ifelse(grepl("40_PMBlock_A_C3_H_", matrix1$Sample), "40_A_C3_H",
                                   ifelse(grepl("40_PMBlock_B_C3_H_", matrix1$Sample), "40_B_C3_H",
                                          ifelse(grepl("40_PMBlock_A_dT_H_", matrix1$Sample), "40_A_dT_H",
                                                 ifelse(grepl("40_PMBlock_B_dT_H_", matrix1$Sample), "40_B_dT_H",
                                                        ifelse(grepl("40_PMBlock_A_C3", matrix1$Sample), "40_A_C3",
                                                               ifelse(grepl("40_PMBlock_B_C3", matrix1$Sample), "40_B_C3",
                                                                      ifelse(grepl("40_PMBlock_A_dT", matrix1$Sample), "40_A_dT",
                                                                             ifelse(grepl("40_PMBlock_B_dT", matrix1$Sample), "40_B_dT",
                                                                                    ifelse(grepl("25_PMBlock_B_C3_H_", matrix1$Sample), "25_B_C3_H",
                                                                                           ifelse(grepl("25_PMBlock_B_C3", matrix1$Sample), "25_B_C3", "NB"))))))))))

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




### MISSING DATA?
table(matrix1$sample)

# missing 2 CA14 and 6 NTC samples
# --- likely filtered out and not a big deal (?)





########################################################################
# comparing blocking primer effectiveness ratios
########################################################################


# make a new column with total reads 
matrix1$total_reads <- rowSums(matrix1[2:6])

# combine lake trout and Salmonidae for total salmonids (probably all lake trout)
matrix1$total_salmonids <- matrix1$Salmonidae_unclassified + matrix1$Salvelinus_namaycush

# creating a function to compare the increase in lake trout amplification per sample
# for each Sample as compared to the "no block" of that same sample. For example,
# if the "no block" version of HP3 gives x amount of lake trout reads, by how 
# much does that increase when the different blocking primers are added?

ratio_calc <- function(sample_name) {
  # start by subsampling the community matrix for just the given sample and at 40 cycles
  m2 <- filter(matrix1, sample == sample_name, cycles == "40")
  # initialize the total salmonid reads for each "no block" sample
  nb_total <- NULL
  # define the amount of salmonid reads for the "no block" sample
  for (i in 1:nrow(m2)) {
    if (m2$block[i] == "No") {
      nb_total <- m2$total_salmonids[i]
    }
  }
  # create a new column to store ratios
  m2$ratio <- 1:nrow(m2)
  # then calculate the ratio of salmonid reads for each Sample as compared to the 
  # number of salmonid reads when no block was included
  for (i in 1:nrow(m2)) {
    if (m2$block[i] == "Yes") {
      m2$ratio[i] <- m2$total_salmonids[i]/nb_total
    }
  }
  # print out the ratios
  print(select(m2, Sample, ratio))
}  

# trying with a sample
ratio_calc("HP3")
# this gives the factor by which each Sample increased the amount of salmonid reads
# so, for sample HP3, the PMBlock_A_C3 Sample increased the amount of salmonid reads 
# by a factor of 19.8.

# so, let's compare the average increase factors (ratios) for all the Samples over
# each of the samples

all_samples <- c(unique(matrix1$sample))  # organize all the sample names into an object
all_samples_noNTC <- all_samples[all_samples != "NTC"]  # exclude NTC cause no DNA

# set up a data frame to store the data
factor_df <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(factor_df) <- c("sample_name", "ratio")

# loop to run all samples
for (i in all_samples_noNTC) {
  result <- ratio_calc(i)
  factor_df <- rbind(factor_df, result)  # stores each tibble 
}

# make a single data frame with the sample names (again, just the 40 cycle samples)
matrix2 <- matrix1 %>%
  select("Sample", "specific_blocker", "sample") %>%  # only need these rows
  full_join(factor_df, by = "Sample") %>%  # link the ratio to the Sample
  filter(!grepl("NTC", sample) & !grepl("M5", sample)) %>%  # filter out the NTC samples, and M5 (WS only)
  slice(22:n())  # there were some NAs from the 25-cycle rows, so sliced those out

# with our data frame set up, plot the mean ratio (factor by which adding a blocking primer
# increases salmonid reads) for each blocking primer across the difference samples
# for which the blocking primer was included. In other words, how well, on average,
# did each blocking primer improve the amount of salmonid (lake trout) reads as 
# compared to no blocking primer?


matrix2 %>% 
  filter(!grepl("CA14", sample)) %>%  # CA14 seems to be producing some outliers, so can try without that
  ggplot(aes(x = specific_blocker, y = ratio, fill = "Same Color")) +
  geom_boxplot(alpha = 0.8) +
  geom_jitter(width = 0.2, height = 0, color = "black", alpha = 0.5) +
  geom_text_repel(aes(label = sample, alpha=0.8), vjust = 2, color = "black", size = 3, fontface = "bold", show.legend = F) +
  #scale_y_log10() +  # had log scale for previous graph that included CA14
  scale_fill_manual(values = "#69b3a2", guide = "none") +
  labs(x="Blocking Primer", y="Lake Trout Read Increase Ratio") +
  theme_classic()

# adding a mean value to the chart as well, first as a small table
mean_values <- matrix2 %>%
  filter(!grepl("CA14", sample)) %>%  # removed again for the outliers
  group_by(specific_blocker) %>%
  summarize(mean_ratio = mean(ratio))

# white dots will represent the mean
matrix2 %>% 
  filter(!grepl("CA14", sample)) %>%  # CA14 seems to be producing some outliers, so can try without that
  ggplot(aes(x = specific_blocker, y = ratio, fill = "Same Color")) +
  geom_boxplot(alpha = 0.8) +
  geom_jitter(width = 0.2, height = 0, color = "black", alpha = 0.5) +
  geom_text_repel(aes(label = sample, alpha=0.8), vjust = 2, color = "black", size = 3, fontface = "bold", show.legend = F) +
  geom_text(data = mean_values, aes(x = specific_blocker, y = mean_ratio, label = round(mean_ratio, 2)),
            vjust = -1, color = "white", size = 3, fontface = "bold") +
  geom_point(data = mean_values, aes(x = specific_blocker, y = mean_ratio),
             shape = 16, color = "white", size = 2) +
  scale_fill_manual(values = "#69b3a2", guide = "none") +
  labs(x="Blocking Primer", y="Lake Trout Read Increase Ratio") +
  theme_classic()




################################################################
# Stacked bar charts
################################################################

# using the full matrix for this, then will simplify
full_matrix <- read_excel("./12S_SL_BP_Micro.mothur.community.matrix.xls")

# change name of first column to indicate the sample
colnames(full_matrix)[1] <- "Sample_Name"

# combine similar taxonomic groups due to known feeding/concentrations
# salmonids
full_matrix$salmonids <- full_matrix$Salvelinus_namaycush + full_matrix$Salmonidae_unclassified + full_matrix$Salvelinus_unclassified + full_matrix$Salmo_unclassified
# lamprey
full_matrix$lamprey <- full_matrix$Petromyzontidae_unclassified + full_matrix$Lampetra_appendix
# sucker
full_matrix$white_sucker <- full_matrix$Catostomus_commersonii + full_matrix$Catostomidae_unclassified + full_matrix$Catostomus_unclassified
# total between these groups
full_matrix$total_grouped <- full_matrix$salmonids + full_matrix$lamprey + full_matrix$white_sucker

# Calculate proportions of each species within each sample
full_matrix <- full_matrix %>% mutate(
  lake_trout_prop = salmonids / total_grouped,
  lamprey_prop = lamprey / total_grouped,
  white_sucker_prop = white_sucker / total_grouped
)

full_matrix_long <- full_matrix %>% 
  select(Sample_Name, salmonids, lamprey, white_sucker) %>%  # Exclude the 'total' column
  pivot_longer(cols = -c(Sample_Name), names_to = "species", values_to = "reads")

# remove the ".12S" from the names
full_matrix_long$Sample_Name <- gsub(".12S","",full_matrix_long$Sample_Name)

ggplot(full_matrix_long, aes(x = Sample_Name, y = reads, fill = species)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Sample", y = "Sequence Reads", fill = "Species") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))

# reorder
cycles <- sapply(full_matrix_long$Sample_Name, FUN = function(x){strsplit(x,split = "_")[[1]][1]})
extraction <- sapply(full_matrix_long$Sample_Name, FUN = function(x){tail(strsplit(x,split = "_")[[1]],1)})

full_matrix_long2 <- full_matrix_long[order(extraction),]


pattern_order <- c("CA14", "HP3", "HP5", "HP15", "M1", "M4", "M5", "NTC")



# plot again
ggplot(full_matrix_long2, aes(x = reorder(Sample_Name, gsub(".*?(CA14|HP[0-9]+|M[0-9]+|NTC).*", "\\1", Sample_Name), FUN = function(x) match(x, pattern_order)), y = reads, fill = species)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Sample", y = "Sequence Reads", fill = "Species") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))





##### random stats ##########

matrix1 %>% 
  filter(specific_blocker == "NB" & cycles == "25") %>%
  summarize(mean = mean(Petromyzontidae_unclassified))  # 12211

matrix1 %>% 
  filter(specific_blocker == "NB" & cycles == "40") %>%
  summarize(mean = mean(Petromyzontidae_unclassified))  # 19387
  

# comparing the two cycle groups for number of lake trout reads

group25 <- subset(matrix1, cycles == "25" & specific_blocker != "NB" & sample != "NTC" & sample != "M5")
group40 <- subset(matrix1, cycles == "40" & specific_blocker != "NB" & sample != "NTC" & sample != "M5")
groupNB <- subset(matrix1, cycles == "40" & specific_blocker == "NB" & sample != "NTC" & sample != "M5")

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

reordered_long$Digestive_Sample <- ifelse(grepl("CA14", reordered_long$Sample), "CA14",
                                          ifelse(grepl("HP3", reordered_long$Sample), "HP3",
                                                 ifelse(grepl("HP15", reordered_long$Sample), "HP15",
                                                        ifelse(grepl("HP5", reordered_long$Sample), "HP5",
                                                               ifelse(grepl("M1", reordered_long$Sample), "M1",
                                                                      ifelse(grepl("M4", reordered_long$Sample), "M4",
                                                                             ifelse(grepl("M5", reordered_long$Sample), "M5", "NTC")))))))

# taking ".12S" out of the sample names
reordered_long$Sample <- gsub("\\.12S$", "", reordered_long$Sample)

#colors
custom_colors <- c("#79a4ed", "#6554AF","#001C30","#DAFFFB","#64CCC5")

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
    reordered_long$specific_blocker[i] <- "NB_25"
  }
}

# new dataframe that does not include the 25-cycle PCR
seq_df_40 <- reordered_long %>%
  filter(!grepl("25", Sample))

#--- new graph with new df
library(ggthemes)

ggplot(seq_df_40, aes(x = specific_blocker, y = Reads, fill = Species)) +
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

# ----------------------------------------------------------------
# some stats

# ANOVA

# block vs no block
aov1 <- aov(Petromyzontidae_unclassified ~ block, data = matrix1)
summary(aov)

aov2 <- aov(Petromyzontidae_unclassified ~ specific_blocker, data = matrix1)
summary(aov2)

# excluding NB samples
filtered_noNB <- matrix1 %>% filter(specific_blocker != "NB")
aov3 <- aov(Petromyzontidae_unclassified ~ specific_blocker, data = filtered_noNB)
summary(aov3)

tukey1_sl <- TukeyHSD(aov3) # no significant differences between sl reads and blockers

# looking at lake trout reads 
aov4 <- aov(Salvelinus_namaycush ~ specific_blocker, data = filtered_noNB)
summary(aov4) # no significant differences between lake trout reads and blockers

TukeyHSD(aov4)

# testing for normality
shapiro.test(filtered_noNB$Salvelinus_namaycush)
hist(filtered_noNB$Salvelinus_namaycush)

kruskal.test(Salvelinus_namaycush ~ specific_blocker, data = filtered_noNB) # p = 0.7467
wilcox.test()

# ----------------------------------------------------------------------

### average percent decrease in reads per blocking primer --- IN PROGRESS
blocker_names_list <- c("40_B_dT_H", "40_B_dT", "40_B_C3_H", "40_B_C3", "40_A_dT_H", "40_A_dT", "40_A_C3_H", "40_A_C3")

df1 <- data.frame(blocker = blocker_names_list, reads = NA)


avg_reads <- function(blocker) {
  # subset the data for each blocker
  block_data <- seq_df_40[
    seq_df_40$specific_blocker == blocker &
    seq_df_40$Species == "Petromyzontidae_unclassified" &
    seq_df_40$Digestive_Sample != "NTC",
  ]
  # subset the data for "no block"
  nb_data <- seq_df_40[
    seq_df_40$specific_blocker == "NB" &
    seq_df_40$Species == "Petromyzontidae_unclassified" &
    seq_df_40$Digestive_Sample != "NTC",
  ]
  # calculate the sum of reads for the blocker
  block_sum <- sum(block_data$Reads)
  nb_sum <- sum(nb_data$Reads)
  # calculate the mean
  mean_blockp <- 100 - ((block_sum / nb_sum) * 100)
  # calculate
  # print
  print(mean_blockp)
  print()
}

mean_block_df <- data.frame(blocker = blocker_names_list, mean_blocking_per = NA)
mean_block_df$mean_blocking_per <- sapply(mean_block_df$blocker, avg_reads)

mean(mean_block_df$mean_blocking_per)
sd(mean_block_df$mean_blocking_per)
max(mean_block_df$mean_blocking_per)
min(mean_block_df$mean_blocking_per)

write.csv(mean_block_df, "mean_block_df.csv")
