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

# change facet labels for clarity
levels(matrix1_long$Digestive_Sample) <- c("M1 (1:1, SL:LT)",
                                           "M4 (9:1, SL:LT)",
                                           "M5 (9:1, SL:WS)",
                                           "HP3 (Wild Juvenile)",
                                           "HP5 (Wild Juvenile)",
                                           "HP15 (Wild Juvenile)",
                                           "CA14 (Wild Adult)",
                                           "NTC (No-Template Control)")


#plot
matrix1_long_plot <- ggplot(matrix1_long, aes(x = specific_blocker, y = Reads, fill = Species)) +
  geom_col() +
  scale_fill_manual(values = custom_colors, 
                    labels = c("Salmonid",  "Lake Trout", "White Sucker", "Sea Lamprey"))+
  theme_few()+
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1, size = 15),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("Sequence Read Counts per OTU per Blocking Primer") +
  labs(x = "Blocking Primer", y = "Number of Sequence Reads", fill = "OTU") +
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
  filter(Digestive_Sample %in% c("M4 (9:1, SL:LT)",
                                 "M5 (9:1, SL:WS)",
                                 "HP15 (Wild Juvenile)",
                                 "CA14 (Wild Adult)"))

#drop unused levels
matrix1_simple$Digestive_Sample <- droplevels(matrix1_simple$Digestive_Sample)

# change facet labels for clarity
levels(matrix1_simple$Digestive_Sample) <- c("M4 (9:1, SL:LT)",
                                             "M5 (9:1, SL:WS)",
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





# looking at CA14 specifically
ca14 <- matrix1_long %>% 
  filter(sample == "CA14") %>%
  filter(specific_blocker == "Blocker2" | specific_blocker == "Blocker6")

blocker6_ratio <- 3063 / (3063 + 617) # lake trout proportion of total
blocker2_ratio <- 15563 / (15563 + 2819) # lake trout porportion of total


# comparison of sea lamprey blocked vs unblocked read counts

# set up data frame
read_count_comparisons <- matrix1 %>%
  filter(specific_blocker == "No_Blocker") %>%
  select(sample, Petromyzontidae_unclassified) %>%
  rename(unblocked = Petromyzontidae_unclassified)

# get blocker names
blockers <- unique(matrix1$specific_blocker)
blockers <- blockers[blockers != "No_Blocker"]

# create list to store data
blocker_list <- list()

# loop through and set read counts to blockers
for (i in 1:length(blockers)) {
  blocker_name <- paste0("blocker", i)
  blocker_df <- matrix1 %>%
    filter(specific_blocker == blockers[i]) %>%
    select(sample, Petromyzontidae_unclassified) %>%
    rename(!!blocker_name := Petromyzontidae_unclassified)
  blocker_list[[i]] <- blocker_df
}

# reduce list by sample name
all_data <- reduce(blocker_list, full_join, by = "sample")

# combine data to comparison df
final_read_count_comparisons <- full_join(read_count_comparisons, all_data, by = "sample")
final_read_count_comparisons <- filter(final_read_count_comparisons, sample != "NTC") # remove NTC row

# statistical analysis
# extract column names for blockers
blocker_columns <- names(final_read_count_comparisons)[grepl("^blocker", names(final_read_count_comparisons))]

# perform t test on each column
t_test_results <- lapply(blocker_columns, function(blocker) {
  t_test <- t.test(final_read_count_comparisons$unblocked, final_read_count_comparisons[[blocker]], paired = TRUE, na.rm = TRUE)
  data.frame(
    blocker = blocker,
    p_value = t_test$p.value,
    statistic = t_test$statistic,
    mean_unblocked = mean(final_read_count_comparisons$unblocked, na.rm = TRUE),
    mean_blocker = mean(final_read_count_comparisons[[blocker]], na.rm = TRUE)
  )
})

t_test_results_df <- do.call(rbind, t_test_results)
t_test_results_df



# comparison of lake trout blocked vs unblocked read counts
read_count_comparisons_lt <- matrix1 %>%
  filter(specific_blocker == "No_Blocker") %>%
  select(sample, Salvelinus_namaycush) %>%
  rename(unblocked = Salvelinus_namaycush)

# reset blocker list
blocker_list_lt <- list()

# loop through each blocker to assign read counts to blockers
for (i in 1:length(blockers)) {
  blocker_name <- paste0("blocker", i)
  blocker_df_lt <- matrix1 %>%
    filter(specific_blocker == blockers[i]) %>%
    select(sample, Salvelinus_namaycush) %>%
    rename(!!blocker_name := Salvelinus_namaycush)
  blocker_list_lt[[i]] <- blocker_df_lt
}

#combine
all_data_lt <- reduce(blocker_list_lt, full_join, by = "sample")

# combine data to comparison df
final_read_count_comparisons_lt <- full_join(read_count_comparisons_lt, all_data_lt, by = "sample")
final_read_count_comparisons_lt <- filter(final_read_count_comparisons_lt, sample != "NTC" & sample != "M5") # remove NTC row and M5 (white sucker)

# statistical analysis

# extract column names for blockers
blocker_columns_lt <- names(final_read_count_comparisons_lt)[grepl("^blocker", names(final_read_count_comparisons_lt))]

# perform t test on each column
t_test_results_lt <- lapply(blocker_columns_lt, function(blocker) {
  t_test <- t.test(final_read_count_comparisons_lt$unblocked, final_read_count_comparisons_lt[[blocker]], paired = TRUE, na.rm = TRUE)
  data.frame(
    blocker = blocker,
    p_value = t_test$p.value,
    statistic = t_test$statistic,
    mean_unblocked = mean(final_read_count_comparisons_lt$unblocked, na.rm = TRUE),
    mean_blocker = mean(final_read_count_comparisons_lt[[blocker]], na.rm = TRUE)
  )
})

t_test_results_lt_df <- do.call(rbind, t_test_results_lt)
print(t_test_results_df[5], row.names =FALSE)




# Check change in read counts across blocking primers and species
# This groups and summarizes the average delta in read count for all species
# across all samples when a blocking primer was included vs without a blocker

no_blocker_data <- matrix1_long %>%
  filter(specific_blocker == "No_Blocker") %>%
  select(sample, Species, Reads, Digestive_Sample) %>%
  rename(No_Blocker_Reads = Reads)

df_with_no_blocker <- matrix1_long %>%
  left_join(no_blocker_data, by = c("sample", "Species", "Digestive_Sample"))

df_with_percent_decrease <- df_with_no_blocker %>%
  mutate(Read_Count_Delta = ifelse(No_Blocker_Reads > 0, 
                                   100 * (Reads - No_Blocker_Reads) / No_Blocker_Reads, 
                                   NA))

summary_df <- df_with_percent_decrease %>%
  filter(specific_blocker != "No_Blocker") %>%  
  group_by(Digestive_Sample, Species) %>%
  summarise(Average_Read_Delta = mean(Read_Count_Delta, na.rm = TRUE))

print(summary_df)





