###'''
####Process qPCR data for pcrass analysis 
###Date: Mar 4 2024
##'''

library("tidyverse")
library(ggplot2)
library(dplyr)
library(reshape2)
library(here)
library(paletteer)
library(readxl)
library(tidyr)
library(scales)
library(here)
library(stringr)

#open raw qPCR results files
raw_results_plate_1 <- read_csv("~/Documents/240306_isolate_qPCR_plate_1_final_export.csv", col_names = TRUE)
raw_results_plate_2 <- read_csv("~/Documents/240306_isolate_qPCR_plate_2_final_export.csv", col_names = TRUE)
raw_results_plate_3 <- read_csv("~/Documents/240306_isolate_qPCR_plate_3_final_export.csv", col_names = TRUE)
raw_results_plate_4 <- read_csv("~/Documents/240306_isolate_qPCR_plate_4_final_export.csv", col_names = TRUE)
CFU_data <- read_csv("~/Documents/240306_final_isolate_CFU_data.csv", col_names = TRUE)
raw_results_rerun <- read_csv("~/Documents/240314_qPCR_liquid_colony_timepoints_rerun_isolate_teimpoints.csv", col_names = TRUE)
# Read qPCR_data where bulk ferm and isolate are in the same file
qPCR_data <- read_excel("~/Documents/231116_qPCR_results_plate_1.xlsx", sheet = "Sheet5")

####### Calculate Cq_mean, enter Std. curve, calculate copies for plate 1 #####
#only take the columns you are interested in, change Undetermined values to NA
trimmed_raw_results <- raw_results_plate_1 %>% 
  select("Sample", "Cq") %>% 
  mutate(Cq = if_else(Cq == "Undetermined", NA_character_, Cq)) %>% 
  mutate(Cq = as.numeric(Cq))

#Calculate the mean Cq value for each sample
mean_values <- trimmed_raw_results %>% 
  group_by(Sample) %>% 
  summarize(Cq_mean = mean(Cq))

#merge the raw data with the calculated mean, denote tech reps, 
#filter for samples with tech reps within 2 Cq difference
results_with_mean <- trimmed_raw_results %>% 
  left_join(mean_values, by = "Sample") %>% 
  group_by(Sample) %>% 
  mutate(tech_replicate=row_number(Sample)) %>% 
  mutate(Cq_diff = case_when(abs(combn(Cq, 2, function(x) abs(diff(x)))) > 3 ~ "Y", .default = "N")) %>% 
  filter(Cq_diff == "N") %>% 
  filter(!any(is.na(Cq)))


#Filter for the standard curve 
std_curve_plate_1 <- results_with_mean %>% 
  filter(str_detect(Sample, "std"))

#Use NEBiocalculator to look at standard curve and add the slope and intercept here: 
NEB_intercept = 35.78 #includes -2 to -9 
NEB_slope = -3.53 #includes -2 to -9
Cq_cutoff = 30
Cq_cutoff


# plot standard curve for supplementary figure#
# #Calculate the standard curve in R 
std_curve_plate_1 <- results_with_mean %>% 
  filter(str_detect(Sample, "std")) %>%
  group_by(Sample) %>%
  mutate("copies_ul" = case_when(
    Sample == "std -1" ~ 1000000000,
    Sample == "std -2" ~ 100000000,
    Sample == "std -3" ~ 10000000,
    Sample == "std -4" ~ 1000000,
    Sample == "std -5" ~ 100000,
    Sample == "std -6" ~ 10000,
    Sample == "std -7" ~ 1000,
    Sample == "std -8" ~ 100,
    Sample == "std -9" ~ 10,
  )) %>% 
  filter(!is.na(Cq_mean)) 

# Fit a linear regression model
lm_model_plate_1 <- lm(Cq_mean ~ log10(copies_ul), data = std_curve_plate_1)

# Extract coefficients and R-squared value
lm_summary_plate_1 <- summary(lm_model_plate_1)
intercept_plate_1 <- coef(lm_model_plate_1)[1]
slope_plate_1 <- coef(lm_model_plate_1)[2]
r_squared_plate_1 <- lm_summary_plate_1$r.squared

# Print intercept, slope, and R-squared value
cat("Intercept:", intercept_plate_1, "\n")
cat("Slope:", slope_plate_1, "\n")
cat("R-squared:", r_squared_plate_1, "\n")

# Plot with linear regression line
std_curve_plate_1 %>%
  ggplot(aes(x = copies_ul, y = Cq_mean)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  theme_bw() + 
  scale_x_continuous(trans = "log10") +
  labs(x = "Log10 Copies/ul", y = "Cq mean") +
  geom_text(aes(x = 1e6, y = 28, 
                label = paste("y =", round(intercept_plate_1, 2), "+", 
                              round(slope_plate_1, 2), "x", "\nR-squared =", 
                              round(r_squared_plate_1, 4))), color = "black", size = 6)


#Use the standard curve coefficients to calculate the copies per ul for samples
plate1_copies <- results_with_mean %>%
  filter(Cq_mean < Cq_cutoff) %>% 
  mutate(copies_ul = 10^((Cq_mean - NEB_intercept)/NEB_slope)) %>% #Change here if you are calculating std. curve from R
  mutate(copies_ml = copies_ul * 1000) %>% 
  separate(Sample, into = c("Sample", "Timepoint"), sep = "_T") %>% 
  separate(Sample, into = c("Sample", "bio_rep"), sep = "_(?=\\d+$)", extra = "merge", convert = TRUE) %>% 
  group_by(Sample, Timepoint) %>% 
  mutate(Avg_copies_ml = mean(copies_ml, na.rm = TRUE)) %>% 
  mutate(stdev_avg_copies_ml = sd(copies_ml, na.rm = TRUE))


####### Calculate Cq_mean, enter Std. curve, calculate copies for plate 2 #####
#only take the columns you are interested in, change Undetermined values to NA
trimmed_raw_results <- raw_results_plate_2 %>% 
  select("Sample", "Cq") %>% 
  mutate(Cq = if_else(Cq == "Undetermined", NA_character_, Cq)) %>% 
  mutate(Cq = as.numeric(Cq))

#Calculate the mean Cq value for each sample
mean_values <- trimmed_raw_results %>% 
  group_by(Sample) %>% 
  summarize(Cq_mean = mean(Cq))

#merge the raw data with the calculated mean, denote tech reps, 
#filter for samples with tech reps within 2 Cq difference
results_with_mean <- trimmed_raw_results %>% 
  left_join(mean_values, by = "Sample") %>% 
  group_by(Sample) %>% 
  mutate(tech_replicate=row_number(Sample)) %>% 
  mutate(Cq_diff = case_when(abs(combn(Cq, 2, function(x) abs(diff(x)))) > 3 ~ "Y", .default = "N")) %>% 
  filter(Cq_diff == "N") %>% 
  filter(!any(is.na(Cq)))


#Filter for the standard curve 
std_curve <- results_with_mean %>% 
  filter(str_detect(Sample, "std"))

#Use NEBiocalculator to look at standard curve and add the slope and intercept here: 
NEB_intercept = 37.41 #includes -2 to -9 
NEB_slope = -3.6 #includes -2 to -9
Cq_cutoff <- std_curve$Cq_mean[std_curve$Sample == "std -9"][1]
Cq_cutoff

# plot standard curve for supplementary figure#
# #Calculate the standard curve in R 
std_curve_plate_2 <- results_with_mean %>% 
  filter(str_detect(Sample, "std")) %>%
  group_by(Sample) %>%
  mutate("copies_ul" = case_when(
    Sample == "std -1" ~ 1000000000,
    Sample == "std -2" ~ 100000000,
    Sample == "std -3" ~ 10000000,
    Sample == "std -4" ~ 1000000,
    Sample == "std -5" ~ 100000,
    Sample == "std -6" ~ 10000,
    Sample == "std -7" ~ 1000,
    Sample == "std -8" ~ 100,
    Sample == "std -9" ~ 10,
  )) %>% 
  filter(!is.na(Cq_mean)) 

# Fit a linear regression model
lm_model_plate_2 <- lm(Cq_mean ~ log10(copies_ul), data = std_curve_plate_2)

# Extract coefficients and R-squared value
lm_summary_plate_2 <- summary(lm_model_plate_2)
intercept_plate_2 <- coef(lm_model_plate_2)[1]
slope_plate_2 <- coef(lm_model_plate_2)[2]
r_squared_plate_2 <- lm_summary_plate_2$r.squared

# Print intercept, slope, and R-squared value
cat("Intercept:", intercept_plate_2, "\n")
cat("Slope:", slope_plate_2, "\n")
cat("R-squared:", r_squared_plate_2, "\n")

# Plot with linear regression line
std_curve_plate_2 %>%
  ggplot(aes(x = copies_ul, y = Cq_mean)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  theme_bw() + 
  scale_x_continuous(trans = "log10") +
  labs(x = "Log10 Copies/ul", y = "Cq mean") +
  geom_text(aes(x = 1e6, y = 28, 
                label = paste("y =", round(intercept_plate_2, 2), "+", 
                              round(slope_plate_2, 2), "x", "\nR-squared =", 
                              round(r_squared_plate_2, 4))), color = "black", size = 6)


#Use the standard curve coefficients to calculate the copies per ul for samples
plate2_copies <- results_with_mean %>%
  filter(Cq_mean < Cq_cutoff) %>% 
  mutate(copies_ul = 10^((Cq_mean - NEB_intercept)/NEB_slope)) %>% #Change here if you are calculating std. curve from R
  mutate(copies_ml = copies_ul * 1000) %>% 
  separate(Sample, into = c("Sample", "Timepoint"), sep = "_T") %>% 
  separate(Sample, into = c("Sample", "bio_rep"), sep = "_(?=\\d+$)", extra = "merge", convert = TRUE) %>% 
  group_by(Sample, Timepoint) %>% 
  mutate(Avg_copies_ml = mean(copies_ml, na.rm = TRUE)) %>% 
  mutate(stdev_avg_copies_ml = sd(copies_ml, na.rm = TRUE))


####### Calculate Cq_mean, enter Std. curve, calculate copies for plate 3 #####
#only take the columns you are interested in, change Undetermined values to NA
trimmed_raw_results <- raw_results_plate_3 %>% 
  select("Sample", "Cq") %>% 
  mutate(Cq = if_else(Cq == "Undetermined", NA_character_, Cq)) %>% 
  mutate(Cq = as.numeric(Cq))

#Calculate the mean Cq value for each sample
mean_values <- trimmed_raw_results %>% 
  group_by(Sample) %>% 
  summarize(Cq_mean = mean(Cq))

#merge the raw data with the calculated mean, denote tech reps, 
#filter for samples with tech reps within 2 Cq difference
results_with_mean <- trimmed_raw_results %>% 
  left_join(mean_values, by = "Sample") %>% 
  group_by(Sample) %>% 
  mutate(tech_replicate=row_number(Sample)) %>% 
  mutate(Cq_diff = case_when(abs(combn(Cq, 2, function(x) abs(diff(x)))) > 3 ~ "Y", .default = "N")) %>% 
  filter(Cq_diff == "N") %>% 
  filter(!any(is.na(Cq)))


#Filter for the standard curve 
std_curve_plate_3 <- results_with_mean %>% 
  filter(str_detect(Sample, "std"))

#Use NEBiocalculator to look at standard curve and add the slope and intercept here: 
NEB_intercept = 36.69 #includes -2 to -9 
NEB_slope = -3.6 #includes -2 to -9
Cq_cutoff <- std_curve_plate_3$Cq_mean[std_curve_plate_3$Sample == "std -9"][1]
Cq_cutoff

# plot standard curve for supplementary figure#
# #Calculate the standard curve in R 
std_curve_plate_3 <- results_with_mean %>% 
  filter(str_detect(Sample, "std")) %>%
  group_by(Sample) %>%
  mutate("copies_ul" = case_when(
    Sample == "std -1" ~ 1000000000,
    Sample == "std -2" ~ 100000000,
    Sample == "std -3" ~ 10000000,
    Sample == "std -4" ~ 1000000,
    Sample == "std -5" ~ 100000,
    Sample == "std -6" ~ 10000,
    Sample == "std -7" ~ 1000,
    Sample == "std -8" ~ 100,
    Sample == "std -9" ~ 10,
  )) %>% 
  filter(!is.na(Cq_mean)) 

# Fit a linear regression model
lm_model_plate_3 <- lm(Cq_mean ~ log10(copies_ul), data = std_curve_plate_3)

# Extract coefficients and R-squared value
lm_summary_plate_3 <- summary(lm_model_plate_3)
intercept_plate_3 <- coef(lm_model_plate_3)[1]
slope_plate_3 <- coef(lm_model_plate_3)[2]
r_squared_plate_3 <- lm_summary_plate_3$r.squared

# Print intercept, slope, and R-squared value
cat("Intercept:", intercept_plate_3, "\n")
cat("Slope:", slope_plate_3, "\n")
cat("R-squared:", r_squared_plate_3, "\n")

# Plot with linear regression line
std_curve_plate_3 %>%
  ggplot(aes(x = copies_ul, y = Cq_mean)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  theme_bw() + 
  scale_x_continuous(trans = "log10") +
  labs(x = "Log10 Copies/ul", y = "Cq mean") +
  geom_text(aes(x = 1e6, y = 28, 
                label = paste("y =", round(intercept_plate_3, 2), "+", 
                              round(slope_plate_3, 2), "x", "\nR-squared =", 
                              round(r_squared_plate_3, 4))), color = "black", size = 6)


#Use the standard curve coefficients to calculate the copies per ul for samples
plate3_copies <- results_with_mean %>%
  filter(Cq_mean < Cq_cutoff) %>% 
  mutate(copies_ul = 10^((Cq_mean - NEB_intercept)/NEB_slope)) %>% #Change here if you are calculating std. curve from R
  mutate(copies_ml = copies_ul * 1000) %>% 
  separate(Sample, into = c("Sample", "Timepoint"), sep = "_T") %>% 
  separate(Sample, into = c("Sample", "bio_rep"), sep = "_(?=\\d+$)", extra = "merge", convert = TRUE) %>% 
  group_by(Sample, Timepoint) %>% 
  mutate(Avg_copies_ml = mean(copies_ml, na.rm = TRUE)) %>% 
  mutate(stdev_avg_copies_ml = sd(copies_ml, na.rm = TRUE))



####### Calculate Cq_mean, enter Std. curve, calculate copies for plate 4 #####
#only take the columns you are interested in, change Undetermined values to NA
trimmed_raw_results <- raw_results_plate_4 %>% 
  select("Sample", "Cq") %>% 
  mutate(Cq = if_else(Cq == "Undetermined", NA_character_, Cq)) %>% 
  mutate(Cq = as.numeric(Cq))

#Calculate the mean Cq value for each sample
mean_values <- trimmed_raw_results %>% 
  group_by(Sample) %>% 
  summarize(Cq_mean = mean(Cq))

#merge the raw data with the calculated mean, denote tech reps, 
#filter for samples with tech reps within 3 Cq difference
results_with_mean <- trimmed_raw_results %>% 
  left_join(mean_values, by = "Sample") %>% 
  group_by(Sample) %>% 
  mutate(tech_replicate=row_number(Sample)) %>% 
  mutate(Cq_diff = case_when(abs(combn(Cq, 2, function(x) abs(diff(x)))) > 3 ~ "Y", .default = "N")) %>% 
  filter(Cq_diff == "N") %>% 
  filter(!any(is.na(Cq)))


#Filter for the standard curve 
std_curve_plate_4 <- results_with_mean %>% 
  filter(str_detect(Sample, "std"))

#Use NEBiocalculator to look at standard curve and add the slope and intercept here: 
NEB_intercept = 36.78 #includes -2 to -9 excludes -7 because tech reps are not similar enough for our cutoff
NEB_slope = -3.47 #includes -2 to -9 excludes -7 because tech reps are not similar enough for our cutoff
Cq_cutoff <- std_curve_plate_4$Cq_mean[std_curve_plate_4$Sample == "std -9"][1]
Cq_cutoff


# plot standard curve for supplementary figure#
# #Calculate the standard curve in R 
std_curve_plate_4 <- results_with_mean %>% 
  filter(str_detect(Sample, "std")) %>%
  group_by(Sample) %>%
  mutate("copies_ul" = case_when(
    Sample == "std -1" ~ 1000000000,
    Sample == "std -2" ~ 100000000,
    Sample == "std -3" ~ 10000000,
    Sample == "std -4" ~ 1000000,
    Sample == "std -5" ~ 100000,
    Sample == "std -6" ~ 10000,
    Sample == "std -8" ~ 100,
    Sample == "std -9" ~ 10,
  )) %>% 
  filter(!is.na(Cq_mean)) 

# Fit a linear regression model
lm_model_plate_4 <- lm(Cq_mean ~ log10(copies_ul), data = std_curve_plate_4)

# Extract coefficients and R-squared value
lm_summary_plate_4 <- summary(lm_model_plate_4)
intercept_plate_4 <- coef(lm_model_plate_4)[1]
slope_plate_4 <- coef(lm_model_plate_4)[2]
r_squared_plate_4 <- lm_summary_plate_4$r.squared

# Print intercept, slope, and R-squared value
cat("Intercept:", intercept_plate_4, "\n")
cat("Slope:", slope_plate_4, "\n")
cat("R-squared:", r_squared_plate_4, "\n")

# Plot with linear regression line
std_curve_plate_4 %>%
  ggplot(aes(x = copies_ul, y = Cq_mean)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  theme_bw() + 
  scale_x_continuous(trans = "log10") +
  labs(x = "Log10 Copies/ul", y = "Cq mean") +
  geom_text(aes(x = 1e6, y = 28, 
                label = paste("y =", round(intercept_plate_4, 2), "+", 
                              round(slope_plate_4, 2), "x", "\nR-squared =", 
                              round(r_squared_plate_4, 4))), color = "black", size = 6)

#Use the standard curve coefficients to calculate the copies per ul for samples
plate4_copies <- results_with_mean %>%
  filter(Cq_mean < Cq_cutoff) %>% 
  mutate(copies_ul = 10^((Cq_mean - NEB_intercept)/NEB_slope)) %>% #Change here if you are calculating std. curve from R
  mutate(copies_ml = copies_ul * 1000) %>% 
  separate(Sample, into = c("Sample", "Timepoint"), sep = "_T") %>% 
  separate(Sample, into = c("Sample", "bio_rep"), sep = "_(?=\\d+$)", extra = "merge", convert = TRUE) %>% 
  group_by(Sample, Timepoint) %>% 
  mutate(Avg_copies_ml = mean(copies_ml, na.rm = TRUE)) %>% 
  mutate(stdev_avg_copies_ml = sd(copies_ml, na.rm = TRUE))


####### Calculate Cq_mean, enter Std. curve, calculate copies for rerun #####
#only take the columns you are interested in, change Undetermined values to NA
trimmed_raw_results <- raw_results_rerun %>% 
  select("Sample", "Cq") %>% 
  filter(Sample != 0) %>% 
  mutate(Cq = if_else(Cq == "Undetermined", NA_character_, Cq)) %>% 
  mutate(Cq = as.numeric(Cq))

#Calculate the mean Cq value for each sample
mean_values <- trimmed_raw_results %>% 
  group_by(Sample) %>% 
  summarize(Cq_mean = mean(Cq))

mean_values <- mean_values %>% filter(Sample != 0)

#merge the raw data with the calculated mean, denote tech reps, 
#filter for samples with tech reps within 2 Cq difference
results_with_mean <- trimmed_raw_results %>% 
  left_join(mean_values, by = "Sample") %>% 
  group_by(Sample) %>% 
  mutate(tech_replicate=row_number(Sample)) %>% 
  mutate(Cq_diff = case_when(abs(combn(Cq, 2, function(x) abs(diff(x)))) > 3 ~ "Y", .default = "N")) %>% 
  filter(Cq_diff == "N") %>% 
  filter(!any(is.na(Cq)))


#Filter for the standard curve 
std_curve_plate_1 <- results_with_mean %>% 
  filter(str_detect(Sample, "std"))

#Use NEBiocalculator to look at standard curve and add the slope and intercept here: 
NEB_intercept = 34.76 #includes -2 to -10 
NEB_slope = -3.41 #includes -2 to -10
Cq_cutoff <- std_curve_plate_1$Cq_mean[std_curve_plate_1$Sample == "std -10"][1]
Cq_cutoff

#Use the standard curve coefficients to calculate the copies per ul for samples
plate_rerun_copies <- results_with_mean %>%
  filter(Cq_mean < Cq_cutoff) %>% 
  mutate(copies_ul = 10^((Cq_mean - NEB_intercept)/NEB_slope)) %>% #Change here if you are calculating std. curve from R
  mutate(copies_ml = copies_ul * 1000) %>% 
  separate(Sample, into = c("Sample", "Timepoint"), sep = "_T") %>% 
  separate(Sample, into = c("Sample", "bio_rep"), sep = "_(?=\\d+$)", extra = "merge", convert = TRUE) %>% 
  group_by(Sample, Timepoint) %>% 
  mutate(Avg_copies_ml = mean(copies_ml, na.rm = TRUE)) %>% 
  mutate(stdev_avg_copies_ml = sd(copies_ml, na.rm = TRUE))




####### Aggregate the dataframes from all plates above #######

plate1_copies <- plate1_copies %>% 
  filter(!(Sample == "media_phage" & Timepoint == 0)) %>%
  filter(!(Sample == "media_phage" & Timepoint == 4)) %>%
  filter(!(Sample == "media_phage" & Timepoint == 8)) 


plate2_copies <- plate2_copies %>% 
  filter(!(Sample == "media_phage" & Timepoint == 0)) %>%
  filter(!(Sample == "media_phage" & Timepoint == 4)) %>%
  filter(!(Sample == "media_phage" & Timepoint == 8))


plate_rerun_copies <- plate_rerun_copies %>% 
  filter((Sample == "media_phage")) 

all_isolate_copies <- bind_rows(plate1_copies, plate2_copies,
                                plate3_copies, plate4_copies, 
                                plate_rerun_copies)


##### plot CFU and phage copies and qPCR replicates ####

all_qpcr_CFU_data <- all_isolate_copies %>% 
  mutate(Timepoint = as.numeric(Timepoint)) %>% 
  full_join(CFU_data)

all_qpcr_CFU_data$tech_replicate <- factor(all_qpcr_CFU_data$tech_replicate)

tech_rep_colors <- c("1" = "#009292FF", "2" = "black")

#plot supplementary figure panel E
isolate_replicates_Bv_phage <- ggplot(all_qpcr_CFU_data %>% 
                                        filter(!is.na(Cq), Sample %in% c("Bv_phage")),
                                      aes(x = interaction(Sample, Timepoint, bio_rep), y = Cq, color = tech_replicate)) +
  geom_point(size = 3, alpha = 0.5) +
  labs(x = "Sample.timepoint.biorep", y = "Cq") +
  scale_color_manual(values = tech_rep_colors) +  # Apply color scale
  theme_bw() +
  theme(text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        legend.background = element_rect(linetype = 1, 
                                         linewidth = 0.5, 
                                         color = "gray"), 
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  coord_cartesian(ylim = c(0, 30))

isolate_replicates_Bv_phage


#plot supplementary figure panel D
isolate_replicates_Bd_phage <- ggplot(all_qpcr_CFU_data %>% 
                                        filter(!is.na(Cq), Sample %in% c("Bd_phage")),
                                      aes(x = interaction(Sample, Timepoint, bio_rep), y = Cq, color = tech_replicate)) +
  geom_point(size = 3, alpha = 0.5) +
  labs(x = "Sample.timepoint.biorep", y = "Cq") +
  scale_color_manual(values = tech_rep_colors) +  # Apply color scale
  theme_bw() +
  theme(text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        legend.background = element_rect(linetype = 1, 
                                         linewidth = 0.5, 
                                         color = "gray"), 
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  coord_cartesian(ylim = c(0, 30))

isolate_replicates_Bd_phage

#plot supplementary figure panel C
isolate_replicates_Bs_phage <- ggplot(all_qpcr_CFU_data %>% 
                                        filter(!is.na(Cq), Sample %in% c("Bs_phage")),
                                      aes(x = interaction(Sample, Timepoint, bio_rep), y = Cq, color = tech_replicate)) +
  geom_point(size = 3, alpha = 0.5) +
  labs(x = "Sample.timepoint.biorep", y = "Cq") +
  scale_color_manual(values = tech_rep_colors) +  # Apply color scale
  theme_bw() +
  theme(text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        legend.background = element_rect(linetype = 1, 
                                         linewidth = 0.5, 
                                         color = "gray"), 
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  coord_cartesian(ylim = c(0, 35))

isolate_replicates_Bs_phage


#plot supplementary figure panel B
isolate_replicates_Pm_phage <- ggplot(all_qpcr_CFU_data %>% 
                                        filter(!is.na(Cq), Sample %in% c("Pm_phage")),
                                      aes(x = interaction(Sample, Timepoint, bio_rep), y = Cq, color = tech_replicate)) +
  geom_point(size = 3, alpha = 0.5) +
  labs(x = "Sample.timepoint.biorep", y = "Cq") +
  scale_color_manual(values = tech_rep_colors) +  # Apply color scale
  theme_bw() +
  theme(text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        legend.background = element_rect(linetype = 1, 
                                         linewidth = 0.5, 
                                         color = "gray"), 
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  coord_cartesian(ylim = c(0, 35))

isolate_replicates_Pm_phage


#plot supplementary figure panel A
isolate_replicates_media_phage <- ggplot(all_qpcr_CFU_data %>% 
                                           filter(!is.na(Cq), Sample %in% c("media_phage")),
                                         aes(x = interaction(Sample, Timepoint, bio_rep), y = Cq, color = tech_replicate)) +
  geom_point(size = 3, alpha = 0.5) +
  labs(x = "Sample.timepoint.biorep", y = "Cq") +
  scale_color_manual(values = tech_rep_colors) +  # Apply color scale
  theme_bw() +
  theme(text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        legend.background = element_rect(linetype = 1, 
                                         linewidth = 0.5, 
                                         color = "gray"), 
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  coord_cartesian(ylim = c(0, 35))

isolate_replicates_media_phage



for_plotting_df <- all_qpcr_CFU_data %>% 
  select(Sample, Timepoint, Avg_copies_ml, stdev_avg_copies_ml, CFU_ml, CFU_stdev) %>% 
  group_by(Sample, Timepoint) %>% 
  distinct(Sample, .keep_all = TRUE) %>% 
  ungroup() %>%
  filter(stdev_avg_copies_ml != 0 | is.na(stdev_avg_copies_ml))


#make plot for figure 7 Panel A
Bv_phage_and_CFU <- ggplot(for_plotting_df %>% 
                             filter((Sample %in% c("Bv", "media_phage") | (Sample == "Bv_phage" & !is.na(CFU_ml))) & Timepoint != 4), 
                           aes(x = Timepoint, group = Sample)) +
  geom_point(aes(y = Avg_copies_ml, color = Sample), shape = 16, size = 4) +
  geom_line(aes(y = Avg_copies_ml, color = Sample), size = 1.5) +
  geom_point(aes(y = CFU_ml, color = Sample), shape = 16, size = 4) +
  geom_line(aes(y = CFU_ml, color = Sample), size = 1.5, linetype = "dashed") +
  scale_color_manual(values = c("Bv_phage" = "#009292FF", 
                                "Bv" = "black", 
                                "media_phage" = "grey"),
                     name = "Sample", 
                     labels = c("Bv_phage" = "B. vulgauts + phage", 
                                "Bv" = "B. vulgatus - phage", 
                                "media_phage" = "Media + phage")) +
  labs(title = "", x = "Time (Hours)", y = "CFU or phage per mL") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        text = element_text(size = 40),
        panel.grid.major = element_line(linewidth = 1.5),
        panel.grid.minor = element_line(linewidth = 1.5),
        axis.line = element_line(linetype = 1, linewidth = 1),
        panel.border = element_rect(linetype = 1, linewidth = 3),
        legend.background = element_rect(linetype = 1, linewidth = 0.5, color = "gray")) +
  scale_y_log10() +
  coord_cartesian(ylim = c(1e4, 1e10)) + 
  geom_errorbar(aes(ymin = Avg_copies_ml - stdev_avg_copies_ml, 
                    ymax = Avg_copies_ml + stdev_avg_copies_ml, color = Sample), width = 0.5) +
  geom_errorbar(aes(ymin = CFU_ml - CFU_stdev, ymax = CFU_ml + CFU_stdev, 
                    color = Sample), width = 0.5)

Bv_phage_and_CFU

pdf(here("~/Documents/240307_Bv_phage_CFU.pdf"), width = 14, height = 14)
print(Bv_phage_and_CFU)
dev.off()

#make plot for figure 7 Panel B
Bd_phage_and_CFU <- ggplot(for_plotting_df %>% filter(Sample %in% c("Bd_phage", "Bd", "media_phage")), 
                           aes(x = Timepoint)) +
  geom_point(aes(y = Avg_copies_ml, color = Sample), shape = 16, size = 4) +
  geom_line(aes(y = Avg_copies_ml, color = Sample), size = 1.5) +
  geom_point(aes(y = CFU_ml, color = Sample), shape = 16, size = 4) +
  geom_line(aes(y = CFU_ml, color = Sample), size = 1.5, linetype = "dashed") +
  scale_color_manual(values = c("Bd_phage" = "#009292FF", "Bd" = "black", "media_phage" = "grey"),
                     name = "Sample", labels = c("Bd_phage" = "B. dorei + phage", "Bd" = "B. dorei - phage")) +
  labs(title = "", x = "Time (Hours)", y = "CFU or phage per mL") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        text = element_text(size = 40),
        panel.grid.major = element_line(linewidth = 1.5),
        panel.grid.minor = element_line(linewidth = 1.5),
        axis.line = element_line(linetype = 1, linewidth = 1),
        panel.border = element_rect(linetype = 1, linewidth = 3),
        legend.background = element_rect(linetype = 1, linewidth = 0.5, color = "gray")) +
  scale_y_log10() +
  coord_cartesian(ylim = c(1e4, 1e10)) + 
  geom_errorbar(aes(ymin = Avg_copies_ml - stdev_avg_copies_ml, 
                    ymax = Avg_copies_ml + stdev_avg_copies_ml, color = Sample), width = 0.5) +
  geom_errorbar(aes(ymin = CFU_ml - CFU_stdev, ymax = CFU_ml + CFU_stdev, 
                    color = Sample), width = 0.5)

Bd_phage_and_CFU
pdf(here("~/Documents/240307_Bd_phage_CFU.pdf"), width = 14, height = 14)
print(Bd_phage_and_CFU)
dev.off()


#supplementary panel E
combined_phage <- ggplot(for_plotting_df %>% filter(Sample %in% c("Bs_phage", "Pm_phage", "media_phage")), 
                         aes(x = Timepoint)) +
  geom_point(aes(y = Avg_copies_ml, color = Sample), shape = 16, size = 4) +
  geom_line(aes(y = Avg_copies_ml, color = Sample), size = 1.5) +
  geom_point(aes(y = CFU_ml, color = Sample), shape = 16, size = 4) +
  geom_line(aes(y = CFU_ml, color = Sample), size = 1.5, linetype = "dashed") +
  scale_color_manual(values = c("Bs_phage" = "#009292FF", "Pm_phage" = "purple", "media_phage" = "grey"),
                     name = "Sample", labels = c("Bs_phage" = "B. stercoris + phage", 
                                                 "Pm_phage" = "P.merdae + phage",
                                                 "media_phage" = "media + phage")) +
  labs(title = "", x = "Time (Hours)", y = "Phage per mL") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        text = element_text(size = 40),
        panel.grid.major = element_line(linewidth = 1.5),
        panel.grid.minor = element_line(linewidth = 1.5),
        axis.line = element_line(linetype = 1, linewidth = 1),
        panel.border = element_rect(linetype = 1, linewidth = 3),
        legend.background = element_rect(linetype = 1, linewidth = 0.5, color = "gray")) +
  scale_y_log10() +
  geom_errorbar(aes(ymin = Avg_copies_ml - stdev_avg_copies_ml, ymax = Avg_copies_ml + stdev_avg_copies_ml, color = Sample), width = 0.5) +
  geom_errorbar(aes(ymin = CFU_ml - CFU_stdev, ymax = CFU_ml + CFU_stdev, color = Sample), width = 0.5) # Applying log scale to y-axis

combined_phage
pdf(here("~/Documents/240307_Bs_Pm_media.pdf"), width = 14, height = 14)
print(combined_phage)
dev.off()

######## perform t tests on CFU values ########
#function from stack overflow to do t test 

# m1, m2: the sample means
# s1, s2: the sample standard deviations
# n1, n2: the same sizes
# m0: the null value for the difference in means to be tested for. Default is 0. 
# equal.variance: whether or not to assume equal variance. Default is FALSE. 
t.test2 <- function(m1,m2,s1,s2,n1,n2,m0=0,equal.variance=FALSE)
{
  if( equal.variance==FALSE ) 
  {
    se <- sqrt( (s1^2/n1) + (s2^2/n2) )
    # welch-satterthwaite df
    df <- ( (s1^2/n1 + s2^2/n2)^2 )/( (s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1) )
  } else
  {
    # pooled standard deviation, scaled by the sample sizes
    se <- sqrt( (1/n1 + 1/n2) * ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2) ) 
    df <- n1+n2-2
  }      
  t <- (m1-m2-m0)/se 
  dat <- c(m1-m2, se, t, 2*pt(-abs(t),df))    
  names(dat) <- c("Difference of means", "Std Error", "t", "p-value")
  return(dat) 
}

##For loop to iterate over all timepoints and calculate p-value for CFU with or without phage

### this runs for Bd and Bd phage
# Initialize an empty list to store the t-test results
t_test_results <- list()

# Iterate over each timepoint
for (tp in c(0, 4, 8, 16, 20, 28, 32, 40, 44)) {
  # Filter the data for the current timepoint and samples
  filtered_data <- CFU_data %>% 
    filter(Sample %in% c("Bd_phage", "Bd") & Timepoint == tp)
  
  filtered_data$CFU_ml[1]
  filtered_data$CFU_ml[2]
  # Perform the t-test
  t_test <- t.test2(filtered_data$CFU_ml[1], filtered_data$CFU_ml[2], filtered_data$CFU_stdev[1], filtered_data$CFU_stdev[2], 3, 3)
  
  
  # Store the t-test result
  t_test_results[[as.character(tp)]] <- t_test
}

# Access and print the results
for (tp in names(t_test_results)) {
  cat("Timepoint:", tp, "\n")
  print(t_test_results[[tp]])
  cat("\n")
}

### this runs for Bv and Bv phage 
# Initialize an empty list to store the t-test results
t_test_results_Bv <- list()

# Iterate over each timepoint
for (tp in c(0, 4, 8, 12, 16, 20, 32, 36, 40, 44)) {
  # Filter the data for the current timepoint and samples
  filtered_data <- CFU_data %>% 
    filter(Sample %in% c("Bv_phage", "Bv") & Timepoint == tp)
  
  filtered_data$CFU_ml[1]
  filtered_data$CFU_ml[2]
  # Perform the t-test
  t_test_Bv <- t.test2(filtered_data$CFU_ml[1], filtered_data$CFU_ml[2], filtered_data$CFU_stdev[1], filtered_data$CFU_stdev[2], 3, 3)
  
  # Store the t-test result
  t_test_results_Bv[[as.character(tp)]] <- t_test_Bv
}

# Access and print the results
for (tp in names(t_test_results_Bv)) {
  cat("Timepoint:", tp, "\n")
  print(t_test_results_Bv[[tp]])
  cat("\n")
}

# m1, m2: the sample means
# s1, s2: the sample standard deviations
# n1, n2: the same sizes
# m0: the null value for the difference in means to be tested for. Default is 0. 
# equal.variance: whether or not to assume equal variance. Default is FALSE. 
t.test2 <- function(m1,m2,s1,s2,n1,n2,m0=0,equal.variance=FALSE)
{
  if( equal.variance==FALSE ) 
  {
    se <- sqrt( (s1^2/n1) + (s2^2/n2) )
    # welch-satterthwaite df
    df <- ( (s1^2/n1 + s2^2/n2)^2 )/( (s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1) )
  } else
  {
    # pooled standard deviation, scaled by the sample sizes
    se <- sqrt( (1/n1 + 1/n2) * ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2) ) 
    df <- n1+n2-2
  }      
  t <- (m1-m2-m0)/se 
  dat <- c(m1-m2, se, t, 2*pt(-abs(t),df))    
  names(dat) <- c("Difference of means", "Std Error", "t", "p-value")
  return(dat) 
}

######plot phage:host #####
##### plot phage:host ratio###


# plopt panel C for figure 7
phage_host_ratio_plot <- ggplot(qPCR_data %>% filter(sample == "Bd_phage_1" | sample == "Bv_phage_1" | sample == "bulk_ferm"), 
                              aes(x = time, y = phage_per_CFU)) +
  geom_point(aes(color = sample, shape = sample), size = 4) +
  scale_color_manual(values = c("Bd_phage_1" = "#009292FF", "Bv_phage_1" = "black", "bulk_ferm" = "honeydew4")) +
  labs(title = expression(paste("")),
       x = "Time (Hours)", y = "C. communis copies:Bacteria (CFU or copies)") +
  geom_line(aes(color = sample), size = 1.1) +
  geom_hline(yintercept = 10, linetype = "dashed", color = "grey59", size = 1) +  # Add horizontal line at y = 10
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey59", size = 1) +
  theme_bw() +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    panel.grid.major = element_line(linewidth = 1),
    panel.grid.minor = element_line(linewidth = 1),
    axis.line = element_line(linetype = 1, linewidth = 1),
    panel.border = element_rect(linetype = 1, linewidth = 2),
    legend.background = element_rect(linetype = 1, linewidth = 0.5, color = "gray")
  ) +
  scale_y_log10(expand = expand_scale(mult = c(0.01, 0.1)),
                labels = scales::label_number())

phage_host_ratio_plot





########plot supernatant vs. pellet##########
plate_1 <- read_csv("~/Documents/240709_T0_T12_sup_v_pel_nanopore.csv", col_names = TRUE)
plate_2 <- read_csv("~/Documents/240709_plate_2_results_T24_T36.csv", col_names = TRUE)
plate_3 <- read_csv("~/Documents/240709_plate_3_T44.csv", col_names = TRUE)
plate_4 <- read_csv("~/Documents/240709_T12_sup_pel_results.csv", col_names = TRUE)

##Calculate Cq_mean, enter Std. curve, calculate copies for plate 1 communis copies
#only take the columns you are interested in, change Undetermined values to NA

trimmed_raw_results_P1 <- plate_1 %>% 
  select("Sample", "Cq") %>% 
  filter(Sample != 0) %>% 
  mutate(Cq = if_else(Cq == "Undetermined", NA_character_, Cq)) %>% 
  mutate(Cq = as.numeric(Cq))

#Calculate the mean Cq value for each sample
mean_values_P1 <- trimmed_raw_results_P1 %>% 
  group_by(Sample) %>% 
  summarize(Cq_mean = mean(Cq))

mean_values_P1 <- mean_values_P1 %>% filter(Sample != 0)

#merge the raw data with the calculated mean, denote tech reps, 
#filter for samples with tech reps within 2 Cq difference
results_with_mean_P1 <- trimmed_raw_results_P1 %>% 
  left_join(mean_values_P1, by = "Sample") %>% 
  group_by(Sample) %>% 
  mutate(tech_replicate=row_number(Sample)) %>% 
  mutate(Cq_diff = case_when(abs(combn(Cq, 2, function(x) abs(diff(x)))) > 3 ~ "Y", .default = "N")) %>% 
  filter(Cq_diff == "N") %>% 
  filter(!any(is.na(Cq)))


#Filter for the standard curve 
std_curve_plate_1 <- results_with_mean_P1 %>% 
  filter(str_detect(Sample, "std"))

#Use NEBiocalculator to look at standard curve and add the slope and intercept here: 
NEB_intercept = 37.74 #includes -2 to -9 
NEB_slope = -3.54 #includes -2 to -9
Cq_cutoff <- 32.704367
Cq_cutoff


# plot standard curve for supplementary figure#
# #Calculate the standard curve in R 
std_curve_plate_1 <- results_with_mean_P1 %>% 
  filter(str_detect(Sample, "std") & !Sample %in% c("std_-10", "std_-11", "std_-12")) %>%
  group_by(Sample) %>%
  mutate("copies_ul" = case_when(
    Sample == "std_-2" ~ 100000000,
    Sample == "std_-3" ~ 10000000,
    Sample == "std_-4" ~ 1000000,
    Sample == "std_-5" ~ 100000,
    Sample == "std_-6" ~ 10000,
    Sample == "std_-7" ~ 1000,
    Sample == "std_-8" ~ 100,
    Sample == "std_-9" ~ 10
  )) %>% 
  filter(!is.na(Cq_mean)) 



# Fit a linear regression model
lm_model_plate_1 <- lm(Cq_mean ~ log10(copies_ul), data = std_curve_plate_1)

# Extract coefficients and R-squared value
lm_summary_plate_1 <- summary(lm_model_plate_1)
intercept_plate_1 <- coef(lm_model_plate_1)[1]
slope_plate_1 <- coef(lm_model_plate_1)[2]
r_squared_plate_1 <- lm_summary_plate_1$r.squared

# Print intercept, slope, and R-squared value
cat("Intercept:", intercept_plate_1, "\n")
cat("Slope:", slope_plate_1, "\n")
cat("R-squared:", r_squared_plate_1, "\n")

# Plot with linear regression line
std_curve_plate_1 %>%
  ggplot(aes(x = copies_ul, y = Cq_mean)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  theme_bw() + 
  scale_x_continuous(trans = "log10") +
  labs(x = "Log10 Copies/ul", y = "Cq mean") +
  geom_text(aes(x = 1e6, y = 28, 
                label = paste("y =", round(intercept_plate_1, 2), "+", 
                              round(slope_plate_1, 2), "x", "\nR-squared =", 
                              round(r_squared_plate_1, 4))), color = "black", size = 6)


#Use the standard curve coefficients to calculate the copies per ul for samples
plate1_copies <- results_with_mean_P1 %>%
  filter(Cq_mean < Cq_cutoff) %>%
  mutate(copies_ul = 10^((Cq_mean - NEB_intercept)/NEB_slope)) %>%
  mutate(copies_ml = copies_ul * 1000) %>%
  separate(Sample, into = c("Sample", "Timepoint"), sep = "_T") %>%
  separate(Sample, into = c("Sample", "bio_rep"), sep = "_(?=\\d+$)", extra = "merge", convert = TRUE) %>%
  group_by(Sample, Timepoint) %>%
  mutate(Avg_copies_ml = mean(copies_ml, na.rm = TRUE)) %>%
  mutate(stdev_avg_copies_ml = sd(copies_ml, na.rm = TRUE))




###Calculate Cq_mean, enter Std. curve, calculate copies for plate 2
#only take the columns you are interested in, change Undetermined values to NA
trimmed_raw_results_P2 <- plate_2 %>% 
  select("Sample", "Cq") %>% 
  filter(Sample != 0) %>% 
  mutate(Cq = if_else(Cq == "Undetermined", NA_character_, Cq)) %>% 
  mutate(Cq = as.numeric(Cq))

#Calculate the mean Cq value for each sample
mean_values_P2 <- trimmed_raw_results_P2 %>% 
  group_by(Sample) %>% 
  summarize(Cq_mean = mean(Cq))

mean_values_P2 <- mean_values_P2 %>% filter(Sample != 0)

#merge the raw data with the calculated mean, denote tech reps, 
#filter for samples with tech reps within 2 Cq difference
results_with_mean_P2 <- trimmed_raw_results_P2 %>% 
  left_join(mean_values_P2, by = "Sample") %>% 
  group_by(Sample) %>% 
  mutate(tech_replicate=row_number(Sample)) %>% 
  mutate(Cq_diff = case_when(abs(combn(Cq, 2, function(x) abs(diff(x)))) > 3 ~ "Y", .default = "N")) %>% 
  filter(Cq_diff == "N") %>% 
  filter(!any(is.na(Cq)))


#Filter for the standard curve 
std_curve_plate_2 <- results_with_mean_P2 %>% 
  filter(str_detect(Sample, "std"))

#Use NEBiocalculator to look at standard curve and add the slope and intercept here: 
NEB_intercept = 40.65 #includes -2 to -9 
NEB_slope = -3.97 #includes -2 to -9
Cq_cutoff <- 32.514924
Cq_cutoff


# plot standard curve for supplementary figure#
# #Calculate the standard curve in R 
std_curve_plate_2 <- results_with_mean_P2 %>% 
  filter(str_detect(Sample, "std") & !Sample %in% c("std_-10", "std_-11", "std_-12")) %>%
  group_by(Sample) %>%
  mutate("copies_ul" = case_when(
    Sample == "std_-2" ~ 100000000,
    Sample == "std_-3" ~ 10000000,
    Sample == "std_-4" ~ 1000000,
    Sample == "std_-5" ~ 100000,
    Sample == "std_-6" ~ 10000,
    Sample == "std_-7" ~ 1000,
    Sample == "std_-8" ~ 100,
    Sample == "std_-9" ~ 10
  )) %>% 
  filter(!is.na(Cq_mean)) 

# Fit a linear regression model
lm_model_plate_2 <- lm(Cq_mean ~ log10(copies_ul), data = std_curve_plate_2)

# Extract coefficients and R-squared value
lm_summary_plate_2 <- summary(lm_model_plate_2)
intercept_plate_2 <- coef(lm_model_plate_2)[1]
slope_plate_2 <- coef(lm_model_plate_2)[2]
r_squared_plate_2 <- lm_summary_plate_2$r.squared

# Print intercept, slope, and R-squared value
cat("Intercept:", intercept_plate_2, "\n")
cat("Slope:", slope_plate_2, "\n")
cat("R-squared:", r_squared_plate_2, "\n")

# Plot with linear regression line
std_curve_plate_2 %>%
  ggplot(aes(x = copies_ul, y = Cq_mean)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  theme_bw() + 
  scale_x_continuous(trans = "log10") +
  labs(x = "Log10 Copies/ul", y = "Cq mean") +
  geom_text(aes(x = 1e6, y = 28, 
                label = paste("y =", round(intercept_plate_2, 2), "+", 
                              round(slope_plate_2, 2), "x", "\nR-squared =", 
                              round(r_squared_plate_2, 4))), color = "black", size = 6)


#Use the standard curve coefficients to calculate the copies per ul for samples
plate2_copies <- results_with_mean_P2 %>%
  filter(Cq_mean < Cq_cutoff) %>%
  mutate(copies_ul = 10^((Cq_mean - NEB_intercept)/NEB_slope)) %>%
  mutate(copies_ml = copies_ul * 1000) %>%
  separate(Sample, into = c("Sample", "Timepoint"), sep = "_T") %>%
  separate(Sample, into = c("Sample", "bio_rep"), sep = "_(?=\\d+$)", extra = "merge", convert = TRUE) %>%
  group_by(Sample, Timepoint) %>%
  mutate(Avg_copies_ml = mean(copies_ml, na.rm = TRUE)) %>%
  mutate(stdev_avg_copies_ml = sd(copies_ml, na.rm = TRUE))



####### Calculate Cq_mean, enter Std. curve, calculate copies for plate 3 
#only take the columns you are interested in, change Undetermined values to NA
trimmed_raw_results_P3 <- plate_3 %>% 
  select("Sample", "Cq") %>% 
  filter(Sample != 0) %>% 
  mutate(Cq = if_else(Cq == "Undetermined", NA_character_, Cq)) %>% 
  mutate(Cq = as.numeric(Cq))

#Calculate the mean Cq value for each sample
mean_values_P3 <- trimmed_raw_results_P3 %>% 
  group_by(Sample) %>% 
  summarize(Cq_mean = mean(Cq))

mean_values_P3 <- mean_values_P3 %>% filter(Sample != 0)

#merge the raw data with the calculated mean, denote tech reps, 
#filter for samples with tech reps within 2 Cq difference
results_with_mean_P3 <- trimmed_raw_results_P3 %>% 
  left_join(mean_values_P3, by = "Sample") %>% 
  group_by(Sample) %>% 
  mutate(tech_replicate=row_number(Sample)) %>% 
  mutate(Cq_diff = case_when(abs(combn(Cq, 2, function(x) abs(diff(x)))) > 4 ~ "Y", .default = "N")) %>% 
  filter(Cq_diff == "N") %>% 
  filter(!any(is.na(Cq)))


#Filter for the standard curve 
std_curve_plate_3 <- results_with_mean_P3 %>% 
  filter(str_detect(Sample, "std"))

#Use NEBiocalculator to look at standard curve and add the slope and intercept here: 
NEB_intercept = 36.97 #includes -2 to -7 and exclused 6 because tech reps are not similar enough for our cutoff
NEB_slope = -3.48 #includes -2 to -7
Cq_cutoff <- 32.378567
Cq_cutoff

# plot standard curve for supplementary figure#
# #Calculate the standard curve in R 
std_curve_plate_3 <- results_with_mean_P3 %>% 
  filter(str_detect(Sample, "std") & !Sample %in% c("std_-6", "std_-8", "std_-9", "std_-10", "std_-11", "std_-12")) %>%
  group_by(Sample) %>%
  mutate("copies_ul" = case_when(
    Sample == "std_-2" ~ 100000000,
    Sample == "std_-3" ~ 10000000,
    Sample == "std_-4" ~ 1000000,
    Sample == "std_-5" ~ 100000,
    Sample == "std_-7" ~ 1000,
  )) %>% 
  filter(!is.na(Cq_mean)) 

# Fit a linear regression model
lm_model_plate_3 <- lm(Cq_mean ~ log10(copies_ul), data = std_curve_plate_3)

# Extract coefficients and R-squared value
lm_summary_plate_3 <- summary(lm_model_plate_3)
intercept_plate_3 <- coef(lm_model_plate_3)[1]
slope_plate_3 <- coef(lm_model_plate_3)[2]
r_squared_plate_3 <- lm_summary_plate_3$r.squared

# Print intercept, slope, and R-squared value
cat("Intercept:", intercept_plate_3, "\n")
cat("Slope:", slope_plate_3, "\n")
cat("R-squared:", r_squared_plate_3, "\n")

# Plot with linear regression line
std_curve_plate_3 %>%
  ggplot(aes(x = copies_ul, y = Cq_mean)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  theme_bw() + 
  scale_x_continuous(trans = "log10") +
  labs(x = "Log10 Copies/ul", y = "Cq mean") +
  geom_text(aes(x = 1e6, y = 28, 
                label = paste("y =", round(intercept_plate_3, 2), "+", 
                              round(slope_plate_3, 2), "x", "\nR-squared =", 
                              round(r_squared_plate_3, 4))), color = "black", size = 6)


#Use the standard curve coefficients to calculate the copies per ul for samples
plate3_copies <- results_with_mean_P3 %>%
  filter(Cq_mean < Cq_cutoff) %>%
  mutate(copies_ul = 10^((Cq_mean - NEB_intercept)/NEB_slope)) %>%
  mutate(copies_ml = copies_ul * 1000) %>%
  separate(Sample, into = c("Sample", "Timepoint"), sep = "_T") %>%
  separate(Sample, into = c("Sample", "bio_rep"), sep = "_(?=\\d+$)", extra = "merge", convert = TRUE) %>%
  group_by(Sample, Timepoint) %>%
  mutate(Avg_copies_ml = mean(copies_ml, na.rm = TRUE)) %>%
  mutate(stdev_avg_copies_ml = sd(copies_ml, na.rm = TRUE))


#Use the standard curve coefficients to calculate the copies per ul for samples
plate3_copies <- results_with_mean_P3 %>%
  filter(Cq_mean < Cq_cutoff) %>%
  mutate(copies_ul = 10^((Cq_mean - NEB_intercept)/NEB_slope)) %>%
  mutate(copies_ml = copies_ul * 1000) %>%
  separate(Sample, into = c("Sample", "Timepoint"), sep = "_T") %>%
  separate(Sample, into = c("Sample", "bio_rep"), sep = "_(?=\\d+$)", extra = "merge", convert = TRUE) %>%
  group_by(Sample, Timepoint) %>%
  mutate(Avg_copies_ml = mean(copies_ml, na.rm = TRUE)) %>%
  mutate(stdev_avg_copies_ml = sd(copies_ml, na.rm = TRUE))


####### Calculate Cq_mean, enter Std. curve, calculate copies for plate 4 communis copies 
#only take the columns you are interested in, change Undetermined values to NA

trimmed_raw_results_P4 <- plate_4 %>% 
  select("Sample", "Cq") %>% 
  filter(Sample != 0) %>% 
  mutate(Cq = if_else(Cq == "Undetermined", NA_character_, Cq)) %>% 
  mutate(Cq = as.numeric(Cq))

#Calculate the mean Cq value for each sample
mean_values_P4 <- trimmed_raw_results_P4 %>% 
  group_by(Sample) %>% 
  summarize(Cq_mean = mean(Cq))

mean_values_P4 <- mean_values_P4 %>% filter(Sample != 0)

#merge the raw data with the calculated mean, denote tech reps, 
#filter for samples with tech reps within 2 Cq difference
results_with_mean_P4 <- trimmed_raw_results_P4 %>% 
  left_join(mean_values_P4, by = "Sample") %>% 
  group_by(Sample) %>% 
  mutate(tech_replicate=row_number(Sample)) %>% 
  mutate(Cq_diff = case_when(abs(combn(Cq, 2, function(x) abs(diff(x)))) > 3 ~ "Y", .default = "N")) %>% 
  filter(Cq_diff == "N") %>% 
  filter(!any(is.na(Cq)))


#Filter for the standard curve 
std_curve_plate_4 <- results_with_mean_P4 %>% 
  filter(str_detect(Sample, "std"))

#Use NEBiocalculator to look at standard curve and add the slope and intercept here: 
NEB_intercept = 41.18 #includes -2 to -9 
NEB_slope = -3.93 #includes -2 to -9
Cq_cutoff <- 36.426077
Cq_cutoff

# plot standard curve for supplementary figure#
# #Calculate the standard curve in R 
std_curve_plate_4 <- results_with_mean_P4 %>% 
  filter(str_detect(Sample, "std") & !Sample %in% c("std_-10", "std_-11", "std_-12")) %>%
  group_by(Sample) %>%
  mutate("copies_ul" = case_when(
    Sample == "std_-2" ~ 100000000,
    Sample == "std_-3" ~ 10000000,
    Sample == "std_-4" ~ 1000000,
    Sample == "std_-5" ~ 100000,
    Sample == "std_-7" ~ 1000,
    Sample == "std_-8" ~ 100,
    Sample == "std_-9" ~ 10,
  )) %>% 
  filter(!is.na(Cq_mean)) 

# Fit a linear regression model
lm_model_plate_4 <- lm(Cq_mean ~ log10(copies_ul), data = std_curve_plate_4)

# Extract coefficients and R-squared value
lm_summary_plate_4 <- summary(lm_model_plate_4)
intercept_plate_4 <- coef(lm_model_plate_4)[1]
slope_plate_4 <- coef(lm_model_plate_4)[2]
r_squared_plate_4 <- lm_summary_plate_4$r.squared

# Print intercept, slope, and R-squared value
cat("Intercept:", intercept_plate_4, "\n")
cat("Slope:", slope_plate_4, "\n")
cat("R-squared:", r_squared_plate_4, "\n")

# Plot with linear regression line
std_curve_plate_4 %>%
  ggplot(aes(x = copies_ul, y = Cq_mean)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  theme_bw() + 
  scale_x_continuous(trans = "log10") +
  labs(x = "Log10 Copies/ul", y = "Cq mean") +
  geom_text(aes(x = 1e6, y = 28, 
                label = paste("y =", round(intercept_plate_4, 2), "+", 
                              round(slope_plate_4, 2), "x", "\nR-squared =", 
                              round(r_squared_plate_4, 4))), color = "black", size = 6)


#Use the standard curve coefficients to calculate the copies per ul for samples
plate4_copies <- results_with_mean_P4 %>%
  filter(Cq_mean < Cq_cutoff) %>%
  mutate(copies_ul = 10^((Cq_mean - NEB_intercept)/NEB_slope)) %>%
  mutate(copies_ml = copies_ul * 1000) %>%
  separate(Sample, into = c("Sample", "Timepoint"), sep = "_T") %>%
  separate(Sample, into = c("Sample", "bio_rep"), sep = "_(?=\\d+$)", extra = "merge", convert = TRUE) %>%
  group_by(Sample, Timepoint) %>%
  mutate(Avg_copies_ml = mean(copies_ml, na.rm = TRUE)) %>%
  mutate(stdev_avg_copies_ml = sd(copies_ml, na.rm = TRUE))


combined_copies_1 <- bind_rows(plate4_copies, plate3_copies, plate2_copies, plate1_copies) %>%
  group_by(Sample, Timepoint) %>% 
  filter(n_distinct(bio_rep) > 1) %>% 
  ungroup()


combined_copies <- combined_copies_1 %>%
  mutate(Timepoint = as.numeric(Timepoint))
combined_copies <- combined_copies %>% filter(Timepoint != 12)

combined_copies$tech_replicate <- factor(combined_copies$tech_replicate)

tech_rep_colors <- c("1" = "#009292FF", "2" = "black")

#plot supplementary figure 11G
Bd_P_sup_v_pel_replicates <- ggplot(combined_copies %>% 
                                        filter(!is.na(Cq), Sample %in% c("Bd+P_sup_v_pel", "Bd+P_sup", "Bd+P_pel")),
                                      aes(x = interaction(Sample, Timepoint, bio_rep), y = Cq, color = tech_replicate)) +
  geom_point(size = 3, alpha = 0.5) +
  labs(x = "Sample.timepoint.biorep", y = "Cq") +
  scale_color_manual(values = tech_rep_colors) +  # Apply color scale
  theme_bw() +
  theme(text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        legend.background = element_rect(linetype = 1, 
                                         linewidth = 0.5, 
                                         color = "gray"), 
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  coord_cartesian(ylim = c(0, 35))

Bd_P_sup_v_pel_replicates


#plot supplementary figure 11H
Bv_P_sup_v_pel_replicates <- ggplot(combined_copies %>% 
                                      filter(!is.na(Cq), Sample %in% c("Bv+P_sup_v_pel", "Bv+P_sup", "Bv+P_pel")),
                                    aes(x = interaction(Sample, Timepoint, bio_rep), y = Cq, color = tech_replicate)) +
  geom_point(size = 3, alpha = 0.5) +
  labs(x = "Sample.timepoint.biorep", y = "Cq") +
  scale_color_manual(values = tech_rep_colors) +  # Apply color scale
  theme_bw() +
  theme(text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        legend.background = element_rect(linetype = 1, 
                                         linewidth = 0.5, 
                                         color = "gray"), 
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  coord_cartesian(ylim = c(0, 35))

Bv_P_sup_v_pel_replicates

#plot supplementary figure 11I
media_P_sup_v_pel_replicates <- ggplot(combined_copies %>% 
                                      filter(!is.na(Cq), Sample %in% c("media+P")),
                                    aes(x = interaction(Sample, Timepoint, bio_rep), y = Cq, color = tech_replicate)) +
  geom_point(size = 3, alpha = 0.5) +
  labs(x = "Sample.timepoint.biorep", y = "Cq") +
  scale_color_manual(values = tech_rep_colors) +  # Apply color scale
  theme_bw() +
  theme(text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        legend.background = element_rect(linetype = 1, 
                                         linewidth = 0.5, 
                                         color = "gray"), 
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  coord_cartesian(ylim = c(0, 35))

media_P_sup_v_pel_replicates



#generate plot for figure 7 panel F
bd_sup_pel <- ggplot(combined_copies %>% filter(Sample %in% c("Bd+P_sup_v_pel", "Bd+P_sup", "Bd+P_pel", "media+P")), 
                     aes(x = Timepoint)) +
  geom_point(aes(y = Avg_copies_ml, color = Sample), shape = 16, size = 2) +
  geom_line(aes(y = Avg_copies_ml, color = Sample, group = Sample)) +
  scale_color_manual(values = c("Bd+P_sup_v_pel" = "darkgreen", 
                                "Bd+P_sup" = "black", 
                                "Bd+P_pel" = "#009292FF", 
                                "media+P" = "grey"),
                     labels = c("Bd+P_sup_v_pel" = "P. dorei + phage - whole cell suspension", 
                                "Bd+P_sup" = "P. dorei + phage - supernatant", 
                                "Bd+P_pel" = "P. dorei + phage - pellet", 
                                "media+P" = "Media + phage")) +
  labs(title = "", x = "Time (Hours)", y = "copies per mL") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none",  # Adjust legend position (x, y)
        legend.background = element_rect(fill = "white", color = "gray"), # Optional: Background color and border
        legend.title = element_blank(),
        text = element_text(size = 10),
        axis.line = element_line(linetype = 1, linewidth = 1),
        panel.border = element_rect(linetype = 1, linewidth = 1),
        panel.grid.minor = element_blank()) +
  geom_errorbar(aes(ymin = Avg_copies_ml - stdev_avg_copies_ml, 
                    ymax = Avg_copies_ml + stdev_avg_copies_ml, 
                    color = Sample), width = 0.5) +
  scale_y_continuous(trans = 'log10', breaks = 10^(3:10), limits = c(10^3, 10^10),
                     labels = scales::trans_format("log10", scales::math_format(10^.x)))


bd_sup_pel

#generate plot for figure 7 panel E
bv_sup_pel <- ggplot(combined_copies %>% filter(Sample %in% c("Bv+P_sup_v_pel", "Bv+P_sup", "Bv+P_pel", "media+P")), 
                     aes(x = Timepoint)) +
  geom_point(aes(y = Avg_copies_ml, color = Sample), shape = 16, size = 2) +
  geom_line(aes(y = Avg_copies_ml, color = Sample, group = Sample)) +
  scale_color_manual(values = c("Bv+P_sup_v_pel" = "darkgreen", 
                                "Bv+P_sup" = "black", 
                                "Bv+P_pel" = "#009292FF", 
                                "media+P" = "grey"),
                     labels = c("Bv+P_sup_v_pel" = "whole cell suspension", 
                                "Bv+P_sup" = "supernatant", 
                                "Bv+P_pel" = "pellet", 
                                "media+P" = "Media")) +
  labs(title = "", x = "Time (Hours)", y = "copies per mL") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "top",  # Adjust legend position (x, y)
        legend.background = element_rect(fill = "white", color = "gray"), # Optional: Background color and border
        legend.title = element_blank(),
        text = element_text(size = 10),
        axis.line = element_line(linetype = 1, linewidth = 1),
        panel.border = element_rect(linetype = 1, linewidth = 1),
        panel.grid.minor = element_blank()) +
  geom_errorbar(aes(ymin = Avg_copies_ml - stdev_avg_copies_ml, 
                    ymax = Avg_copies_ml + stdev_avg_copies_ml, 
                    color = Sample), width = 0.5) +
  scale_y_continuous(trans = 'log10', breaks = 10^(3:10), limits = c(10^3, 10^10),
                     labels = scales::trans_format("log10", scales::math_format(10^.x)))


bv_sup_pel


#######plot DNase #######

plate_1 <- read_csv("~/Documents/240614_T0_T40_T44_part_T24.csv", col_names = TRUE)
plate_2 <- read_csv("~/Documents/240618_sup_v_pel_qPCR_results.csv", col_names = TRUE)


####### Calculate Cq_mean, enter Std. curve, calculate copies for plate 1 
#only take the columns you are interested in, change Undetermined values to NA
trimmed_raw_results_P1 <- plate_1 %>% 
  select("Sample", "Cq") %>% 
  filter(Sample != 0) %>% 
  mutate(Cq = if_else(Cq == "Undetermined", NA_character_, Cq)) %>% 
  mutate(Cq = as.numeric(Cq))

#Calculate the mean Cq value for each sample
mean_values_P1 <- trimmed_raw_results_P1 %>% 
  group_by(Sample) %>% 
  summarize(Cq_mean = mean(Cq))

mean_values_P1 <- mean_values_P1 %>% filter(Sample != 0)

#merge the raw data with the calculated mean, denote tech reps, 
#filter for samples with tech reps within 2 Cq difference
results_with_mean_P1 <- trimmed_raw_results_P1 %>% 
  left_join(mean_values_P1, by = "Sample") %>% 
  group_by(Sample) %>% 
  mutate(tech_replicate=row_number(Sample)) %>% 
  mutate(Cq_diff = case_when(abs(combn(Cq, 2, function(x) abs(diff(x)))) > 3 ~ "Y", .default = "N")) %>% 
  filter(Cq_diff == "N") %>% 
  filter(!any(is.na(Cq)))


#Filter for the standard curve 
std_curve_plate_1 <- results_with_mean_P1 %>% 
  filter(str_detect(Sample, "std"))

#Use NEBiocalculator to look at standard curve and add the slope and intercept here: 
NEB_intercept = 34.34 #includes -2 to -9 
NEB_slope = -3.4 #includes -2 to -9
Cq_cutoff <- 31.236182
Cq_cutoff

#Use the standard curve coefficients to calculate the copies per ul for samples
plate1_copies <- results_with_mean_P1 %>%
  filter(Cq_mean < Cq_cutoff) %>%
  mutate(copies_ul = 10^((Cq_mean - NEB_intercept)/NEB_slope)) %>%
  mutate(copies_ml = copies_ul * 1000) %>%
  separate(Sample, into = c("Sample", "Timepoint"), sep = "_T") %>%
  separate(Sample, into = c("Sample", "bio_rep"), sep = "_(?=\\d+$)", extra = "merge", convert = TRUE) %>%
  filter(!(Sample == "Bd_P__nanopore" & bio_rep == 4)) %>%  # Remove the specific sample before calculation
  filter(!(Sample == "Bd_P__nanopore" & bio_rep == 5)) %>%
  group_by(Sample, Timepoint) %>%
  mutate(Avg_copies_ml = mean(copies_ml, na.rm = TRUE)) %>%
  mutate(stdev_avg_copies_ml = sd(copies_ml, na.rm = TRUE))


####### Calculate Cq_mean, enter Std. curve, calculate copies for plate 1 
#only take the columns you are interested in, change Undetermined values to NA
trimmed_raw_results_P2 <- plate_2 %>% 
  select("Sample", "Cq") %>% 
  filter(Sample != 0) %>% 
  mutate(Cq = if_else(Cq == "Undetermined", NA_character_, Cq)) %>% 
  mutate(Cq = as.numeric(Cq))

#Calculate the mean Cq value for each sample
mean_values_P2 <- trimmed_raw_results_P2 %>% 
  group_by(Sample) %>% 
  summarize(Cq_mean = mean(Cq))

mean_values_P2 <- mean_values_P2 %>% filter(Sample != 0)

#merge the raw data with the calculated mean, denote tech reps, 
#filter for samples with tech reps within 2 Cq difference
results_with_mean_P2 <- trimmed_raw_results_P2 %>% 
  left_join(mean_values_P2, by = "Sample") %>% 
  group_by(Sample) %>% 
  mutate(tech_replicate=row_number(Sample)) %>% 
  mutate(Cq_diff = case_when(abs(combn(Cq, 2, function(x) abs(diff(x)))) > 3 ~ "Y", .default = "N")) %>% 
  filter(Cq_diff == "N") %>% 
  filter(!any(is.na(Cq)))


#Filter for the standard curve 
std_curve_plate_2 <- results_with_mean_P2 %>% 
  filter(str_detect(Sample, "std"))

#Use NEBiocalculator to look at standard curve and add the slope and intercept here: 
NEB_intercept = 36.07 #includes -2 to -9 
NEB_slope = -3.5 #includes -2 to -9
Cq_cutoff <- 32.747552
Cq_cutoff

#Use the standard curve coefficients to calculate the copies per ul for samples
plate2_copies <- results_with_mean_P2 %>%
  filter(Cq_mean < Cq_cutoff) %>%
  mutate(copies_ul = 10^((Cq_mean - NEB_intercept)/NEB_slope)) %>%
  mutate(copies_ml = copies_ul * 1000) %>%
  separate(Sample, into = c("Sample", "Timepoint"), sep = "_T") %>%
  separate(Sample, into = c("Sample", "bio_rep"), sep = "_(?=\\d+$)", extra = "merge", convert = TRUE) %>%
  filter(!(Sample == "Bd_P__nanopore" & bio_rep == 4)) %>%  # Remove the specific sample before calculation
  filter(!(Sample == "Bd_P__nanopore" & bio_rep == 5)) %>%
  group_by(Sample, Timepoint) %>%
  mutate(Avg_copies_ml = mean(copies_ml, na.rm = TRUE)) %>%
  mutate(stdev_avg_copies_ml = sd(copies_ml, na.rm = TRUE))



combined_copies_1 <- bind_rows(plate2_copies, plate1_copies)


combined_copies <- combined_copies_1 %>%
  mutate(Timepoint = as.numeric(Timepoint))

dnase_bd <- ggplot(combined_copies %>% filter(Sample %in% c("Bd_P__nanopore", "Bd_P_dnase", "media_P_dnase", "media_P")), 
                   aes(x = Timepoint)) +
  geom_point(aes(y = Avg_copies_ml, color = Sample), shape = 16, size = 2) +
  geom_line(aes(y = Avg_copies_ml, color = Sample, group = Sample)) +
  scale_color_manual(values = c("Bd_P__nanopore" = "#009292FF", "Bd_P_dnase" = "black", "media_P_dnase" = "lightgrey", "media_P" = "darkgrey")) +
  labs(title = "", x = "Time (Hours)", y = "copies per mL") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "right",
        text = element_text(size = 10),
        panel.grid.major = element_line(linewidth = 1),
        panel.grid.minor = element_line(linewidth = 1),
        axis.line = element_line(linetype = 1, linewidth = 1),
        panel.border = element_rect(linetype = 1, linewidth = 1),
        legend.background = element_rect(linetype = 1, linewidth = 0.5, color = "gray")) +
  geom_errorbar(aes(ymin = Avg_copies_ml - stdev_avg_copies_ml, 
                    ymax = Avg_copies_ml + stdev_avg_copies_ml, 
                    color = Sample), width = 1) +
  scale_y_continuous(trans = 'log10', breaks = 10^(3:8), limits = c(10^3, 10^8))


dnase_bd

#######plot single cell ratio########

# Read the Excel file
data <- read_excel("~/Documents/single_cell_seq_results.xlsx", sheet = "Sheet10")

# Sort the data frame based on the `pcrAss_bv_ratio` column
data <- data[order(data$pcrAss_bv_ratio), ]

# Remove duplicates from sample names and reorder the factor
unique_samples <- unique(data$sample)
data$sample <- factor(data$sample, levels = unique_samples[order(match(unique_samples, data$sample))])

p <- ggplot(data, aes(x = sample, y = pcrAss_bv_ratio)) +
  geom_point(size=3) +
  geom_hline(yintercept = 10, linetype = "dashed", color = "grey", size = 1) +
  labs(title = "",
       x = "sample",
       y = "C. communis:B. vulgatus coverage ratio") +
  theme_bw() +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    panel.grid.major = element_line(linewidth = 1),
    panel.grid.minor = element_line(linewidth = 1),
    axis.line = element_line(linetype = 1, linewidth = 1),
    panel.border = element_rect(linetype = 1, linewidth = 2),
    legend.background = element_rect(linetype = 1, linewidth = 0.5, color = "gray")
  ) +
  scale_y_log10(breaks = c(0.1, 1, 10, 100)) +
  annotate("text", x = 1, y = 1, label = "1", vjust = -0.5, hjust = 1, size = 5)
p


pdf(here("~/Documents/single_cell_ratio.pdf"), width = 10, height = 10)
print(p)
dev.off()

#######plots OD600 supplementary figures ######

indata <- read_excel("~/Documents/231118_od_combined.xlsx")


bv <- ggplot(indata %>% filter(sample == "Bv_phage_1" | sample == "Bv_1" | sample == "media_phage_1" | sample == "media_1"), aes(x = time, y = avg_OD)) + 
  geom_point(aes(color = sample, shape = sample), size = 2) + 
  labs(title = "P. vulgatus", x = "Time (Hours)", y = "OD600") +
  geom_line(aes(color = sample)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "italic"), legend.position = "right") +
  scale_color_manual(values = c("Bv_phage_1" = "#009292FF", "Bv_1" = "purple", "media_phage_1" = "darkgrey", "media_1" = "black"),
                     name = "Sample", labels = c("Bv_phage_1" = "B. vulgatus + phage", "Bv_1" = "B. vulgatus", "media_phage_1" = "media + phage", "media_1" = "media")) +
  scale_shape_manual(values = c("Bv_phage_1" = 0, "Bv_1" = 1, "media_phage_1" = 2, "media_1" = 3),
                     name = "Sample", labels = c("Bv_phage_1" = "B. vulgatus + phage", "Bv_1" = "B. vulgatus", "media_phage_1" = "media + phage", "media_1" = "media")) +
  geom_errorbar(aes(x = time, ymin = avg_OD - stdev, ymax = avg_OD + stdev, color = sample), width = 0.2) +
  scale_y_log10()  # Adding log scale for the y-axis
bv
ggsave(here("~/Documents/240308_OD600_Bv.pdf"), dpi=300, w=10, h=9)




bd <- ggplot(indata %>% filter(sample == "Bd_phage_1" | sample == "Bd_1" | sample == "media_phage_1" | sample == "media_1"), aes(x = time, y = avg_OD)) + 
  geom_point(aes(color = sample, shape = sample), size = 2) + 
  labs(title = "P. dorei", x = "Time (Hours)", y = "OD600") +
  geom_line(aes(color = sample)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "italic"), legend.position = "right") +
  scale_color_manual(values = c("Bd_phage_1" = "#009292FF", "Bd_1" = "purple", "media_phage_1" = "darkgrey", "media_1" = "black"),
                     name = "Sample", labels = c("Bd_phage_1" = "B. dorei + phage", "Bd_1" = "B. dorei", "media_phage_1" = "media + phage", "media_1" = "media")) +
  scale_shape_manual(values = c("Bd_phage_1" = 0, "Bd_1" = 1, "media_phage_1" = 2, "media_1" = 3),
                     name = "Sample", labels = c("Bd_phage_1" = "B. dorei + phage", "Bd_1" = "B. dorei", "media_phage_1" = "media + phage", "media_1" = "media")) +
  geom_errorbar(aes(x = time, ymin = avg_OD - stdev, ymax = avg_OD + stdev, color = sample), width = 0.2) +
  scale_y_log10()  # Adding log scale for the y-axis
bd
ggsave(here("~/Documents/240308_OD600_Bd.pdf"), dpi=300, w=10, h=9)



Pm <- ggplot(indata %>% filter(sample == "Pm_phage_1" | sample == "Pm_1" | sample == "media_phage_1" | sample == "media_1"), aes(x = time, y = avg_OD)) + 
  geom_point(aes(color = sample, shape = sample), size = 2) + 
  labs(title = "P. merdae", x = "Time (Hours)", y = "OD600") +
  geom_line(aes(color = sample)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "italic"), legend.position = "right") +
  scale_color_manual(values = c("Pm_phage_1" = "#009292FF", "Pm_1" = "purple", "media_phage_1" = "darkgrey", "media_1" = "black"),
                     name = "Sample", labels = c("Pm_phage_1" = "P. merdae + phage", "Pm_1" = "P. merdae", "media_phage_1" = "media + phage", "media_1" = "media")) +
  scale_shape_manual(values = c("Pm_phage_1" = 0, "Pm_1" = 1, "media_phage_1" = 2, "media_1" = 3),
                     name = "Sample", labels = c("Pm_phage_1" = "P. merdae + phage", "Pm_1" = "P. merdae", "media_phage_1" = "media + phage", "media_1" = "media")) +
  geom_errorbar(aes(x = time, ymin = avg_OD - stdev, ymax = avg_OD + stdev, color = sample), width = 0.2) +
  scale_y_log10()  # Adding log scale for the y-axis
Pm
ggsave(here("~/Documents/240308_OD600_Pm.pdf"), dpi=300, w=10, h=9)


Bs <- ggplot(indata %>% filter(sample == "Bs_phage_1" | sample == "Bs_1" | sample == "media_phage_1" | sample == "media_1"), aes(x = time, y = avg_OD)) + 
  geom_point(aes(color = sample, shape = sample), size = 2) + 
  labs(title = "B. stercoris", x = "Time (Hours)", y = "OD600") +
  geom_line(aes(color = sample)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "italic"), legend.position = "right") +
  scale_color_manual(values = c("Bs_phage_1" = "#009292FF", "Bs_1" = "purple", "media_phage_1" = "darkgrey", "media_1" = "black"),
                     name = "Sample", labels = c("Bs_phage_1" = "B. stercoris + phage", "Bs_1" = "B. stercoris", "media_phage_1" = "media + phage", "media_1" = "media")) +
  scale_shape_manual(values = c("Bs_phage_1" = 0, "Bs_1" = 1, "media_phage_1" = 2, "media_1" = 3),
                     name = "Sample", labels = c("Bs_phage_1" = "B. stercoris + phage", "Bs_1" = "B. stercoris", "media_phage_1" = "media + phage", "media_1" = "media")) +
  geom_errorbar(aes(x = time, ymin = avg_OD - stdev, ymax = avg_OD + stdev, color = sample), width = 0.2) +
  scale_y_log10()  # Adding log scale for the y-axis
Bs
ggsave(here("~/Documents/240308_OD600_Bs.pdf"), dpi=300, w=10, h=9)

OD600_data <- read_excel("~/Documents/231001_OD600_bulk_ferm_combined.xlsx")
OD600_data$sample <- factor(OD600_data$sample, levels = c("1.1", 
                                                          "2.1", 
                                                          "3.1", 
                                                          "4.1", 
                                                          "5.1", 
                                                          "6.1"))

OD600_data <- OD600_data %>% filter(time != 48)
OD600_plot_abx_p4100 <- ggplot(OD600_data %>% filter(sample %in% c("4.1", 
                                                                   "6.1")), 
                               aes(x = time, y = `average OD600`)) + 
  geom_point(aes(color = sample), shape = 16, size = 4) +
  labs(title = expression(paste("")), 
       x = "Time (Hours)", y = "OD600") +
  geom_line(aes(color = sample), size = 0.75) +
  scale_color_manual(values = c("4.1" = "#009292FF", 
                                "6.1" = "#276A91"), 
                     name = "Sample", labels = c("p-crAss positive stool", "p-crAss negative stool")) +
  theme_bw() +
  scale_y_log10() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        text = element_text(size = 40),
        panel.grid.major = element_line(linewidth = 1.5),
        panel.grid.minor = element_line(linewidth = 1.5),
        axis.line = element_line(linetype = 1, linewidth = 1),
        panel.border = element_rect(linetype = 1, linewidth = 2.5),
        legend.background = element_rect(linetype = 1, linewidth = 0.5, color = "gray")) +
  geom_errorbar(aes(x = time, ymin = `average OD600` - `stdev OD600`, 
                    ymax = `average OD600` + `stdev OD600`, 
                    color = sample), 
                size = 0.75, width = 2)

OD600_plot_abx_p4100

ggsave(here("~/Documents/240308_OD600_bulk_ferm_p4100.pdf"), dpi=300, w=10, h=9)


#######plot combined qPCR standards for supplementary figure (individual standards are plotted in each plate analysis section)#####

#Create dataframe with linear regressions
x <- seq(0, 15, by = 1)
data <- data.frame(x = x, 
                   iplate1 = 35.78 - 3.53 * x, 
                   iplate2 = 37.41 - 3.6 * x, 
                   iplate3 = 36.69 - 3.6 * x,
                   iplate4 = 36.78 - 3.74 * x,
                   bfplate1 = 36.7 - 3.54 * x,
                   bfplate2 = 36.03 - 3.42 * x,
                   bfplate3 = 34.47 - 3.43 * x,
                   suppelplate1 = 37.72 - 3.54 * x,
                   suppelplate2 = 40.65 - 3.97 * x,
                   suppelplate3 = 36.97 - 3.48 * x,
                   suppelplate4 = 41.18 - 3.93 * x,
                   average = 37.31 - 3.62 * x,
                   average_plus_3 = 40.31 - 3.62 * x,
                   average_minus_3 = 34.31 - 3.62 * x)


# Create plot with all p-crass qPCR standard curves
ggplot(data, aes(x = x)) +
  geom_line(aes(y = iplate1, color = "isolate Plate 1 standard")) +
  geom_line(aes(y = iplate2, color = "isolate Plate 2 standard")) +
  geom_line(aes(y = iplate3, color = "isolate Plate 3 standard")) +
  geom_line(aes(y = iplate4, color = "isolate Plate 4 standard")) +
  geom_line(aes(y = bfplate1, color = "bulk ferm Plate 1 standard")) +
  geom_line(aes(y = bfplate2, color = "bulk ferm Plate 2 standard")) +
  geom_line(aes(y = bfplate3, color = "bulk ferm Plate 3 standard")) +
  geom_line(aes(y = suppelplate1, color = "sup v pel Plate 1 standard")) +
  geom_line(aes(y = suppelplate2, color = "sup v pel Plate 2 standard")) +
  geom_line(aes(y = suppelplate3, color = "sup v pel Plate 3 standard")) +  
  geom_line(aes(y = suppelplate4, color = "sup v pel Plate 4 standard")) +
  geom_line(aes(y = average, color = "average"), linetype = "dotted") +
  geom_line(aes(y = average_plus_3, color = "average_plus_3"), linetype = "dotted") +
  geom_line(aes(y = average_minus_3, color = "average_minus_3"), linetype = "dotted") +
  scale_color_manual(values = c("isolate Plate 1 standard" = "#009292FF",
                                "isolate Plate 2 standard" = "#35B779FF",
                                "isolate Plate 3 standard" = "darkgreen",
                                "isolate Plate 4 standard" = "#404688FF",
                                "bulk ferm Plate 1 standard" = "#5DC863FF",
                                "bulk ferm Plate 2 standard" = "#440154FF",
                                "bulk ferm Plate 3 standard" = "#018060",
                                "sup v pel Plate 1 standard" = "green",
                                "sup v pel Plate 2 standard" = "lightgreen",
                                "sup v pel Plate 3 standard" = "yellowgreen",
                                "sup v pel Plate 4 standard" = "green3",
                                "average" = "black",
                                "average_plus_3" = "black",
                                "average_minus_3" = "black")) +
  labs(x = "log Copy Number", y = "Cq", color="") +
  theme_bw()




labs(x = "logCopyNumber", y = "Cq", color="") +
  xlim(0,15) +
  ylim(-5, 30) +
  theme_bw() 


#Create dataframe with linear regressions
x <- seq(0, 15, by = 1)
data <- data.frame(x = x, 
                   bfplate1_bv = 32.21 - 3.54 * x, 
                   bfplate2_bv = 31.02 - 3.39 * x,
                   average = 31.62 - 3.47 * x,
                   average_plus_3 = 34.62 - 3.47 * x,
                   average_minus_3 = 28.62 - 3.47 * x)


# Create plot with all P. vulgatus qPCR standard curves
ggplot(data, aes(x = x)) +
  geom_line(aes(y = bfplate1_bv, color = "bulk ferm Plate 1")) +
  geom_line(aes(y = bfplate2_bv, color = "bulk ferm Plate 2")) +
  geom_line(aes(y = average, color = "average"), linetype = "dotted") +
  geom_line(aes(y = average_plus_3, color = "average_plus_3"), linetype = "dotted") +
  geom_line(aes(y = average_minus_3, color = "average_minus_3"), linetype = "dotted") +
  scale_color_manual(values = c("bulk ferm Plate 1" = "#009292FF",
                                "bulk ferm Plate 2" = "#35B779FF",
                                "average" = "black",
                                "average_plus_3" = "black",
                                "average_minus_3" = "black")) +
  labs(x = "log Copy Number", y = "Cq", color="") +
  theme_bw()




#open files with raw qPCR data
raw_results_plate_1 <- read_csv("~/Documents/240306_final_bulk_ferm_export_plate_1.csv", col_names = TRUE)
raw_results_plate_2 <- read_csv("~/Documents/240306_final_bulk_ferm_export_plate_2.csv", col_names = TRUE)
raw_results_vulgatus_plate_1 <- read_csv("~/Documents/240306_vulgatus_final_export_bulk_ferm_plate_1.csv", col_names = TRUE)
raw_results_vulgatus_plate_2 <- read_csv("~/Documents/240306_vulgatus_final_export_bulk_ferm_plate_2.csv", col_names = TRUE)
bulk_ferm_rerun <- read_csv("~/Documents/240304_redo_bulk_ferm.csv", col_names = TRUE)


####### Calculate Cq_mean, enter Std. curve, calculate copies for plate 1 communis copies #####
#only take the columns you are interested in, change Undetermined values to NA
trimmed_raw_results <- raw_results_plate_1 %>% 
  select("Sample", "Cq") %>% 
  mutate(Cq = if_else(Cq == "Undetermined", NA_character_, Cq)) %>% 
  mutate(Cq = as.numeric(Cq))

# Assuming "Sample" column is character type, you can use str_split to separate the parts
trimmed_raw_results <- trimmed_raw_results %>%
  mutate(Sample = str_split(Sample, " ") %>%
           map_chr(~ paste(rev(str_replace_all(.x, "\\.", "_")), collapse = "_")))

#Calculate the mean Cq value for each sample
mean_values <- trimmed_raw_results %>% 
  group_by(Sample) %>% 
  summarize(Cq_mean = mean(Cq))

#merge the raw data with the calculated mean, denote tech reps, 
#filter for samples with tech reps within 3 Cq difference


results_with_mean <- trimmed_raw_results %>%
  filter(!is.na(Cq)) %>%
  left_join(mean_values, by = "Sample") %>% 
  group_by(Sample) %>% 
  mutate(tech_replicate = row_number(Sample),
         Cq_diff = if (n() >= 2) {
           any(combn(Cq, 2, function(x) abs(diff(x))) > 3)
         } else {
           FALSE
         }) %>%
  filter(!Cq_diff)




#Filter for the standard curve 
std_curve_plate_1 <- results_with_mean %>% 
  filter(str_detect(Sample, "std"))

#Use NEBiocalculator to look at standard curve and add the slope and intercept here: 
NEB_intercept = 36.7 #includes -1 to -9 
NEB_slope = -3.54 #includes -1 to -9
Cq_cutoff <- std_curve_plate_1$Cq_mean[std_curve_plate_1$Sample == "-9_std"][1]
Cq_cutoff

# if you want to calculate std. curve in R #
# #Calculate the standard curve in R 
std_curve_plate_1 <- results_with_mean %>% 
  filter(str_detect(Sample, "std")) %>%
  group_by(Sample) %>%
  mutate("copies_ul" = case_when(
    Sample == "-1_std" ~ 1000000000,
    Sample == "-2_std" ~ 100000000,
    Sample == "-3_std" ~ 10000000,
    Sample == "-4_std" ~ 1000000,
    Sample == "-5_std" ~ 100000,
    Sample == "-6_std" ~ 10000,
    Sample == "-7_std" ~ 1000,
    Sample == "-8_std" ~ 100,
    Sample == "-9_std" ~ 10,
  )) %>% 
  filter(!is.na(Cq_mean)) 
std_curve_plate_1 %>%
  ggplot(aes(x = copies_ul, y = Cq_mean )) +
  geom_point() +
  theme_bw() + 
  scale_x_continuous(trans = "log10") +
  geom_smooth(method = "lm", se = FALSE) 

# Fit a linear regression model
lm_model_plate_1 <- lm(Cq_mean ~ log10(copies_ul), data = std_curve_plate_1)

# Extract coefficients and R-squared value
lm_summary_plate_1 <- summary(lm_model_plate_1)
intercept_plate_1 <- coef(lm_model_plate_1)[1]
slope_plate_1 <- coef(lm_model_plate_1)[2]
r_squared_plate_1 <- lm_summary_plate_1$r.squared

# Print intercept, slope, and R-squared value
cat("Intercept:", intercept_plate_1, "\n")
cat("Slope:", slope_plate_1, "\n")
cat("R-squared:", r_squared_plate_1, "\n")

# Plot with linear regression line
std_curve_plate_1 %>%
  ggplot(aes(x = copies_ul, y = Cq_mean)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  theme_bw() + 
  scale_x_continuous(trans = "log10") +
  labs(x = "Log10 Copies/ul", y = "Cq mean") +
  geom_text(x = 1e3, y = 35, label = paste("y =", round(intercept_plate_1, 2), "+", round(slope_plate_1, 2), "* log10(x)", "\nR-squared =", round(r_squared_plate_1, 4)), color = "blue")



#Use the standard curve coefficients to calculate the copies per ul for samples
plate1_copies <- results_with_mean %>%
  filter(Cq_mean < Cq_cutoff) %>% 
  mutate(copies_ul = 10^((Cq_mean - NEB_intercept)/NEB_slope)) %>% #Change here if you are calculating std. curve from R
  mutate(copies_ml = copies_ul * 1000) %>% 
  separate(Sample, into = c("Sample", "Timepoint"), sep = "_T") %>% 
  separate(Sample, into = c("Sample", "bio_rep"), sep = "_(?=\\d+$)", extra = "merge", convert = TRUE) %>% 
  group_by(Sample, Timepoint) %>% 
  mutate(Avg_copies_ml = mean(copies_ml, na.rm = TRUE)) %>% 
  mutate(stdev_avg_copies_ml = sd(copies_ml, na.rm = TRUE))



####### Calculate Cq_mean, enter Std. curve, calculate copies for plate 2 communis copies #####
#only take the columns you are interested in, change Undetermined values to NA
trimmed_raw_results <- raw_results_plate_2 %>% 
  select("Sample", "Cq") %>% 
  mutate(Cq = if_else(Cq == "Undetermined", NA_character_, Cq)) %>% 
  mutate(Cq = as.numeric(Cq))

# Assuming "Sample" column is character type, you can use str_split to separate the parts
trimmed_raw_results <- trimmed_raw_results %>%
  mutate(Sample = str_split(Sample, " ") %>%
           map_chr(~ paste(rev(str_replace_all(.x, "\\.", "_")), collapse = "_")))

#Calculate the mean Cq value for each sample
mean_values <- trimmed_raw_results %>% 
  group_by(Sample) %>% 
  summarize(Cq_mean = mean(Cq))

#merge the raw data with the calculated mean, denote tech reps, 
#filter for samples with tech reps within 2 Cq difference


results_with_mean <- trimmed_raw_results %>%
  filter(!is.na(Cq)) %>%
  left_join(mean_values, by = "Sample") %>% 
  group_by(Sample) %>% 
  mutate(tech_replicate = row_number(Sample),
         Cq_diff = if (n() >= 2) {
           any(combn(Cq, 2, function(x) abs(diff(x))) > 3)
         } else {
           FALSE
         }) %>%
  filter(!Cq_diff)




#Filter for the standard curve 
std_curve_plate_2 <- results_with_mean %>% 
  filter(str_detect(Sample, "std"))

#Use NEBiocalculator to look at standard curve and add the slope and intercept here: 
NEB_intercept = 36.03 #includes -6 to -9 
NEB_slope = -3.42 #includes -6 to -9
Cq_cutoff <- std_curve_plate_2$Cq_mean[std_curve_plate_2$Sample == "-8_std"][1]
Cq_cutoff

# if you want to calculate std. curve in R #
# #Calculate the standard curve in R 
std_curve_plate_2 <- results_with_mean %>% 
  filter(str_detect(Sample, "std")) %>%
  group_by(Sample) %>%
  mutate("copies_ul" = case_when(
    Sample == "-5_std" ~ 100000,
    Sample == "-6_std" ~ 10000,
    Sample == "-7_std" ~ 1000,
    Sample == "-8_std" ~ 100,
  )) %>% 
  filter(!is.na(Cq_mean)) 
std_curve_plate_2 %>%
  ggplot(aes(x = copies_ul, y = Cq_mean )) +
  geom_point() +
  theme_bw() + 
  scale_x_continuous(trans = "log10") +
  geom_smooth(method = "lm", se = FALSE) 

# Fit a linear regression model
lm_model_plate_2 <- lm(Cq_mean ~ log10(copies_ul), data = std_curve_plate_2)

# Extract coefficients and R-squared value
lm_summary_plate_2 <- summary(lm_model_plate_2)
intercept_plate_2 <- coef(lm_model_plate_2)[1]
slope_plate_2 <- coef(lm_model_plate_2)[2]
r_squared_plate_2 <- lm_summary_plate_2$r.squared

# Print intercept, slope, and R-squared value
cat("Intercept:", intercept_plate_2, "\n")
cat("Slope:", slope_plate_2, "\n")
cat("R-squared:", r_squared_plate_2, "\n")

# Plot with linear regression line
std_curve_plate_2 %>%
  ggplot(aes(x = copies_ul, y = Cq_mean)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  theme_bw() + 
  scale_x_continuous(trans = "log10") +
  labs(x = "Log10 Copies/ul", y = "Cq mean") +
  geom_text(x = 1e3, y = 35, label = paste("y =", round(intercept_plate_2, 2), "+", round(slope_plate_2, 2), "* log10(x)", "\nR-squared =", round(r_squared_plate_2, 4)), color = "blue")




#Use the standard curve coefficients to calculate the copies per ul for samples
plate2_copies <- results_with_mean %>%
  filter(Cq_mean < Cq_cutoff) %>% 
  mutate(copies_ul = 10^((Cq_mean - NEB_intercept)/NEB_slope)) %>% #Change here if you are calculating std. curve from R
  mutate(copies_ml = copies_ul * 1000) %>% 
  separate(Sample, into = c("Sample", "Timepoint"), sep = "_T") %>% 
  separate(Sample, into = c("Sample", "bio_rep"), sep = "_(?=\\d+$)", extra = "merge", convert = TRUE) %>% 
  group_by(Sample, Timepoint) %>% 
  mutate(Avg_copies_ml = mean(copies_ml, na.rm = TRUE)) %>% 
  mutate(stdev_avg_copies_ml = sd(copies_ml, na.rm = TRUE))


####### Calculate Cq_mean, enter Std. curve, calculate copies for communis samples that were reun #####
#only take the columns you are interested in, change Undetermined values to NA
trimmed_raw_results <- bulk_ferm_rerun %>% 
  select("Sample", "Cq") %>% 
  mutate(Cq = if_else(Cq == "Undetermined", NA_character_, Cq)) %>% 
  mutate(Cq = as.numeric(Cq)) %>%
  filter(Sample != 0)

#Calculate the mean Cq value for each sample
mean_values <- trimmed_raw_results %>% 
  group_by(Sample) %>% 
  summarize(Cq_mean = mean(Cq)) 

#merge the raw data with the calculated mean, denote tech reps, 
#filter for samples with tech reps within 2 Cq difference
results_with_mean <- trimmed_raw_results %>% 
  left_join(mean_values, by = "Sample") %>% 
  group_by(Sample) %>% 
  mutate(tech_replicate=row_number(Sample)) %>% 
  mutate(Cq_diff = case_when(abs(combn(Cq, 2, function(x) abs(diff(x)))) > 3 ~ "Y", .default = "N")) %>% 
  filter(Cq_diff == "N") %>% 
  filter(!any(is.na(Cq)))


#Filter for the standard curve 
std_curve_rerun_plate <- results_with_mean %>% 
  filter(str_detect(Sample, "std"))

#Use NEBiocalculator to look at standard curve and add the slope and intercept here: 
NEB_intercept = 34.47 #includes -2 to -9 
NEB_slope = -3.43 #includes -2 to -9
Cq_cutoff <- std_curve_rerun_plate$Cq_mean[std_curve_rerun_plate$Sample == "std -9"][1]
Cq_cutoff

# if you want to calculate std. curve in R #
# #Calculate the standard curve in R 
std_curve_rerun_plate <- results_with_mean %>% 
  filter(str_detect(Sample, "std")) %>%
  group_by(Sample) %>%
  mutate("copies_ul" = case_when(
    Sample == "std -1" ~ 1000000000,
    Sample == "std -2" ~ 100000000,
    Sample == "std -3" ~ 10000000,
    Sample == "std -4" ~ 1000000,
    Sample == "std -5" ~ 100000,
    Sample == "std -6" ~ 10000,
    Sample == "std -7" ~ 1000,
    Sample == "std -8" ~ 100,
    Sample == "std -9" ~ 10,
  )) %>% 
  filter(!is.na(Cq_mean)) 
std_curve_rerun_plate %>%
  ggplot(aes(x = copies_ul, y = Cq_mean )) +
  geom_point() +
  theme_bw() + 
  scale_x_continuous(trans = "log10") +
  geom_smooth(method = "lm", se = FALSE) 

# Fit a linear regression model
lm_model_plate_re <- lm(Cq_mean ~ log10(copies_ul), data = std_curve_rerun_plate)

# Extract coefficients and R-squared value
lm_summary_plate_re <- summary(lm_model_plate_re)
intercept_plate_re <- coef(lm_model_plate_re)[1]
slope_plate_re <- coef(lm_model_plate_re)[2]
r_squared_plate_re <- lm_summary_plate_re$r.squared

# Print intercept, slope, and R-squared value
cat("Intercept:", intercept_plate_re, "\n")
cat("Slope:", slope_plate_re, "\n")
cat("R-squared:", r_squared_plate_re, "\n")

# Plot with linear regression line
std_curve_rerun_plate %>%
  ggplot(aes(x = copies_ul, y = Cq_mean)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  theme_bw() + 
  scale_x_continuous(trans = "log10") +
  labs(x = "Log10 Copies/ul", y = "Cq mean") +
  geom_text(x = 1e3, y = 35, label = paste("y =", round(intercept_plate_re, 2), "+", round(slope_plate_re, 2), "* log10(x)", "\nR-squared =", round(r_squared_plate_re, 4)), color = "blue")





#Use the standard curve coefficients to calculate the copies per ul for samples
rerun_plate_copies <- results_with_mean %>%
  filter(Cq_mean < Cq_cutoff) %>% 
  mutate(copies_ul = 10^((Cq_mean - NEB_intercept)/NEB_slope)) %>% #Change here if you are calculating std. curve from R
  mutate(copies_ml = copies_ul * 1000) %>% 
  separate(Sample, into = c("Timepoint", "Sample", "bio_rep"), sep = "_T|_", extra = "merge", convert = TRUE) %>%
  mutate(Timepoint = str_remove(Timepoint, "T")) %>% # Remove "T" from the Timepoint column
  group_by(Sample, Timepoint) %>% 
  mutate(Avg_copies_ml = mean(copies_ml, na.rm = TRUE)) %>% 
  mutate(stdev_avg_copies_ml = sd(copies_ml, na.rm = TRUE))



####### Calculate Cq_mean, enter Std. curve, calculate copies for plate 1 vulgatus copies #####
#only take the columns you are interested in, change Undetermined values to NA
trimmed_raw_results <- raw_results_vulgatus_plate_1 %>% 
  select("Sample", "Cq") %>% 
  mutate(Cq = if_else(Cq == "Undetermined", NA_character_, Cq)) %>% 
  mutate(Cq = as.numeric(Cq))

# Assuming "Sample" column is character type, you can use str_split to separate the parts
trimmed_raw_results <- trimmed_raw_results %>%
  mutate(Sample = str_split(Sample, " ") %>%
           map_chr(~ paste(rev(str_replace_all(.x, "\\.", "_")), collapse = "_")))

#Calculate the mean Cq value for each sample
mean_values <- trimmed_raw_results %>% 
  group_by(Sample) %>% 
  summarize(Cq_mean = mean(Cq))

#merge the raw data with the calculated mean, denote tech reps, 
#filter for samples with tech reps within 2 Cq difference


results_with_mean <- trimmed_raw_results %>%
  filter(!is.na(Cq)) %>%
  left_join(mean_values, by = "Sample") %>% 
  group_by(Sample) %>% 
  mutate(tech_replicate = row_number(Sample),
         Cq_diff = if (n() >= 2) {
           any(combn(Cq, 2, function(x) abs(diff(x))) > 3)
         } else {
           FALSE
         }) %>%
  filter(!Cq_diff)




#Filter for the standard curve 
std_curve_vulgatus_plate_1 <- results_with_mean %>% 
  filter(str_detect(Sample, "std"))

#Use NEBiocalculator to look at standard curve and add the slope and intercept here: 
NEB_intercept = 32.21 #includes -3 to -6 
NEB_slope = -3.54 #includes -3 to -6
Cq_cutoff <- std_curve_vulgatus_plate_1$Cq_mean[std_curve_vulgatus_plate_1$Sample == "-8_std"][1]
Cq_cutoff

# if you want to calculate std. curve in R #
# #Calculate the standard curve in R 
std_curve_vulgatus_plate_1 <- results_with_mean %>% 
  filter(str_detect(Sample, "std")) %>%
  group_by(Sample) %>%
  mutate("copies_ul" = case_when(
    Sample == "-3_std" ~ 10000000,
    Sample == "-4_std" ~ 1000000,
    Sample == "-5_std" ~ 100000,
    Sample == "-6_std" ~ 10000,
  )) %>% 
  filter(!is.na(Cq_mean)) 
std_curve_vulgatus_plate_1 %>%
  ggplot(aes(x = copies_ul, y = Cq_mean )) +
  geom_point() +
  theme_bw() + 
  scale_x_continuous(trans = "log10") +
  geom_smooth(method = "lm", se = FALSE) 

# Fit a linear regression model
lm_model_plate_bv_1 <- lm(Cq_mean ~ log10(copies_ul), data = std_curve_vulgatus_plate_1)

# Extract coefficients and R-squared value
lm_summary_plate_bv_1 <- summary(lm_model_plate_bv_1)
intercept_plate_bv_1 <- coef(lm_model_plate_bv_1)[1]
slope_plate_bv_1 <- coef(lm_model_plate_bv_1)[2]
r_squared_plate_bv_1 <- lm_summary_plate_bv_1$r.squared

# Print intercept, slope, and R-squared value
cat("Intercept:", intercept_plate_bv_1, "\n")
cat("Slope:", slope_plate_bv_1, "\n")
cat("R-squared:", r_squared_plate_bv_1, "\n")

# Plot with linear regression line
std_curve_vulgatus_plate_1 %>%
  ggplot(aes(x = copies_ul, y = Cq_mean)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  theme_bw() + 
  scale_x_continuous(trans = "log10") +
  labs(x = "Log10 Copies/ul", y = "Cq mean") +
  geom_text(x = 1e3, y = 35, label = paste("y =", round(intercept_plate_bv_1, 2), "+", round(slope_plate_bv_1, 2), "* log10(x)", "\nR-squared =", round(r_squared_plate_bv_1, 4)), color = "blue")


#Use the standard curve coefficients to calculate the copies per ul for samples
plate1_vulgatus_copies <- results_with_mean %>%
  filter(Cq_mean < Cq_cutoff) %>% 
  mutate(copies_ul = 10^((Cq_mean - NEB_intercept)/NEB_slope)) %>% #Change here if you are calculating std. curve from R
  mutate(copies_ml = copies_ul * 1000) %>% 
  separate(Sample, into = c("Sample", "Timepoint"), sep = "_T") %>% 
  separate(Sample, into = c("Sample", "bio_rep"), sep = "_(?=\\d+$)", extra = "merge", convert = TRUE) %>% 
  group_by(Sample, Timepoint) %>% 
  mutate(Avg_copies_vulgatus_ml = mean(copies_ml, na.rm = TRUE)) %>% 
  mutate(stdev_avg_vulgatus_copies_ml = sd(copies_ml, na.rm = TRUE))
####### Calculate Cq_mean, enter Std. curve, calculate copies for plate 2 vulgatus copies #####
#only take the columns you are interested in, change Undetermined values to NA
trimmed_raw_results <- raw_results_vulgatus_plate_2 %>% 
  select("Sample", "Cq") %>% 
  mutate(Cq = if_else(Cq == "Undetermined", NA_character_, Cq)) %>% 
  mutate(Cq = as.numeric(Cq))

# Assuming "Sample" column is character type, you can use str_split to separate the parts
trimmed_raw_results <- trimmed_raw_results %>%
  mutate(Sample = str_split(Sample, " ") %>%
           map_chr(~ paste(rev(str_replace_all(.x, "\\.", "_")), collapse = "_")))

#Calculate the mean Cq value for each sample
mean_values <- trimmed_raw_results %>% 
  group_by(Sample) %>% 
  summarize(Cq_mean = mean(Cq))

#merge the raw data with the calculated mean, denote tech reps, 
#filter for samples with tech reps within 2 Cq difference


results_with_mean <- trimmed_raw_results %>%
  filter(!is.na(Cq)) %>%
  left_join(mean_values, by = "Sample") %>% 
  group_by(Sample) %>% 
  mutate(tech_replicate = row_number(Sample),
         Cq_diff = if (n() >= 2) {
           any(combn(Cq, 2, function(x) abs(diff(x))) > 3)
         } else {
           FALSE
         }) %>%
  filter(!Cq_diff)




#Filter for the standard curve 
std_curve_vulgatus_plate_2 <- results_with_mean %>% 
  filter(str_detect(Sample, "std"))

#Use NEBiocalculator to look at standard curve and add the slope and intercept here: 
NEB_intercept = 32.20 #includes -2 to -8 
NEB_slope = -3.59 #includes -2 to -8
Cq_cutoff = 30
Cq_cutoff

# if you want to calculate std. curve in R #
# #Calculate the standard curve in R 
std_curve_vulgatus_plate_2 <- results_with_mean %>% 
  filter(str_detect(Sample, "std")) %>%
  group_by(Sample) %>%
  mutate("copies_ul" = case_when(
    Sample == "-2_std" ~ 100000000,
    Sample == "-3_std" ~ 10000000,
    Sample == "-4_std" ~ 1000000,
    Sample == "-5_std" ~ 100000,
    Sample == "-6_std" ~ 10000,
    Sample == "-7_std" ~ 1000,
    Sample == "-8_std" ~ 100,
  )) %>% 
  filter(!is.na(Cq_mean)) 
std_curve_vulgatus_plate_2 %>%
  ggplot(aes(x = copies_ul, y = Cq_mean )) +
  geom_point() +
  theme_bw() + 
  scale_x_continuous(trans = "log10") +
  geom_smooth(method = "lm", se = FALSE) 

# Fit a linear regression model
lm_model_plate_bv_2 <- lm(Cq_mean ~ log10(copies_ul), data = std_curve_vulgatus_plate_2)

# Extract coefficients and R-squared value
lm_summary_plate_bv_2 <- summary(lm_model_plate_bv_2)
intercept_plate_bv_2 <- coef(lm_model_plate_bv_2)[1]
slope_plate_bv_2 <- coef(lm_model_plate_bv_2)[2]
r_squared_plate_bv_2 <- lm_summary_plate_bv_2$r.squared

# Print intercept, slope, and R-squared value
cat("Intercept:", intercept_plate_bv_2, "\n")
cat("Slope:", slope_plate_bv_2, "\n")
cat("R-squared:", r_squared_plate_bv_2, "\n")

# Plot with linear regression line
std_curve_vulgatus_plate_2 %>%
  ggplot(aes(x = copies_ul, y = Cq_mean)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  theme_bw() + 
  scale_x_continuous(trans = "log10") +
  labs(x = "Log10 Copies/ul", y = "Cq mean") +
  geom_text(x = 1e3, y = 35, label = paste("y =", round(intercept_plate_bv_2, 2), "+", round(slope_plate_bv_2, 2), "* log10(x)", "\nR-squared =", round(r_squared_plate_bv_2, 4)), color = "blue")


#Use the standard curve coefficients to calculate the copies per ul for samples
plate2_vulgatus_copies <- results_with_mean %>%
  filter(Cq_mean < Cq_cutoff) %>% 
  mutate(copies_ul = 10^((Cq_mean - NEB_intercept)/NEB_slope)) %>% #Change here if you are calculating std. curve from R
  mutate(copies_ml = copies_ul * 1000) %>% 
  separate(Sample, into = c("Sample", "Timepoint"), sep = "_T") %>% 
  separate(Sample, into = c("Sample", "bio_rep"), sep = "_(?=\\d+$)", extra = "merge", convert = TRUE) %>% 
  group_by(Sample, Timepoint) %>% 
  mutate(Avg_copies_vulgatus_ml = mean(copies_ml, na.rm = TRUE)) %>% 
  mutate(stdev_avg_vulgatus_copies_ml = sd(copies_ml, na.rm = TRUE))



####### Aggregate the dataframes from all plates above #######
plate1_copies <- plate1_copies %>% 
  mutate(Timepoint = as.numeric(Timepoint)) %>% 
  mutate(Sample = as.numeric(Sample)) %>%
  filter((Sample == "4")) %>% 
  filter(!(Sample == "4" & Timepoint == 20)) %>% 
  filter(!(Sample == "4" & Timepoint == 24)) %>% 
  filter(!(Sample == "4" & Timepoint == 36)) %>%
  filter(!(Sample == "4" & Timepoint == 48)) 
plate2_copies <- plate2_copies %>% 
  mutate(Timepoint = as.numeric(Timepoint)) %>% 
  mutate(Sample = as.numeric(Sample)) %>%
  filter((Sample == "4")) %>% 
  filter(!(Sample == "4" & Timepoint == 20)) %>% 
  filter(!(Sample == "4" & Timepoint == 24)) %>% 
  filter(!(Sample == "4" & Timepoint == 36)) %>%
  filter(!(Sample == "4" & Timepoint == 48)) 
rerun_plate_copies_24 <- rerun_plate_copies %>% 
  mutate(Timepoint = as.numeric(Timepoint)) %>% 
  mutate(Sample = as.numeric(Sample)) %>%
  filter((Sample == "4" & Timepoint == 24))
rerun_plate_copies_20 <- rerun_plate_copies %>% 
  mutate(Timepoint = as.numeric(Timepoint)) %>% 
  mutate(Sample = as.numeric(Sample)) %>%
  filter((Sample == "4" & Timepoint == 20))

plate1_copies$Cq_diff <- as.character(plate1_copies$Cq_diff)
plate2_copies$Cq_diff <- as.character(plate2_copies$Cq_diff)
rerun_plate_copies_20$Cq_diff <- as.character(rerun_plate_copies_20$Cq_diff)
rerun_plate_copies_24$Cq_diff <- as.character(rerun_plate_copies_24$Cq_diff)

##########plot sup fig 11F ###############

crass_bulk_ferm_copies <- bind_rows(plate1_copies, plate2_copies, rerun_plate_copies_20, rerun_plate_copies_24)

crass_bulk_ferm_copies_tech_rep_plot <- na.omit(crass_bulk_ferm_copies)

crass_bulk_ferm_copies_tech_rep_plot$tech_replicate <- factor(crass_bulk_ferm_copies_tech_rep_plot$tech_replicate)

tech_rep_colors <- c("1" = "#009292FF", "2" = "black")

bulk_ferm_replicates_crass <- ggplot(crass_bulk_ferm_copies_tech_rep_plot,
                                     aes(x = interaction(Sample, Timepoint, bio_rep), y = Cq, color = tech_replicate)) +
  geom_point(size = 3, alpha = 0.5) +
  labs(x = "Sample.timepoint.biorep", y = "Cq") +
  theme_bw() +
  theme(text = element_text(size = 10)) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        legend.background = element_rect(linetype = 1, 
                                         linewidth = 0.5, 
                                         color = "gray"), 
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  coord_cartesian(ylim = c(0, 30)) +
  scale_color_manual(values = tech_rep_colors) 

bulk_ferm_replicates_crass



##### plot communis copies only####

all_qpcr_data <- crass_bulk_ferm_copies %>% 
  mutate(Timepoint = as.numeric(Timepoint))

for_plotting_df <- all_qpcr_data %>% 
  select(Sample, Timepoint, Avg_copies_ml, stdev_avg_copies_ml) %>% 
  group_by(Sample, Timepoint) %>% 
  distinct(Sample, .keep_all = TRUE) %>% 
  ungroup() %>%
  filter(stdev_avg_copies_ml != 0 | is.na(stdev_avg_copies_ml))

for_plotting_df$Sample <- factor(for_plotting_df$Sample)


abx_communis <- ggplot(for_plotting_df %>% filter(Sample == "4"), 
                       aes(x = Timepoint, y = Avg_copies_ml)) +
  geom_point(aes(color = Sample, shape = Sample), size = 4) +
  labs(title = expression(paste("")),
       x = "Time (Hours)", y = "copies/mL") +
  geom_line(aes(color = Sample), size = 1.5) +
  geom_hline(yintercept = 10000, linetype = "dotted", color = "black", linewidth = 1.5) +
  annotate("text", x = max(for_plotting_df$Timepoint), y = 10000, 
           label = expression(paste("LOD ", italic("C. communis"))), color = "black", hjust = 1, vjust = 2, size = 10) +
  geom_hline(yintercept = 100000, linetype = "dotted", color = "black", , linewidth = 1.5) +
  annotate("text", x = max(for_plotting_df$Timepoint), y = 100000, 
           label = expression(paste("LOD ", italic("P. vulgatus"))), color = "black", hjust = 1, vjust = 2, size = 10) +
  scale_color_manual(values = c("4" = "#009292FF"),
                     name = "qPCR targets", 
                     labels = c("4" = "p-crAss")) +
  scale_shape_manual(values = c("4" = 16),
                     name = "stool sample", 
                     labels = c("4" = "p-crAss +")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        text = element_text(size = 40),
        panel.grid.major = element_line(linewidth = 1.5),
        panel.grid.minor = element_line(linewidth = 1.5),
        axis.line = element_line(linetype = 1, linewidth = 1),
        panel.border = element_rect(linetype = 1, linewidth = 3),
        legend.background = element_rect(linetype = 1, linewidth = 0.5, color = "gray")) +
  geom_errorbar(aes(x = Timepoint, ymin = Avg_copies_ml - stdev_avg_copies_ml, ymax = Avg_copies_ml + stdev_avg_copies_ml, color = Sample),
                size = 0.75, width = 0.5) +
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
                     breaks = c(1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000),
                     labels = c("1000" = "undetected",
                                "10000" = "10000", "100000" = "100000", "1000000" = "1000000",
                                "10000000" = "10000000", "100000000" = "100000000", "1000000000" = "1000000000"))  +
  annotation_logticks(sides = "l")  # Only log axis ticks on the left side (y-axis) 

abx_communis



###############plot sup fig 12###################

plate1_vulgatus_copies <- plate1_vulgatus_copies %>%
  mutate(Sample = as.character(Sample)) %>%
  filter(!(Timepoint == 28))

plate2_vulgatus_copies <- plate2_vulgatus_copies %>%
  mutate(Sample = as.character(Sample)) 
vulgatus_bulk_ferm_copies <- bind_rows(plate1_vulgatus_copies, plate2_vulgatus_copies)


vulgatus_bulk_ferm_copies_tech_rep_plot <- na.omit(vulgatus_bulk_ferm_copies)

vulgatus_bulk_ferm_copies_tech_rep_plot$tech_replicate <- factor(vulgatus_bulk_ferm_copies_tech_rep_plot$tech_replicate)

vulgatus_bulk_ferm_copies_tech_rep_plot <- vulgatus_bulk_ferm_copies_tech_rep_plot %>%
  filter(Sample == "4")

tech_rep_colors <- c("1" = "#009292FF", "2" = "black")


bulk_ferm_replicates_crass <- ggplot(vulgatus_bulk_ferm_copies_tech_rep_plot,
                                     aes(x = interaction(Sample, Timepoint, bio_rep), y = Cq, color = tech_replicate)) +
  geom_point(size = 3, alpha = 0.5) +
  labs(x = "Sample.timepoint.biorep", y = "Cq") +
  theme_bw() +
  theme(text = element_text(size = 10)) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        legend.background = element_rect(linetype = 1, 
                                         linewidth = 0.5, 
                                         color = "gray"), 
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  coord_cartesian(ylim = c(0, 30)) +
  scale_color_manual(values = tech_rep_colors) 

bulk_ferm_replicates_crass


#######plot vulgatus copies in bulk ferm ############
for_plotting_vulgatus_df <- vulgatus_bulk_ferm_copies %>% 
  select(Sample, Timepoint, Avg_copies_vulgatus_ml, stdev_avg_vulgatus_copies_ml) %>% 
  group_by(Sample, Timepoint) %>% 
  distinct(Sample, .keep_all = TRUE) %>% 
  ungroup() %>% 
  filter(stdev_avg_vulgatus_copies_ml != 0 | is.na(stdev_avg_vulgatus_copies_ml))

for_plotting_vulgatus_df <- for_plotting_vulgatus_df %>%
  mutate(Timepoint = as.numeric(Timepoint))


abx_vulgatus <- ggplot(for_plotting_vulgatus_df %>% filter(Sample == "4"), 
                       aes(x = Timepoint, y = Avg_copies_vulgatus_ml)) +
  geom_point(aes(color = Sample, shape = Sample), size = 4) +
  labs(title = expression(paste("")),
       x = "Time (Hours)", y = "copies/mL") +
  geom_line(aes(color = Sample, group = Sample), size = 1.5) +
  geom_hline(yintercept = 10000, linetype = "dotted", color = "black", linewidth = 1.5) +
  annotate("text", x = max(for_plotting_vulgatus_df$Timepoint), y = 10000, 
           label = expression(paste("LOD ", italic("C. communis"))), color = "black", hjust = 1, vjust = 2, size = 10) +
  geom_hline(yintercept = 100000, linetype = "dotted", color = "black", , linewidth = 1.5) +
  annotate("text", x = max(for_plotting_vulgatus_df$Timepoint), y = 100000, 
           label = expression(paste("LOD ", italic("P. vulgatus"))), color = "black", hjust = 1, vjust = 2, size = 10) +
  scale_color_manual(values = c("4" = "#009292FF"),
                     name = "qPCR targets", 
                     labels = c("4" = "p-crAss")) +
  scale_shape_manual(values = c("4" = 16, "6" = 2),
                     name = "stool sample", 
                     labels = c("4" = "p-crAss +")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        text = element_text(size = 40),
        panel.grid.major = element_line(linewidth = 1.5),
        panel.grid.minor = element_line(linewidth = 1.5),
        axis.line = element_line(linetype = 1, linewidth = 1),
        panel.border = element_rect(linetype = 1, linewidth = 3),
        legend.background = element_rect(linetype = 1, linewidth = 0.5, color = "gray")) +
  geom_errorbar(aes(x = Timepoint, ymin = Avg_copies_vulgatus_ml - stdev_avg_vulgatus_copies_ml, ymax = Avg_copies_vulgatus_ml + stdev_avg_vulgatus_copies_ml, color = Sample),
                size = 0.75, width = 0.5) +
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
                     breaks = c(1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000),
                     labels = c("1000" = "undetected",
                                "10000" = "10000", "100000" = "100000", "1000000" = "1000000",
                                "10000000" = "10000000", "100000000" = "100000000", "1000000000" = "1000000000"))  +
  annotation_logticks(sides = "l")  # Only log axis ticks on the left side (y-axis) 

abx_vulgatus


#######plot sup fig 6B ############

combined_plot <- ggplot() +
  geom_point(data = for_plotting_df %>% filter(Sample == "4"), 
             aes(x = Timepoint, y = Avg_copies_ml, color = Sample, shape = Sample), size = 4) +
  geom_line(data = for_plotting_df %>% filter(Sample == "4"), 
            aes(x = Timepoint, y = Avg_copies_ml, color = Sample), size = 1.5) +
  geom_errorbar(data = for_plotting_df %>% filter(Sample == "4"),
                aes(x = Timepoint, ymin = Avg_copies_ml - stdev_avg_copies_ml, ymax = Avg_copies_ml + stdev_avg_copies_ml, color = Sample),
                size = 0.75, width = 0.5) +
  geom_point(data = for_plotting_vulgatus_df %>% filter(Sample == "4"), 
             aes(x = Timepoint, y = Avg_copies_vulgatus_ml, shape = Sample), size = 4, color = "black") +
  geom_line(data = for_plotting_vulgatus_df %>% filter(Sample == "4"), 
            aes(x = Timepoint, y = Avg_copies_vulgatus_ml), color = "black", size = 1.5, linetype = "solid", group = 1) +
  geom_errorbar(data = for_plotting_vulgatus_df %>% filter(Sample == "4"),
                aes(x = Timepoint, ymin = Avg_copies_vulgatus_ml - stdev_avg_vulgatus_copies_ml, ymax = Avg_copies_vulgatus_ml + stdev_avg_vulgatus_copies_ml),
                size = 0.75, width = 0.5, color = "black") +
  labs(title = "",
       x = "Time (Hours)", y = "copies/mL",
       color = "qPCR targets", shape = "stool sample") +
  scale_color_manual(values = c("4" = "#009292FF"),
                     labels = c("4" = "p-crAss")) +
  scale_shape_manual(values = c("4" = 16),
                     labels = c("4" = "p-crAss +")) +
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
                     breaks = c(1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000),
                     labels = c("1000" = "1000",
                                "10000" = "10000", "100000" = "100000", "1000000" = "1000000",
                                "10000000" = "10000000", "100000000" = "100000000", "1000000000" = "1000000000")) +
  annotation_logticks(sides = "l") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        text = element_text(size = 40),
        panel.grid.major = element_line(linewidth = 1.5),
        panel.grid.minor = element_line(linewidth = 1.5),
        axis.line = element_line(linetype = 1, linewidth = 1),
        panel.border = element_rect(linetype = 1, linewidth = 3),
        legend.background = element_rect(linetype = 1, linewidth = 0.5, color = "gray"))

combined_plot
pdf(here("~/Documents/240307_p4100_bulk_ferm_abx.pdf"), width = 14, height = 12)
print(combined_plot)
dev.off()



#####calculate phage:host ratio#####

# Merge the dataframes on Sample and Timepoint
merged_data <- merge(for_plotting_df, for_plotting_vulgatus_df, by = c("Sample", "Timepoint"))


# Calculate the ratio column
final_merged_ratio <- merged_data %>%
  mutate(ratio = Avg_copies_ml / Avg_copies_vulgatus_ml)

##### plot phage:host ratio####

final_merged_ratio <- final_merged_ratio %>%
  mutate(Timepoint = as.numeric(Timepoint)) %>% 
  filter(stdev_avg_vulgatus_copies_ml != 0 | is.na(stdev_avg_vulgatus_copies_ml))


phage_host_ratio_plot <- ggplot(final_merged_ratio %>% filter(Sample == "4"), 
                                aes(x = Timepoint, y = ratio)) +
  geom_point(aes(color = Sample, shape = Sample), size = 4) +
  labs(title = expression(paste("")),
       x = "Time (Hours)", y = "copies/mL") +
  geom_line(aes(color = Sample), size = 1.5) +
  scale_color_manual(values = c("4" = "#009292FF"),
                     name = "qPCR targets", 
                     labels = c("4" = "p-crAss")) +
  scale_shape_manual(values = c("4" = 16),
                     name = "stool sample", 
                     labels = c("4" = "p-crAss +")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        text = element_text(size = 40),
        panel.grid.major = element_line(linewidth = 1.5),
        panel.grid.minor = element_line(linewidth = 1.5),
        axis.line = element_line(linetype = 1, linewidth = 1),
        panel.border = element_rect(linetype = 1, linewidth = 3),
        legend.background = element_rect(linetype = 1, linewidth = 0.5, color = "gray")) +
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
                     breaks = c(0.1, 1, 10, 100, 1000))

phage_host_ratio_plot



#######Plot Fig 7C make combined plot of phage:host ratio for vulgatus, dorei, and bulk ferm#####
##MUST RUN ISOLATE CODE FOR THIS FIRST##

combined_plot_phage_host_ratio <- ggplot() +
  geom_point(data = final_ratio %>% filter(Sample == "Bv_phage" | Sample == "Bd_phage"),
             aes(x = Timepoint, y = ratio, color = Sample, shape = Sample), size = 4) +
  geom_line(data = final_ratio %>% filter(Sample == "Bv_phage" | Sample == "Bd_phage"),
            aes(x = Timepoint, y = ratio, color = Sample), size = 1.5) +
  geom_point(data = final_merged_ratio %>% filter(Sample == "4"),
             aes(x = Timepoint, y = ratio, color = Sample, shape = Sample), size = 4) +
  geom_line(data = final_merged_ratio %>% filter(Sample == "4"),
            aes(x = Timepoint, y = ratio, color = Sample), size = 1.5) +
  scale_color_manual(values = c("Bv_phage" = "black", "Bd_phage" = "#009292FF", "4" = "gray"),
                     name = "qPCR targets", 
                     labels = c("Bv_phage" = "Bv_phage", "Bd_phage" = "Bd_phage", "4" = "p-crAss")) +
  scale_shape_manual(values = c("Bv_phage" = 16, "Bd_phage" = 16, "4" = 16),
                     name = "stool sample", 
                     labels = c("Bv_phage" = "Bv_phage", "Bd_phage" = "Bd_phage", "4" = "p-crAss +")) +
  labs(title = expression(paste("")),
       x = "Time (Hours)", y = "copies/mL") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        text = element_text(size = 40),
        panel.grid.major = element_line(linewidth = 1.5),
        panel.grid.minor = element_line(linewidth = 1.5),
        axis.line = element_line(linetype = 1, linewidth = 1),
        panel.border = element_rect(linetype = 1, linewidth = 3),
        legend.background = element_rect(linetype = 1, linewidth = 0.5, color = "gray")) +
  scale_y_log10(expand = expand_scale(mult = c(0.01, 0.1)),
                labels = scales::label_number())

print(combined_plot_phage_host_ratio)
pdf(here("~/Documents/240307_combined_plot_phage_host_ratio.pdf"), width = 13, height = 10)
print(combined_plot_phage_host_ratio)
dev.off()

