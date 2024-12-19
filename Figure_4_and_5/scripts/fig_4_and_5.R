'''
Final script for generating figure 3 panels for Schmidtke et al. 2024
Date: 26 Aug 2024 
Authors: Angela Hickey, Danica Schmidtke 
'''

library(tidyverse)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(gggenes)
library(here)
library(readxl)
library(readr)
##Figure 3A and B are generated in a Google Colab notebook

######## Figure 3C: short read coverage plot##########################
# Read the Excel file
data <- read_excel("~/Documents/B1_reads_align_to_pcrass_assembly_linear_at_DNA_ligase_coverage.xlsx")

# Calculate the median coverage
median_coverage <- median(data$final_coverage, na.rm = TRUE)


# Calculate the fold change in coverage relative to the median
data$coverage_fold_change <- data$final_coverage / median_coverage

#Find the region of the genome where the coverage is high (the bounds of the DTR)
spike_in_cov <- data %>% 
  filter(coverage_fold_change > 5)

# Plot the fold change in coverage
cov_plot <- ggplot(data, aes(x = reverse_position, y = coverage_fold_change)) +
  #geom_point(size = 0.1, color = "lightgrey") +  # Adjust the size of the points
  geom_line() +
  labs(title = "B1",
       x = "Distance from DNA ligase",
       y = "Fold Change in Coverage") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_bw()
cov_plot

pdf(here("~/Documents/SR_coverage_Hi_C.pdf"), width = 8, height = 5)
print(cov_plot)
dev.off()





######## Figure 3 supplement short read, read ends plot##########################

#plot 5 prime ends
ends_5_plot <- ggplot(data, aes(x = reverse_position, y = final_5_ends_over_cov)) +
  geom_point() +
  labs(title = "B1",
       x = "Distance from DNA ligase gene",
       y = "5 prime ends") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  bw_theme
ends_5_plot


#plot 3 prime ends
ends_3_plot <- ggplot(data, aes(x = reverse_position, y = final_3_ends_over_cov)) +
  geom_point() +
  labs(title = "B1",
       x = "Distance from DNA ligase gene",
       y = "5 prime ends") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  bw_theme
ends_3_plot


#plot all ends
read_ends_plot <- ggplot(data, aes(x = reverse_position, y = final_total_ends_over_cov)) +
  geom_point() +
  labs(title = "B1",
       x = "Distance from DNA ligase gene",
       y = "read ends") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  bw_theme
read_ends_plot

pdf(here("~/Documents/SR_read_ends_Hi_C.pdf"), width = 8, height = 5)
print(read_ends_plot)
dev.off()

######## Figure 3D: Publicly available crassphage genome sequencing depth  #########

### Read in Bakta information ###
# Set your working directory to the folder containing the files
setwd("/labs/asbhatt/wirbel/projects/other_projects/collaboration_pcrass")

# List all directories in the main folder
bakta_folder_list <- list.dirs(full.names = TRUE, recursive = FALSE)

# Create an empty list to store data frames
bakta_df_list <- list()

# Loop through each directory
for (folder in bakta_folder_list) {
  file_list <- list.files(path = paste0(folder, "/pcrass_bakta/"), pattern = "*.gff3", full.names = TRUE)
  
  for (file in file_list) {
    # Read the TSV file into a data frame
    files_df <- read.delim(file, header = FALSE, sep = "\t", col.names = c("genome", "origin", "type", "start", "end", "dot", "strand", "idk", "info"), skip = 8)
    
    # Add the data frame to the list
    bakta_df_list[[str_remove(file, '/pcrass_bakta.*')]] <- files_df
  }
}

#parse through the bakta files in the list to get the product name 
for (i in seq_along(bakta_df_list)) {
  df <- bakta_df_list[[i]] %>%
    filter(type == "CDS") %>%
    separate(info, c("ID", "Name", "locus_tag", "product", "other"), sep = ";") %>%
    mutate(product = str_remove(product, ".*=" )) #%>%
  #select(start, end, strand, product) #removing this because with the blastp seq
  #similarity I want to have the locus ID 
  
  bakta_df_list[[i]] <- df
}

bakta_all <- bind_rows(bakta_df_list, .id = "filename") %>% 
  mutate(filename=str_remove(filename, './')) %>% 
  #add here some code to add a column for query that has only the locus ID so I can bind by this 
  separate(ID, c("ID", "query"), sep="=")

### Read in the depth data ###

# Set your working directory to the folder containing the files
setwd("/labs/asbhatt/wirbel/projects/other_projects/collaboration_pcrass")

# List all directories in the main folder
depth_folder_list <- list.dirs(full.names = TRUE, recursive = FALSE)

# Create an empty list to store data frames
depth_df_list <- list()

# Loop through each directory
for (folder in depth_folder_list) {
  file_list <- list.files(path = folder, pattern = "pcrass_contig_depth_2.txt", full.names = TRUE)
  
  for (file in file_list) {
    # Read the TSV file into a data frame
    files_df <- read.delim(file, header = FALSE, sep = "\t", col.names = c("contig", "position", "depth"))
  
    
    files_df$contig <- as.character(files_df$contig)
    files_df$position <- as.numeric(files_df$position)
    files_df$depth <- as.numeric(files_df$depth)
    
    
    # Add the data frame to the list
    depth_df_list[[str_remove(file, '/pcrass.*')]] <- files_df
  }
}

depth_all <- bind_rows(depth_df_list, .id = "filename") %>% 
  mutate(filename=str_remove(filename, './')) 

### Linearize position to the same protein and plot coverage ###

#normalize to the DNA ligase protein called by Bakta  
normalization_protein <- bakta_all %>% 
  filter(product == "DNA ligase") %>% 
  select(filename, start, end, strand)

#Calculate the total contig length for repositioning the genomes
contig_length <- depth_all %>%
  group_by(filename) %>%
  summarise(total_contig_length = max(position))

#calculate the median coverage per sample to normalize 
median_coverage <- depth_all %>% 
  group_by(filename) %>% 
  summarise(median_coverage = median(depth))

# calculate the average depth over 100bp chunks and create new positioning
# based on the normalized protein location and strand 
# also calculate normalized depth to median coverage per contig 
average_depth_100 <- depth_all %>%
  group_by(filename, chunk = ceiling(position/100)) %>%
  summarise(avg_depth = mean(depth)) %>% 
  mutate(position = chunk*100) %>% 
  ungroup() %>% 
  left_join(contig_length) %>% 
  left_join(median_coverage) %>% 
  left_join(normalization_protein) %>% 
  group_by(filename) %>% 
  mutate(offset = case_when(strand == "+" ~ end, strand == "-" ~ start)) %>% 
  mutate(new_position = case_when(
    sign(position - offset) == 1 ~ (position - offset),
    sign(position - offset) == -1 ~ (position - offset + total_contig_length),
    TRUE ~ position
  )) %>% 
  mutate(normalized_depth = avg_depth/median_coverage) %>% 
  mutate(stranded_position = case_when(strand == "+" ~ total_contig_length - new_position, strand == "-" ~ new_position)) %>% 
  mutate(which_strand = case_when(strand == "+" ~ "pos", strand == "-" ~ "neg"))

#Calculate the mean depth across contigs
mean_depth <- average_depth_100 %>%
  filter(filename != "CD95") %>% #filtering because this sample has a misassembled crassphage genome
  group_by(new_chunk = ceiling(stranded_position/100)) %>%
  summarise(normalized_depth = mean(normalized_depth)) %>% 
  mutate(stranded_position = new_chunk*100) %>% 
  mutate(filename = "mean")

#add the mean depth to the existing dataset 
df_with_mean <- average_depth_100 %>% 
  bind_rows(mean_depth)

#putative origin locations: 
  #(21791 , 21674) 
  #(21133 , 21125) 
  #(20983 , 20897) 
  #(20770 , 20468)
  #(17866 , 17054)

# Define the pairs of numbers
pairs <- list(c(21791, 21674), c(21133, 21125), c(20983, 20897), c(20770, 20468), c(17866, 17054))

# Calculate the midpoints
midpoints <- sapply(pairs, function(x) mean(x))

# Print the midpoints
midpoints

#21129,20940,
# Midpoints and their respective colors
midpoints <- c(21732.5, 20619, 17460)
colors <- c("blue", "green", "purple")

# Create the plot
c <- df_with_mean %>% 
  filter(filename != "CD95") %>%
  ggplot(aes(x = stranded_position, y = normalized_depth, color = filename)) +
  geom_line(data = . %>% filter(filename != "mean"), aes(group = filename), color = "lightgrey") +
  geom_line(data = . %>% filter(filename == "mean"), aes(group = filename), color = "black", size = 1.5) +
  theme_bw() +
  xlab("Genome Position") + 
  ylab("Fold change in coverage per contig") +
  scale_color_manual(values = c("mean" = "black")) +
  # Add vertical lines with different colors for each midpoint
  geom_vline(data = data.frame(midpoints, colors), 
             aes(xintercept = midpoints, color = colors), 
             linetype = "dashed", size = 0.5) +
  scale_color_manual(values = setNames(colors, colors)) +  # Manually map colors to the 'color' aesthetic
  guides(size = FALSE)

# Display the plot
c

######### What are the position locations of the double coverage region #####

#calculate the fold change coverage per base 
per_base_depth <- depth_all %>% 
  filter(filename != "CD95") %>% 
  left_join(median_coverage) %>%
  mutate(normalized_coverage = depth/median_coverage)

first_position<- per_base_depth %>%
  group_by(filename) %>% 
  # Find the first position where "depth" exceeds 2
  filter(case_when(filename == "TD47" ~ normalized_coverage > 3, filename != "TD47" ~ normalized_coverage > 2.4)) %>%
  slice(1) %>% 
  ungroup() %>% 
  mutate(dbl_start = position)

last_position <- per_base_depth %>%
  group_by(filename) %>% 
  filter(case_when(filename == "TD47" ~ normalized_coverage > 3, filename != "TD47" ~ normalized_coverage > 2.4)) %>%
  slice(n()) %>% 
  ungroup() %>% 
  mutate(dbl_end = position)

dbl_cov_positions <- first_position %>%
  select(filename, dbl_start) %>% 
  left_join(last_position) %>% 
  select(filename, dbl_start, dbl_end)

dbl_size_df <- dbl_cov_positions %>% mutate(double_cov_size = dbl_end - dbl_start)
dbl_size_df %>% summarise(mean(double_cov_size))

df_with_mean %>% 
  filter(filename != "CD95") %>%
  # filter(filename == "CD73") %>% 
  ggplot(aes(x = position, y = normalized_depth, color = filename)) +
  geom_line(data = . %>% filter(filename != "mean"), aes(group = filename), color = "lightgrey") +
  geom_line(data = . %>% filter(filename == "mean"), aes(group = filename), color = "black", size = 1.5) +
  theme_bw() +
  xlab("Genome Position") + 
  ylab("Fold change in coverage per contig") +
  scale_color_manual(values = c("mean" = "black")) +
  guides(size = FALSE) + 
  geom_vline(data = dbl_cov_positions, aes(xintercept = dbl_start), linetype = "dotted") +
  geom_vline(data = dbl_cov_positions, aes(xintercept = dbl_end), linetype = "dotted") +
  facet_wrap(~filename)# This removes the size legend
# scale_color_identity(guide = "none")



test_adjusted_cov <- average_depth_100 %>% left_join(dbl_cov_positions) %>% mutate(adjusted_cov = case_when((position >= dbl_start & position <= dbl_end) ~ normalized_depth/2, TRUE ~ normalized_depth ))



c <- test_adjusted_cov %>%
       filter(filename != "CD95") %>%
       ggplot(aes(x = stranded_position, y = adjusted_cov, color = filename)) +
       geom_line(data = . %>% filter(filename != "mean"), aes(group = filename), color = "lightgrey") +
       geom_line(data = . %>% filter(filename == "mean"), aes(group = filename), color = "black", size = 1.5) +
       theme_bw() +
       xlab("Genome Position") +
       ylab("Fold change in coverage per contig") +
       scale_color_manual(values = c("mean" = "black")) +
       geom_vline(data = data.frame(midpoints, colors), 
                                   aes(xintercept = midpoints, color = colors), 
                                   linetype = "dashed", size = 0.5) +
       scale_color_manual(values = setNames(colors, colors)) +
       guides(size = FALSE) +
       facet_wrap(~ filename, scales = "free_y")
c


ggsave(c , filename = "/labs/asbhatt/danicats/241107_coverage_pcrass.pdf", width = 10, height = 7)


######## Figure 3E: Read ends for publicly available crassphage genomes #######

## This code relies on the following dataframes generated by the Figure 3D code
# depth_all
# normalization_protein

### MAIN FIGURE CODE ###
# Reading in the read end counts from the samclip processed data 
# this removes the soft clipped reads to count the "native" ends 

# Set your working directory to the folder containing the files
setwd("/labs/asbhatt/angelah/projects/pcrass_project/integration_long_read/assembly_mapping/")

# List all directories in the main folder
readend_folder_list <- list.dirs(full.names = TRUE, recursive = FALSE)

# Create an empty list to store data frames
three_prime_df_list <- list()
five_prime_df_list <- list()

# Loop through each directory
for (folder in readend_folder_list) {
  three_prime_list <- list.files(path = folder, pattern = "3_prime_ends_samclip.txt", full.names = TRUE)
  five_prime_list <- list.files(path = folder, pattern = "5_prime_ends_samclip.txt", full.names = TRUE)
  
  for (file in three_prime_list) {
    # Skip files with "CD95" in the name
    if (grepl("CD95", file)) {
      next
    }
    
    # Read the TSV file into a data frame
    three_prime_df <- read.table(file, header = FALSE, col.names = c("contig", "position", "three_end"), sep = "\t", stringsAsFactors = FALSE)
    
    
    three_prime_df$contig <- as.character(three_prime_df$contig)
    three_prime_df$position <- as.numeric(three_prime_df$position)
    three_prime_df$three_end <- as.numeric(three_prime_df$three_end)
    
    
    # Add the data frame to the list
    three_prime_df_list[[str_remove(file, '/3.*')]] <- three_prime_df
  }
  
  for (file in five_prime_list) {
    # Skip files with "CD95" in the name
    if (grepl("CD95", file)) {
      next
    }
    
    # Read the TSV file into a data frame
    five_prime_df <- read.table(file, header = FALSE, col.names = c("contig", "position", "five_end"), sep = "\t", stringsAsFactors = FALSE)
    
    
    five_prime_df$contig <- as.character(five_prime_df$contig)
    five_prime_df$position <- as.numeric(five_prime_df$position)
    five_prime_df$five_end <- as.numeric(five_prime_df$five_end)
    
    
    # Add the data frame to the list
    five_prime_df_list[[str_remove(file, '/5.*')]] <- five_prime_df
  }
}

three_prime_samclip <- bind_rows(three_prime_df_list, .id = "filename") %>% 
  mutate(filename=str_remove(filename, './')) 

five_prime_samclip <- bind_rows(five_prime_df_list, .id = "filename") %>% 
  mutate(filename=str_remove(filename, './'))

readends_all_notfixed <- three_prime_samclip %>%   
  full_join(five_prime_samclip)

#CD41 pcrass contig is 266 but the grep to generate the read ends file grabs other things that match 266 
#to fix that I am going to merge in the corrected CD41 read ends calculation 

correction_df <- readends_all_notfixed %>% 
  filter(filename == "CD41" & contig == "contig_266")

readends_samclip <- readends_all_notfixed %>% 
  filter(filename != "CD41") %>% 
  bind_rows(correction_df)

contig_length_samclip <- readends_samclip %>%
  group_by(filename) %>%
  summarise(total_contig_length = max(position))


every_position_readends_samclip <- readends_samclip %>% 
  left_join(contig_length_samclip) %>% 
  left_join(normalization_protein) %>% 
  group_by(filename) %>% 
  mutate(offset = case_when(strand == "+" ~ end, strand == "-" ~ start)) %>% 
  mutate(new_position = case_when(
    sign(position - offset) == 1 ~ (position - offset),
    sign(position - offset) == -1 ~ (position - offset + total_contig_length),
    TRUE ~ position
  )) %>% 
  mutate(stranded_position = case_when(strand == "+" ~ total_contig_length - new_position, strand == "-" ~ new_position)) %>% 
  mutate(which_strand = case_when(strand == "+" ~ "pos", strand == "-" ~ "neg"))


every_position_readends_depth_samclip <- every_position_readends_samclip %>% 
  left_join(depth_all) %>% 
  mutate(total_ends = three_end + five_end) %>% 
  mutate(normalized_total_ends = total_ends/depth) %>% 
  mutate(normalized_three_end = three_end/depth) %>% 
  mutate(normalized_five_end = five_end/depth)

#Calculate the mean read ends across contigs at each position 
mean_read_ends <- every_position_readends_depth_samclip %>%
  filter(!(position <= 20 | position >= (total_contig_length - 20))) %>% 
  group_by(stranded_position) %>% 
  filter(filename != "CD95") %>% #filtering because this sample has a misassembled crassphage genome
  summarise(normalized_total_ends = mean(normalized_total_ends)) %>% 
  mutate(filename = "mean") %>% 
  ungroup()

readends_samclip_mean <- every_position_readends_depth_samclip %>% 
  filter(!(position <= 20 | position >= (total_contig_length -20))) %>% 
  bind_rows(mean_read_ends)

plot <- readends_samclip_mean %>% 
  filter(filename != "CD95") %>%
  ggplot(aes(x = stranded_position, y = normalized_total_ends, color = filename)) +
  geom_line(data = . %>% filter(filename != "mean"), aes(group = filename), color = "lightgrey") +
  geom_line(data = . %>% filter(filename == "mean"), aes(group = filename), color = "black", size = 1.5) +
  theme_bw() +
  scale_color_manual(values = c("mean" = "black"))
plot
ggsave( plot,  filename = "/Users/angelahickey/Desktop/240826_mean_read_ends_samclip.pdf", width = 10, height = 7)

### SUPPLEMENTARY FIGURE ###
# This reads in the read ends calculated from the alignment before any filtering of reads 
# Set your working directory to the folder containing the files
setwd("/labs/asbhatt/angelah/projects/pcrass_project/integration_long_read/assembly_mapping/")

# List all directories in the main folder
readend_folder_list <- list.dirs(full.names = TRUE, recursive = FALSE)

# Create an empty list to store data frames
three_prime_df_list <- list()
five_prime_df_list <- list()

# Loop through each directory
for (folder in readend_folder_list) {
  three_prime_list <- list.files(path = folder, pattern = "3_prime_ends.txt", full.names = TRUE)
  five_prime_list <- list.files(path = folder, pattern = "5_prime_ends.txt", full.names = TRUE)
  
  for (file in three_prime_list) {
    # Skip files with "CD95" in the name
    if (grepl("CD95", file)) {
      next
    }
    
    # Read the TSV file into a data frame
    three_prime_df <- read_tsv(file, col_names = c("contig", "position", "three_end"))
    
    three_prime_df$contig <- as.character(three_prime_df$contig)
    three_prime_df$position <- as.numeric(three_prime_df$position)
    three_prime_df$three_end <- as.numeric(three_prime_df$three_end)
    
    
    # Add the data frame to the list
    three_prime_df_list[[str_remove(file, '/3.*')]] <- three_prime_df
  }
  
  for (file in five_prime_list) {
    # Skip files with "CD95" in the name
    if (grepl("CD95", file)) {
      next
    }
    
    # Read the TSV file into a data frame
    five_prime_df <- read_tsv(file, col_names = c("contig", "position", "five_end"))
    
    five_prime_df$contig <- as.character(five_prime_df$contig)
    five_prime_df$position <- as.numeric(five_prime_df$position)
    five_prime_df$five_end <- as.numeric(five_prime_df$five_end)
    
    
    # Add the data frame to the list
    five_prime_df_list[[str_remove(file, '/5.*')]] <- five_prime_df
  }
}

three_prime_all <- bind_rows(three_prime_df_list, .id = "filename") %>% 
  mutate(filename=str_remove(filename, './')) 

five_prime_all <- bind_rows(five_prime_df_list, .id = "filename") %>% 
  mutate(filename=str_remove(filename, './'))

readends_all_notfixed <- three_prime_all %>%   
  full_join(five_prime_all)

#CD41 pcrass contig is 266 but the grep to generate the read ends file grabs other things that match 266 
#to fix that I am going to merge in the corrected CD41 read ends calculation 

correction_df <- readends_all_notfixed %>% 
  filter(filename == "CD41" & contig == "contig_266")

readends_all <- readends_all_notfixed %>% 
  filter(filename != "CD41") %>% 
  bind_rows(correction_df)

every_position_readends <- readends_all %>% 
  left_join(contig_length) %>% 
  left_join(normalization_protein) %>% 
  group_by(filename) %>% 
  mutate(offset = case_when(strand == "+" ~ end, strand == "-" ~ start)) %>% 
  mutate(new_position = case_when(
    sign(position - offset) == 1 ~ (position - offset),
    sign(position - offset) == -1 ~ (position - offset + total_contig_length),
    TRUE ~ position  # You might want to define what to do for other cases
  )) %>% 
  #mutate(normalized_depth = avg_depth/median_coverage) %>% 
  mutate(stranded_position = case_when(strand == "+" ~ total_contig_length - new_position, strand == "-" ~ new_position)) %>% 
  mutate(which_strand = case_when(strand == "+" ~ "pos", strand == "-" ~ "neg"))

every_position_readends_depth <- every_position_readends %>% 
  left_join(depth_all) %>% 
  mutate(total_ends = three_end + five_end) %>% 
  mutate(normalized_total_ends = total_ends/depth) %>% 
  mutate(normalized_three_end = three_end/depth) %>% 
  mutate(normalized_five_end = five_end/depth)

#Calculate the mean read ends across contigs at each position 
mean_read_ends <- every_position_readends_depth %>%
  filter(!(position <= 20 | position >= (total_contig_length - 20))) %>% 
  group_by(stranded_position) %>% 
  filter(filename != "CD95") %>% #filtering because this sample has a misassembled crassphage genome
  summarise(normalized_total_ends = mean(normalized_total_ends)) %>% 
  mutate(filename = "mean") %>% 
  ungroup()

readends_mean <- every_position_readends_depth %>% 
  filter(!(position <= 20 | position >= (total_contig_length -20))) %>% 
  bind_rows(mean_read_ends)

plot <- readends_mean %>% 
  filter(filename != "CD95") %>%
  ggplot(aes(x = stranded_position, y = normalized_total_ends, color = filename)) +
  geom_line(data = . %>% filter(filename != "mean"), aes(group = filename), color = "lightgrey") +
  geom_line(data = . %>% filter(filename == "mean"), aes(group = filename), color = "black", size = 1.5) +
  theme_bw() +
  scale_color_manual(values = c("mean" = "black"))
plot

# plot <- every_position_readends_depth %>% 
#   filter(filename != "CD95") %>%
#   filter(!(position <= 20 | position >= (total_contig_length -20))) %>% 
#   ggplot(aes(x = stranded_position, color = filename)) +
#   geom_point(aes(y = normalized_total_ends), size = 3) +
#   theme_bw() 
# plot
ggsave( plot,  filename = "/Users/angelahickey/Desktop/240826_mean_read_ends_supplement.pdf", width = 10, height = 7)


######## Figure 3F: Reads that span the double coverage region per sample ######
# Read in the span ratio data 

#this file is generated by jupyter notebook analysis of chen et al. long reads aligning to pcrass
span_ratio_counts <- read_tsv("~/Downloads/span_ratio_counts.tsv", 
                              col_names = c("sample", "read_span", "read_ends", "read_starts", 
                                            "total_reads", "ratio_span", "ratio_ends_starts")) %>% 
  mutate(total_reads_without_span = total_reads - read_span)


long_span_ratio_counts <- span_ratio_counts %>% 
  pivot_longer(-sample) %>% 
  mutate(type=factor(name, 
                     levels = c("read_span", "read_starts", "read_ends", "total_reads", "ratio_span", "ratio_ends_starts", "total_reads_without_span")))


# Plot the ratio of reads that span the dbl cov region by sample
long_span_ratio_counts %>% 
  filter(str_detect(type, "ratio_span")) %>% 
  ggplot(aes(x = fct_reorder(sample, value), y = value*100, fill = type)) +
  geom_col() + 
  theme_bw() + 
  scale_fill_manual(values = c("ratio_span" = "#009292FF")) +
  xlab("Sample ID") + 
  ylab("percent reads spanning double coverage region") + 
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none") 

### SUPPLEMENTARY FIGURE ###
# plot the total number of reads aligned to the double coverage - fill them by spanning or not 
long_span_ratio_counts %>% 
  filter(type == "total_reads_without_span" | type == "read_span") %>% 
  ggplot(aes(x = sample, y = value, fill = type)) +
  geom_col() + 
  theme_bw() + 
  scale_fill_manual(values = c("read_span" = "red", "total_reads_without_span" = "#007C92")) +
  xlab("Sample") + 
  ylab("# of reads") + 
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


############## 11/21/24 new linear reference coverage plot ##########

#read in the new depth files generated on the linear reference

# Set your working directory to the folder containing the files
setwd("/labs/asbhatt/angelah/projects/pcrass_project/integration_long_read/pcrass_contigs/new_ref_output/")

# List all directories in the main folder
new_depth_folder_list <- list.dirs(full.names = TRUE, recursive = FALSE)

# Create an empty list to store data frames
new_depth_df_list <- list()

# Loop through each directory
for (folder in new_depth_folder_list) {
  file_list <- list.files(path = folder, pattern = ".*linear_ref_depth\\.txt", full.names = TRUE)
  
  for (file in file_list) {
    # Read the TSV file into a data frame
    files_df <- read.table(file, 
                           header = FALSE, 
                           col.names = c("contig", "position", "depth"), 
                           sep = "\t")
    
    files_df$contig <- as.character(files_df$contig)
    files_df$position <- as.numeric(files_df$position)
    files_df$depth <- as.numeric(files_df$depth)
    
    
    # Add the data frame to the list
    new_depth_df_list[[str_remove(file, '/pcrass.*')]] <- files_df
  }
}

new_depth_all <- bind_rows(new_depth_df_list, .id = "filename") %>% 
  mutate(filename=str_remove(filename, './')) %>% 
  mutate(filename=str_remove(filename, "/.*$"))

#from the non-linear coverage plot, we know which references are in the same orientation based on dna ligase strand information
#here we can extract which sample has the DNA ligase on which strand
strand_info <- average_depth_100 %>% select(filename, which_strand) %>% distinct

#calculate the median coverage per sample to normalize 
new_median_coverage <- new_depth_all %>% 
  group_by(filename) %>% 
  summarise(median_coverage = median(depth))

new_contig_length <- depth_all %>%
  group_by(filename) %>%
  summarise(total_contig_length = max(position))

new_average_depth_100 <- new_depth_all %>%
  group_by(filename, chunk = ceiling(position/100)) %>%
  summarise(avg_depth = mean(depth)) %>% 
  mutate(position = chunk*100) %>% 
  ungroup() %>% 
  left_join(new_median_coverage) %>% 
  left_join(strand_info) %>%
  left_join(new_contig_length) %>% 
  group_by(filename) %>% 
  mutate(normalized_depth = avg_depth/median_coverage) %>% 
  mutate(stranded_position = case_when(which_strand == "neg" ~ total_contig_length - position, which_strand == "pos" ~ position))

#calculate the mean depth
new_mean_depth <- new_average_depth_100 %>%
  filter(filename != "CD95") %>% 
  group_by(new_chunk = ceiling(stranded_position/100)) %>%
  summarise(normalized_depth = mean(normalized_depth)) %>% 
  mutate(stranded_position = new_chunk*100) %>% 
  mutate(filename = "mean")



new_df_with_mean <- new_average_depth_100 %>% 
  bind_rows(new_mean_depth)

linear_ref <- new_df_with_mean %>% 
  filter(filename != "CD95") %>%
  # filter(filename == "CD73") %>% 
  ggplot(aes(x = stranded_position, y = normalized_depth, color = filename)) +
  geom_line(data = . %>% filter(filename != "mean"), aes(group = filename), color = "lightgrey") +
  geom_line(data = . %>% filter(filename == "mean"), aes(group = filename), color = "black", size = 1.5) +
  theme_bw() +
  xlab("Genome Position") + 
  ylab("Fold change in coverage per contig") +
  scale_color_manual(values = c("mean" = "black")) +
  guides(size = FALSE) 
# facet_wrap(~filename)# This removes the size legend
# scale_color_identity(guide = "none")

linear_ref

#find the max mean coverage for the left vs. right origin(s)

# Filter the mean group from new_df_with_mean
mean_data <- new_df_with_mean %>%
  filter(filename == "mean")

# Find the maximum mean value for each range
max_mean_0_12500 <- mean_data %>%
  filter(stranded_position >= 0 & stranded_position <= 12500) %>%
  summarise(max_mean = max(normalized_depth, na.rm = TRUE))

max_mean_75000_100000 <- mean_data %>%
  filter(stranded_position >= 75000 & stranded_position <= 100000) %>%
  summarise(max_mean = max(normalized_depth, na.rm = TRUE))

# Display results
max_mean_0_12500
max_mean_75000_100000

