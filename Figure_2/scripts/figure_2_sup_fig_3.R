#scripts to generate figure 2 plots
library(ggplot2)
library(readxl)
library(dplyr)
library(patchwork)
###figure 2A##############
# Read the Excel file
#data <- read_excel("~/Documents/230619_final_metatranscriptomics_output.xlsx")
data <- read_excel("~/Downloads/241219_final_metatranscriptomics_output.xlsx")

# Create subsets for different x ranges (genes on the forward vs reverse strand)
subset1 <- subset(data, start_position_of_gene >= 0 & start_position_of_gene <= 29300)
subset2 <- subset(data, start_position_of_gene > 29300 & start_position_of_gene <= 95900)

# Assign group labels
subset1$group <- "plasmid"
subset2$group <- "phage"

# Combine the subsets into a single data frame
combined_data <- rbind(subset1, subset2)

# Aggregate data to get average RPKM for each start position
average_data <- aggregate(RPKM ~ start_position_of_gene + group, combined_data, mean)

# Reorder the 'group' factor to ensure proper ordering on the x-axis
average_data$group <- factor(average_data$group, levels = c("phage", "plasmid"))

# Perform an unpaired t-test between the two groups
t_test_result <- t.test(RPKM ~ group, data = combined_data)
t_test_label <- paste("p =", signif(t_test_result$p.value, 3))

# Create the boxplot fig 5A
average_plot <- ggplot(average_data, aes(x = group, y = RPKM, fill = group, color = group)) +
  geom_boxplot(alpha = 0.5, outlier.shape=NA) +  # Adjusted alpha for better visualization
  geom_jitter(width = 0.2, size = 2, alpha = 0.5) +  # Add jittered scatterplot
  labs(title = "A",
       x = "Group",
       y = "Average RPKM") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(trans = "log10") +
  theme_bw() +  # Apply the black and white theme
  
  # Add t-test for means
  annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 1,
           label = t_test_label, color = "black", parse = FALSE) +  
  
  # Customize colors for both points and boxes
  scale_color_manual(values = c("black", "#009292FF"), labels = c("phage", "plasmid")) +
  scale_fill_manual(values = c("black", "#009292FF"), labels = c("phage", "plasmid")) +
  # Remove the legend
  theme(legend.position = "none")



# Display the plot
print(average_plot)

# Save the plot as a PDF file
ggsave("~/Documents/average_box_plot.pdf", average_plot, width = 4, height = 4, units = "in", dpi = 300)


###for figure 2B################
# Read the Excel file
data2 <- read_excel("~/Downloads/241219_pos_neg_express.xlsx")

# Create a black and white theme
bw_theme <- theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "grey"),
    panel.grid.minor = element_blank(),
    text = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.text = element_text(color = "black")
  )

# Sort the data along the x-axis based on ascending y-axis values
data2_sorted <- data2[order(data2$`pos to neg ratio`), ]

# Create the scatter plot
p2 <- ggplot(data2_sorted, aes(x = reorder(`sample_number`, `pos to neg ratio`), y = `pos to neg ratio`, fill = `pos to neg ratio` > 1)) +
  geom_bar(stat = "identity") +  
  labs(title = "B",
       x = "Sample Number",
       y = "Average RPKM plasmid:phage") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey59", size = 1) +
  scale_fill_manual(values = c("black", "#009292FF"), labels = c("Below 1", "Above 1")) +  
  theme_bw() +
  theme(
    text = element_text(size = 20),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.line = element_line(linetype = 1, linewidth = 1),
    legend.position = "none") +  
  scale_y_continuous(
    trans = "log10",
    breaks = c(0.7, 1, 10),
    limits = c(0.7, 10)  # Define the y-axis range explicitly
  ) +  
  scale_x_discrete(labels = function(x) ifelse(as.integer(x) %% 10 == 1, x, ""),
                   breaks = function(x) x[as.integer(x) %% 10 == 1])  

p2


ggsave(path = "~/Documents", width = 10, height = 6, device='pdf', dpi=1000, filename = "metatran_pos_neg_RPKM_num.pdf")


### Combine Plots
combined_plot <- average_plot + p2 + plot_layout(ncol = 2, widths = c(1, 2))

# Display the combined plot
print(combined_plot)

# Save the combined plot as a PDF file
ggsave("241219_figure2_combined.pdf", combined_plot, width = 12, height = 6, units = "in", dpi = 300)

#####for supplementary to figure 2##########
# Read your Excel file (replace "your_file.xlsx" with the actual file path)
data <- readxl::read_excel("~/Downloads/241219_final_metatranscriptomics_output.xlsx", sheet= "RPKM")

# Filter genes with RPKM values above zero
filtered_data <- data %>%
  filter(RPKM > 0)

# Count the number of samples where each gene has RPKM above zero
gene_count <- filtered_data %>%
  group_by(start_position_of_gene) %>%
  summarise(Count = n_distinct(sample))

gene_count$start_position_of_gene <- factor(gene_count$start_position_of_gene, levels = unique(gene_count$start_position_of_gene))

# Create a bar plot
ggplot(gene_count, aes(x = start_position_of_gene, y = Count)) +
  geom_bar(stat = "identity", fill = "#009292FF", color = "black") +
  geom_hline(yintercept = 39, linetype = "dashed", color = "grey", size = 1) +  
  labs(title = "Count of Samples with RPKM > 0 for Each Gene",
       x = "Gene",
       y = "Sample Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better visibility


# Calculate the total RPKM for each gene across all samples
gene_total_rpkm <- data %>%
  group_by(start_position_of_gene) %>%
  summarise(Total_RPKM = sum(RPKM))

# Make the x-axis categorical
gene_total_rpkm$start_position_of_gene <- factor(gene_total_rpkm$start_position_of_gene, levels = unique(gene_total_rpkm$start_position_of_gene))

# Create a bar plot with text labels
ggplot(gene_total_rpkm, aes(x = start_position_of_gene, y = Total_RPKM)) +
  geom_bar(stat = "identity", fill = "#009292FF", color = "black") +
  geom_text(aes(label = round(Total_RPKM, 2)), vjust = -0.5, size = 3) +  # Add text labels above each bar
  labs(title = "Total RPKM for Each Gene Across All Samples",
       x = "Gene",
       y = "Total RPKM") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better visibility



# Calculate the average RPKM for each gene across all samples
gene_avg_rpkm <- data %>%
  group_by(start_position_of_gene) %>%
  summarise(Average_RPKM = mean(RPKM))

# Make the x-axis categorical
gene_avg_rpkm$start_position_of_gene <- factor(gene_avg_rpkm$start_position_of_gene, levels = unique(gene_avg_rpkm$start_position_of_gene))

# Create a bar plot with text labels
sup_express <- ggplot(gene_avg_rpkm, aes(x = start_position_of_gene, y = Average_RPKM)) +
  geom_bar(stat = "identity", fill = "#009292FF", color = "black") +
  labs(title = "Average RPKM for Each Gene Across All Samples",
       x = "Gene",
       y = "Average RPKM") +
  scale_y_log10() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better visibility
sup_express

ggsave("~/Documents/sup_expression.pdf", sup_express, width = 12, height = 6, units = "in", dpi = 300)
