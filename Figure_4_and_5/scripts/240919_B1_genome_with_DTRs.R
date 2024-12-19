# Load the necessary libraries

library(ggplot2)
library(gggenes)
library(scales)
library(stringr)

# Load the GenBank (.gbk) file

gbk_file <- "~/Downloads/phynteny_B1_DTR_final_genome.gbk"

# Read the GenBank file into R as text

genbank_lines <- readLines(gbk_file)

# Initialize empty vectors
locus_tags <- c()
start_positions <- c()
end_positions <- c()
strands <- c()
functions <- c()

# Flag to track feature blocks
in_feature_block <- FALSE

# Loop through the lines and extract gene features
for (i in 1:length(genbank_lines)) {
  line <- genbank_lines[i]
  
  # Identify the start of a CDS (coding sequence) block
  if (grepl("^\\s+CDS", line)) {
    in_feature_block <- TRUE
    
    # Extract positions and strand information
    if (grepl("complement", line)) {
      strand <- "-"
      positions <- gsub(".*complement\\((\\d+)\\.\\.(\\d+)\\).*", "\\1,\\2", line)
    } else {
      strand <- "+"
      positions <- gsub(".*CDS\\s+(\\d+)\\.\\.(\\d+).*", "\\1,\\2", line)
    }
    
    start_positions <- c(start_positions, as.numeric(strsplit(positions, ",")[[1]][1]))
    end_positions <- c(end_positions, as.numeric(strsplit(positions, ",")[[1]][2]))
    strands <- c(strands, strand)
  }
  
  # Extract locus_tag and product fields
  if (in_feature_block) {
    if (grepl("/locus_tag=", line)) {
      locus_tag <- gsub(".*/locus_tag=\"(.*)\"", "\\1", line)
      locus_tags <- c(locus_tags, locus_tag)
    }
    
    if (grepl("/product=", line)) {
      product <- gsub(".*/product=\"(.*)\"", "\\1", line)
      functions <- c(functions, product)
      in_feature_block <- FALSE  # End of feature block
    }
  }
  
  # Handle cases where product field is missing
  if (in_feature_block && grepl("/CDS", line) && length(locus_tags) < length(start_positions)) {
    locus_tags <- c(locus_tags, NA)
    functions <- c(functions, "Unknown function")
    in_feature_block <- FALSE
  }
}


# Ensure all vectors have the same length
if (length(locus_tags) < length(start_positions)) {
  locus_tags <- c(locus_tags, rep(NA, length(start_positions) - length(locus_tags)))
}
if (length(functions) < length(start_positions)) {
  functions <- c(functions, rep("Unknown function", length(start_positions) - length(functions)))
}

# Create the dataframe
gene_data <- data.frame(
  locus_tag = locus_tags,
  start = start_positions,
  end = end_positions,
  strand = strands,
  gene_function = functions,
  stringsAsFactors = FALSE
)
# Define a palette of colors
colors <- c(
  "#0033FF",
  "#0072B2",   
  "#00BEFF", 
  "#56B4E9",
  "#0099FF",
  "#9FBDFF",
  "#00FFFF",
  "#00E5BF",
  "#006600",
  "#5A9D47",
  "#95D840",
  "#94FFB7",
  "#33FF33",
  "#009999",
  "#9C97C1",
  "#7E4F9E",
  "#440154",
  "#CC33FF",
  "#9C5580",
  "#CC79A7",
  "#F564E3",
  "#FFB5BE",
  "#E6347D",
  "#E63B2E",
  "#FF6600",
  "#FF9900",
  "#FFC300",
  "#F0E442",
  "#FFFE91",
  "#F0F921",
  "#DA7372",
  "#993300",
  "#B58B76", 
  "#C2E0D6",
  "#BCCCDC",
  "#999999",
  "#333333"
)


# Create a named vector for colors and assign 'hypothetical protein' to white
color_map <- setNames(colors, unique(gene_data$gene_function))
color_map["hypothetical protein"] <- "white"

# Define positions for the boxes
first_box <- data.frame(xmin = 1, xmax = 1379, ymin = 0.4, ymax = 1.4)  # Adjusted heights
last_box <- data.frame(xmin = 97882 - 1379, xmax = 97882, ymin = 0.4, ymax = 1.4)  # Adjusted heights

# Plot using the pastel color palette
plot <- ggplot(gene_data, aes(xmin = start, xmax = end, y = 1, fill = gene_function)) +
  geom_gene_arrow(aes(forward = strand == "+")) +
  scale_fill_manual(values = color_map) +  # Apply the pastel color palette
  theme_genes() +
  labs(title = "Genome Map from GenBank File", x = "Genomic Position", y = "") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  scale_y_continuous(breaks = NULL) +
  geom_rect(data = first_box, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
            fill = NA, color = "red", linewidth = 1) +  # First box
  geom_rect(data = last_box, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
            fill = NA, color = "red", linewidth = 1)  # Last box

# Print the plot
print(plot)
ggsave("B1_genome_with_DTRs.pdf", plot, width = 20, height = 8)


# Save the color mapping to a CSV file for Python use
color_mapping_df <- data.frame(gene_function = names(color_map), color = color_map)
write.csv(color_mapping_df, "gene_color_mapping.csv", row.names = FALSE)

