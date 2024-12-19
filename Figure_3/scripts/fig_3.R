'''
Plot GC skew for pcrass genome from Skewit data
Date: Mar 18 2024
'''

library("tidyverse")
library("ggplot2")
library("gggenes")
library(readr)
library(stringr)
library(dplyr)
library(forcats)
library(cowplot)


#read in the output from Skewit run on the command line
gc_data <- read_tsv("~/Downloads/pcrass_assembly_gcskew.txt", 
                    col_names = c("Sequence", "Index", "reverse_index", "GC_skew"), skip=1)

# Adjusting the initial data processing steps to use reverse_index
tmp <- gc_data %>% 
  add_row(Index=-1, reverse_index=-1, GC_skew=0) %>% 
  arrange(reverse_index) %>% 
  mutate(GC_skew2=c(GC_skew[-1], 0)) %>% 
  mutate(reverse_index2=c(reverse_index[-1], 0)) %>% 
  mutate(sign=as.factor(case_when(GC_skew <= 0 ~ -1, GC_skew >= 0 ~ 1))) %>% 
  filter(row_number() <= n()-1)

# Update the gc_plot to use reverse_index instead of Index
gc_plot <- gc_data %>% 
  mutate(sign = as.factor(case_when(GC_skew <= 0 ~ -1, GC_skew >= 0 ~ 1))) %>% 
  ggplot(aes(x=reverse_index, y=0)) + 
  geom_rect(aes(xmin = 16400, xmax = 22400, ymin = -1, ymax = 1), 
            fill = 'lightgray', alpha = 0.1) +
  geom_segment(aes(col=sign, yend=GC_skew, xend=reverse_index), alpha=0.3) +
  geom_segment(data=tmp, aes(y=GC_skew, yend=GC_skew2, x=reverse_index, xend=reverse_index2, col=sign)) +
  theme_bw() +
  xlab("Genome position (reverse index)") + 
  ylab("GC Skew") +
  scale_color_manual(values = c("1" = "#009292FF", "-1" = "black"), guide = "none") +
  theme(panel.grid.minor = element_blank()) + 
  geom_vline(xintercept = max(gc_data$reverse_index), linetype = "dotted") +
  geom_vline(xintercept = min(gc_data$reverse_index), linetype = "dotted")

gc_plot



# Updated zoom_plot to zoom in on the x-axis range from 12000 to 15000
zoom_plot <- gc_data %>% 
  mutate(sign = as.factor(case_when(GC_skew <= 0 ~ -1, GC_skew >= 0 ~ 1))) %>% 
  ggplot(aes(x=reverse_index, y=0)) + 
  geom_segment(aes(col=sign, yend=GC_skew, xend=reverse_index), alpha=0.3) +
  geom_segment(data=tmp, aes(y=GC_skew, yend=GC_skew2, x=reverse_index, xend=reverse_index2, col=sign)) +
  theme_bw() +
  xlab("Genome position (reverse index)") + 
  ylab("GC Skew") +
  scale_color_manual(values = c("1" = "#009292FF", "-1" = "black"), guide = "none") +
  theme(panel.grid.minor = element_blank()) +
  xlim(16400, 22400)

zoom_plot


#import gff3 file for pcrass reference genome 

# importing the gff3 file for the pcrasss genome 
pcrass_annotation <- read_tsv("~/Downloads/crass_B1_assembly_dna_ligase_linear_phynteny_annotations.gff",
                              col_names = c("genome","origin", "type", "start", "end", "dot", 
                                            "strand", "idk", "info"), skip = 5)
# Define the new feature details
new_feature_1 <- tibble(
  genome = "genome",      # Adjust if necessary
  origin = "origin",      # Adjust if necessary
  type = "origin_1",      # Type of the new feature
  start = 74713,         # Start position
  end = 74829,           # End position
  dot = ".",             # Dot field, typically "."
  strand = "+",          # Strand field, typically "."
  idk = ".",             # Optional field, typically "."
  info = "origin_1"      # Additional information
)

new_feature_2 <- tibble(
  genome = "genome",      # Adjust if necessary
  origin = "origin",      # Adjust if necessary
  type = "origin_2",      # Type of the new feature
  start = 75370,         # Start position
  end = 75378,           # End position
  dot = ".",             # Dot field, typically "."
  strand = "+",          # Strand field, typically "."
  idk = ".",             # Optional field, typically "."
  info = "origin_2"      # Additional information
)


new_feature_3 <- tibble(
  genome = "genome",      # Adjust if necessary
  origin = "origin",      # Adjust if necessary
  type = "origin_3",      # Type of the new feature
  start = 75520,         # Start position
  end = 75606,           # End position
  dot = ".",             # Dot field, typically "."
  strand = "+",          # Strand field, typically "."
  idk = ".",             # Optional field, typically "."
  info = "origin_3"      # Additional information
)


new_feature_4 <- tibble(
  genome = "genome",      # Adjust if necessary
  origin = "origin",      # Adjust if necessary
  type = "origin_4",      # Type of the new feature
  start = 75733,         # Start position
  end = 76035,           # End position
  dot = ".",             # Dot field, typically "."
  strand = "+",          # Strand field, typically "."
  idk = ".",             # Optional field, typically "."
  info = "origin_4"      # Additional information
)


new_feature_5 <- tibble(
  genome = "genome",      # Adjust if necessary
  origin = "origin",      # Adjust if necessary
  type = "origin_5",      # Type of the new feature
  start = 78637,         # Start position
  end = 79449,           # End position
  dot = ".",             # Dot field, typically "."
  strand = "+",          # Strand field, typically "."
  idk = ".",             # Optional field, typically "."
  info = "origin_5"      # Additional information
)

# Append the new feature to the existing data
updated_annotation <- bind_rows(pcrass_annotation, new_feature_1, new_feature_2, new_feature_3, new_feature_4, new_feature_5)

# Adjust start and end for reverse genes to simulate flipping
updated_annotation <- updated_annotation %>%
  mutate(xmin = ifelse(strand == "-", end, start),
         xmax = ifelse(strand == "-", start, end)) %>% 
  mutate(f=strand=='+')


# Plot gene zoom in
gene_plot <- ggplot(updated_annotation) +
  geom_gene_arrow(aes(xmin = start, xmax = end, forward=f, 
                      y = "C. communis", fill = strand), 
                  arrowhead_height = unit(10, "mm"),  # Shorter arrowhead for a less pronounced tip
                  arrowhead_width = unit(5, "mm"),   # Same width as the body
                  arrow_body_height = unit(10, "mm")) + # Keep the body height wide
  theme_genes() +
  scale_x_reverse(limits=c(80000, 74000)) + 
  scale_fill_manual(values = c("+" = "#009292FF", "-" = "grey32")) + # Custom colors for + and - strands
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none",
        axis.text.y = element_blank()) +
  NULL

gene_plot




# Filter out rows with "origin_" in the info column
updated_annotation_filtered <- updated_annotation %>%
  filter(!grepl("^origin_", info))  # Exclude rows where 'info' starts with 'origin_'

# Mutate to add a new column with reversed strand information
updated_annotation <- updated_annotation %>%
  mutate(
    strand_reverse = ifelse(strand == "+", "-", "+")  # Reverse the strand direction
  )

#plot all genes
gene_plot_2 <- ggplot(updated_annotation) +
  geom_gene_arrow(aes(xmin = start, xmax = end, forward = f,
                      y = "C. communis", fill = strand_reverse), 
                  arrowhead_height = unit(10, "mm"),  # Shorter arrowhead for a less pronounced tip
                  arrowhead_width = unit(5, "mm"),   # Same width as the body
                  arrow_body_height = unit(10, "mm")) + # Keep the body height wide
  theme_genes() +
  scale_fill_manual(values = c("-" = "#009292FF", "+" = "grey32")) + # Custom colors for + and - strands
  scale_x_reverse() +  # Reverse the x-axis to invert the plot
  # theme(axis.title.y = element_blank(), 
        # axis.title.x = element_blank(), 
        # axis.ticks.x = element_blank(), 
        # axis.line.x = element_blank(), 
        # axis.text.x = element_blank(), 
        # legend.position = "none",
        # axis.text.y = element_blank()) + 
  NULL

gene_plot_2

#this combines the gene plot and the zoomed in plot and aligns them 
comb_zoom_plots <- plot_grid(zoom_plot, gene_plot, ncol = 1, align = "v", axis='lr')

#this combines the full gene plot and the full GC skew plot and aligns them 
comb_plots <- plot_grid(gene_plot_2, gc_plot, ncol = 1, align = "v", axis='lr')

# thsi plots both the overall GC skew plot and the combined plot from above 
g <- plot_grid(comb_plots, comb_zoom_plots, nrow = 2)
g

#save as a pdf with a transparent background 
ggsave(g, filename = "~/Documents/241105_fig3.pdf", width = 10, height = 7, useDingbats = FALSE, bg = "transparent")
