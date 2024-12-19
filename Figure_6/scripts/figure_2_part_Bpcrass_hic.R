# Load the required library
library(ggplot2)
library(readxl)
library(dplyr)
library(here)

#read in data
hicdata <- 
  read_excel("~/Documents/Figures_post_reviews/scripts/figure_2/figure_2B_pcrass_hic_plot_data.xlsx")

# Fit a linear regression model (line that fits the data)
model <- lm(links ~ abundance, data = hicdata)

# Calculate residuals (differences between the observed values and the predicted values for each data point)
residuals <- residuals(model)

# Calculate the standard error of residuals
residual_se <- sd(residuals)

# Create a data frame with species, residuals, 
#and z-scores (Z-score of +2 indicates that the data point is 
#2 standard deviations above the mean, while a Z-score of 
#-1.5 indicates that it is 1.5 standard deviations below the mean.)
outlier_data <- data.frame(
  species = hicdata$species,
  links = hicdata$links,
  abundance = hicdata$abundance,
  residuals = residuals,
  z_score = residuals / residual_se
)
# Define a threshold for identifying points with more links than expected by chance
threshold <- 3  # standard threshold for looking for outliers is 3

# Create a plot
pcrass_hic_host_predict <- ggplot(data = outlier_data, aes(x = abundance, y = links)) + 
  geom_point(size = 2, shape = 1, color = "black") + 
  geom_smooth(data = outlier_data, method = lm, se = FALSE, fullrange = TRUE, 
              formula = y ~ 0 + x, linetype = "dashed", 
              color = "black", size = 0.5) + 
  labs(title = "Hi-C links between p-crAss and \n bacteria in the stool sample", 
       x = "Bacterial MAG abundance \n (average read depth/genome length)", 
       y = "number of links between \n bacterial MAG and p-crAss") + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position = "none") +
  geom_text(data = outlier_data, 
            aes(x = abundance, y = links, label = species, color = (abs(z_score) > threshold)), 
            hjust = -0.1, vjust = 0.5, 
            fontface = "italic", size = 3) + 
  #labs(tag = "B") +
  scale_color_manual(values = c("azure4", "#009292FF")) +  # Set label color to black by default, but change color if z score passes threshold
  xlim(c(0, 52)) +
  ylim(c(0,8000))

pcrass_hic_host_predict

outlier_data

#convert z score for B. vulgatus to p-value (4.693218e-11)
pnorm(6.47651650, lower.tail=FALSE)

#convert z score for Prevotella species to p-value (0.7815714)
pnorm(-0.77751134, lower.tail=FALSE)

#convert z score for Blautia species to p-value (0.1050754)
pnorm(1.25315108, lower.tail=FALSE)

#convert z score for Bacteroides stercoris to p-value (0.340737)
pnorm(0.41045247, lower.tail=FALSE)

#convert z score for Ruminococcus species to p-value (0.2782499)
pnorm(0.58804842, lower.tail=FALSE)

#convert z score for Acidaminococcus species to p-value (0.7283435)
pnorm(-0.60781073, lower.tail=FALSE)

#convert z score for Parabacteroides merdae to p-value (0.5188798)
pnorm(-0.04734240, lower.tail=FALSE)

ggsave(here("~/Documents/figure2B.pdf"), dpi=300, w=7, h=6)

