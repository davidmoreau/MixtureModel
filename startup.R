### SETUP ###

# WD
setwd("/Users/davidmoreau/Google Drive/Papers/In Progress/MixtureModel") #change path to yours

# Colors
palette <- c("#58D9BF", "#46ABE3", "#4485C5", "#3276BE", "#2A5766") #(from lighter to darker)

# Load packages (+ install if required)
packages <- c("psych", "ggplot2", "gridExtra", "grid", "compute.es", "metafor", "pwr", "scales", "mixtools", "dplyr")
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
sapply(packages, suppressPackageStartupMessages(require), warn.conflicts=F, quietly=T, character.only=T)

# Set default ggplot parameters
default <-  theme_bw() + 
  #theme(plot.title=element_text(family="Arial", face="bold", size=16)) +
  theme(axis.line.x = element_line(colour = "black", size=.5, linetype = 1),
        axis.line.y = element_line(colour = "black", size=.5, linetype = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 