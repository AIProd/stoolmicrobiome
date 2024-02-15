#library(phyloseq)
library(MicrobiotaProcess) # an R package for analysis, visualization and biomarker discovery of Microbiome.

library(phyloseq) # Handling and analysis of high-throughput microbiome census data.
library(ggplot2) # Create Elegant Data Visualisations Using the Grammar of Graphics.
library(tidyverse) # Easily Install and Load the 'Tidyverse'.
library(vegan) # Community Ecology Package.
library(coin) # Conditional Inference Procedures in a Permutation Test Framework.
library(reshape2) # Flexibly Reshape Data: A Reboot of the Reshape Package.
library(ggnewscale) # Multiple Fill and Colour Scales in 'ggplot2'.
library(phyloseq)
library(dplyr)
library(UpSetR)
library(microbiome)
library(microViz)
library(microbiome)
library(dplyr)

setwd('D:/base/cystfastq/analysis')
otu <- read.csv('tax_cystfastq.csv',row.names = 1,check.names = FALSE)
meta <- read.csv('meta_cystfastq.csv',row.names = 1)
taxa <- read.csv('linm_cystfastq.csv',row.names = 1)

OTU = otu_table(as.matrix(otu), taxa_are_rows = TRUE)
TAX = phyloseq::tax_table(as.matrix(taxa))
samples = sample_data(meta)
ps=phyloseq(OTU, TAX, samples)
classtaxa <- get_taxadf(obj=ps, taxlevel=6)


relab_genera = transform_sample_counts(classtaxa, function(x) x / sum(x)) 
relab_genera = transform_sample_counts(classtaxa, function(x) x / sum(x) * 100) 
rel_abundance_df <- as.data.frame(t(otu_table(relab_genera)))
rel_abundance_df$sample_names <- sample_names(relab_genera)
rel_abundance_df$NAC <- relab_genera@sam_data$NAC
write.csv(rel_abundance_df, file = "relative_abundance_cystfastq_genus1.csv", row.names = FALSE)
