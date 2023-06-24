#load libraries
##########################
library(DESeq2)
library(TCGAbiolinks)
#########################

#Preparing count data
########################
setwd("/home/sam/Documents/Melanoma")

##### LOCAL DATA #############
#Read samples information
info <- read.csv("sample.csv")

# Remove unnecessary samples
sample <- info[!(info$`Sample.Type` == 'Solid Tissue Normal' | 
                   info$`Sample.Type` == 'Additional Metastatic'), ]
# Reset index
sample <- sample[order(row.names(sample)),]
rownames(sample) <- NULL  







