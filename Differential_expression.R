# Load libraries
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


### MATRIX 
# Filling the array files data in datos21 folder
data_store <- lapply(paste0('datos21/', sample$File.Name), function(file) {
  read.table(file, sep = '\t', header = FALSE)[, 2]
})

# Combine the data_store list into a matrix for DEGs
count_matrix <- do.call(cbind, data_store)
# Convert into integers
count_matrix <- round(count_matrix)

# Get genes ensemble IDs without .version (all 60483)
df_genes <- data.frame(
  survey = read.table('datos21/5a067a75-9a3d-4b31-812f-c041620867d7.FPKM.txt',
                      sep='\t', header=FALSE)[, 1]
)
df_genes$survey <- sapply(strsplit(df_genes$survey, "\\."), function(x) x[1])

# Create the df_expression data frame by combining df_genes and expression_matrix
df_expression <- data.frame(count_matrix, stringsAsFactors = FALSE)
rownames(df_expression) <- t(df_genes)


# Perform DESGs Analisys 
######################################
# Construct a DESqDataSet object-------
dds <- DESeqDataSetFromMatrix(countData = df_expression, 
                       colData = sample, 
                       design = ~ Sample.Type)

# Pre-filtering: Removing rows with 0 gene counts
keep <- rowSums(counts(dds)) > 0
dds <- dds[keep,]

# Set the factor level
dds$Sample.Type <- relevel(dds$Sample.Type, ref = "Primary Tumor")

# Run DESeq -----
dds <- DESeq(dds)
res <- results(dds)

res

# Export the matrix to CSV
write.csv(res, file = "DEGs.csv", row.names = FALSE)

####### THE END ##########



