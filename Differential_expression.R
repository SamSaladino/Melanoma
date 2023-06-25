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
# Convert into intergers
count_matrix <- round(count_matrix)

# Get genes ensemble IDs without .version (all 60483)
df_genes <- data.frame(
  survey = read.table('datos21/5a067a75-9a3d-4b31-812f-c041620867d7.FPKM.txt',
                      sep='\t', header=FALSE)[, 1]
)
df_genes$survey <- sapply(strsplit(df_genes$survey, "\\."), function(x) x[1])

# Create the df_expression data frame by combining df_genes and expression_matrix
df_expression <- data.frame(df_genes,expression_matrix, stringsAsFactors = FALSE)


########## DESGs Analisys ####################
# Construct a DESqDataSet object-------
dds <- DESeqDataSetFromMatrix(countData = count_matrix, 
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

# Function to check if a value is an integer
isNonInteger <- function(x) {
  !is.integer(x) || any(x != as.integer(x))
}

# Function to find non-integer values in a matrix
findNonIntegerValues <- function(mat) {
  nonIntegerValues <- mat[apply(mat, c(1, 2), isNonInteger)]
  unique(na.omit(nonIntegerValues))
}


# Check for non-integer values
nonIntegerValues <- findNonIntegerValues(count_matrix)

# Print non-integer values if any
if (length(nonIntegerValues) > 0) {
  cat("The matrix contains non-integer values:\n")
  print(nonIntegerValues)
} else {
  cat("The matrix contains only integers.\n")
 }




# Check if all elements are integers
if (is.integer(count_matrix) && all(count_matrix == as.integer(count_matrix))) {
  print("The matrix contains only integers.")
} else {
  print("The matrix contains non-integer values.")
}
