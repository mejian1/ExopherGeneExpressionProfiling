# Load required libraries
library(data.table)

# Load microarray data from file
data <- fread("microarray_data.csv")

# Background correction using median polish
bg_corrected_data <- data - rowMedians(data) - colMedians(data) + median(data)

# Quantile normalization
sorted_data <- sort(bg_corrected_data)
rank_data <- apply(bg_corrected_data, 2, rank)
norm_data <- rowMeans(sorted_data[rank_data, ])

# Summarization by gene-level expression values
gene_expression <- rowMeans(matrix(norm_data, ncol = 2, byrow = TRUE))

# Create a data frame with gene names and expression values
processed_data <- data.frame(Gene = colnames(data), Expression = gene_expression)

# Save processed data to file
write.csv(processed_data, file = "processed_microarray_data.csv", row.names = FALSE)
