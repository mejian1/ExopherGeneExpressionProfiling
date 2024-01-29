# ImportOurData
# import dataset > from text (base) > click "RawMdt15RNAiData.csv" > headings click "yes" > row names choose "first column" > import
View(RawMdt15RNAiData)
rawdata <- as.matrix(RawMdt15RNAiData)

#data <- fread("RawMdt15RNAiData.csv")
#df <- as.data.frame(RawMdt15RNAiData)

#drops <- c("ID_REF")
#dfdata <- df[ , !(names(df) %in% drops)]

#RGlista <- new("RGList", df)

#MA <- normalizeWithinArrays(df)
#MA <- MA.RG(RGlista, bc.method="subtract", offset=50)

# Background correction using median polish
bg_corrected_data <- backgroundCorrect.matrix(rawdata, Eb = NULL, method = "auto", offset = 0, printer = NULL, normexp.method = "saddle", verbose = TRUE)
#bg_corrected_data <- (rawdata) - (rowMedians(rawdata)) #+ (median(rawdata)) #- (colMedians(rawdata))

dotchart(rawdata)
dotchart(bg_corrected_data)

# Quantile normalization
#sorted_data <- sort(bg_corrected_data)
#rank_data <- apply(bg_corrected_data, 2, rank)
#norm_data <- rowMeans(sorted_data[rank_data, ])

norm_data <- normalizeQuantiles(bg_corrected_data, ties = TRUE)
qqplot(norm_data, bg_corrected_data) # line = normalized = good
qqplot(rawdata, bg_corrected_data) # not a line = not normalized

# Summarization by gene-level expression values
#gene_expression <- rowMeans(norm_data, ncol = 2, byrow = TRUE)

# Create a data frame with gene names and expression values
#processed_data <- data.frame(Gene = colnames(rawdata), Expression = gene_expression)

# Save processed data to file
write.csv(norm_data, file = "processed_microarray_data.csv", row.names = TRUE)
