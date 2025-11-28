############################################################
# Script: compile_ground_coffee_HSI_txt_to_csv.R
# Purpose:
#   - Read multiple HSI-derived .txt files (one per sample)
#   - Extract wavelength-specific mean reflectance
#   - Compile all samples into a single matrix (rows = samples,
#     columns = wavelengths)
#   - Export as a .csv file for further chemometric analysis
############################################################

## 1. Set working directory (update this path as needed)
setwd("C:/Users/abc/OneDrive - UGent/Documenten/Derick Malavi_PhD Docs_UGent/Manuscript 4_Coffee/Compiled_ROI_Ground")

## 2. List all .txt files exported from the HSI software
list_of_files <- list.files(
  path       = ".",
  pattern    = "\\.txt$",
  recursive  = TRUE,
  full.names = TRUE
)

## 3. Derive sample names from file names
#    - Remove the last 5 characters (e.g. "_1.txt")
#    - Then drop the first 2 characters if needed (e.g. "./")
SampleNames <- vector("character", length(list_of_files))
for (i in seq_along(list_of_files)) {
  name_i <- gsub(".{5}$", "", list_of_files[i])   # strip last 5 chars
  name_i <- sub(".{2}", "", name_i)               # strip first 2 chars
  SampleNames[i] <- name_i
}

## 4. Read the first file to:
##    - Identify wavelength column
##    - Extract the "mean" reflectance column
first_sample <- read.delim(list_of_files[1], header = TRUE, sep = "")

# Remove header/meta rows (first 6 rows) and non-spectral columns (7:9)
first_sample <- first_sample[-(1:6), -(7:9)]

# Rename columns for clarity
colnames(first_sample) <- c("wavelength", "min", "mean_minus_SD", "mean", "mean_plus_SD", "max")

# Store wavelength axis (assumed identical for all files)
wavelengths <- first_sample[["wavelength"]]

# Extract mean reflectance for the first sample
rowsMatrix <- vector("list", length(list_of_files))
rowsMatrix[[1]] <- first_sample[["mean"]]

## 5. Loop over remaining files and extract mean reflectance
for (i in 2:length(list_of_files)) {
  x <- read.delim(list_of_files[i], header = TRUE, sep = "")
  x <- x[-(1:6), -(7:9)]          # drop metadata and non-spectral cols
  rowsMatrix[[i]] <- x[, 4]       # 4th column = "mean"
}

## 6. Combine all samples into a single data frame
#    - Each list element is a vector of reflectance values
DATA <- data.frame(rowsMatrix[[1]])
for (i in 2:length(rowsMatrix)) {
  DATA <- cbind(DATA, as.numeric(rowsMatrix[[i]]))
}

## 7. Convert to matrix with the desired orientation:
##    - Rows   = samples
##    - Cols   = wavelengths
MAT <- t(as.matrix(DATA))
colnames(MAT) <- wavelengths
rownames(MAT) <- SampleNames

## 8. Export compiled matrix as .csv
write.csv(MAT, file = "ground_coffee_compiled.csv", row.names = TRUE)
#------------------------------------------------------------------------------