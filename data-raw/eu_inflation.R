library(readxl)
library(dplyr)
library(tibble)

# Read raw file
raw <- read_excel(
  'data-raw/data_inflation_EU.xlsx',
  skip = 8
)

# The first column is the category name
colnames(raw)[1] <- "category"

# Remove category column for now (will store separately)
categories <- raw$category

# Remove first and last unwanted columns
raw <- raw[, -1]
raw <- raw[, -ncol(raw)]

# Convert to matrix
eu_inflation <- as.matrix(raw)

# Standardize each row
eu_inflation <- t(apply(eu_inflation, 1, scale))

# Add row names
rownames(eu_inflation) <- categories

# Save as internal dataset
usethis::use_data(eu_inflation, overwrite = TRUE)
