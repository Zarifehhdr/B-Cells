#####################################################################################################
# Preprocessing gene expressiong matrix ready for SLIDE analysis
#####################################################################################################

# Define the input and output folders
input_path <- ".../Data/Raw_Data/Male_Donor/"
output_folder <- ".../Data/Preprocessed/Male_Donor/"
dir.create(output_folder, recursive = TRUE) 

# Set the zero-filtering threshold
zero_filtering_threshold <- 0.3  # Adjust this value as needed

# Load the datasets
# you can define specific names instead of GC and PB
dataset_paths <- list(
  GC = paste0(input_path, "GC.csv"),
  PB = paste0(input_path, "PB.csv"),
  MBC = paste0(input_path, "MBC.csv"),
  ABC_d2 = paste0(input_path, "ABC_d2.csv"),
  ABC_d4 = paste0(input_path, "ABC_d4.csv")
)

# Function to load datasets
load_datasets <- function(paths) {
  datasets <- lapply(paths, function(path) {
    read.csv(path, row.names = 1)
  })
  names(datasets) <- names(paths)
  return(datasets)
}

# Function to preprocess each dataset
preprocess_data <- function(data, zero_filtering_threshold) {
  # Step 1: Remove mitochondrial genes
  mt_cols1 <- colnames(data)[startsWith(colnames(data), "MT.")]
  mt_cols2 <- colnames(data)[startsWith(colnames(data), "MT-")]
  
  data <- data[, !colnames(data) %in% mt_cols1]
  data <- data[, !colnames(data) %in% mt_cols2]
  
  # Step 2: Remove ribosomal genes
  mt_cols_rp <- grep("^RP[LS]", colnames(data), value = TRUE)
  data <- data[, !colnames(data) %in% mt_cols_rp]
  
  # Step 3: Zero Filtering
  filtered_data <- ZeroFiltering(data, 0, zero_filtering_threshold * dim(data)[1])
  
  return(filtered_data)
}

# Function to filter data for using SLIDE
ZeroFiltering <- function(data, g_thresh, c_thresh) {
  cat("Original dataframe dimension is", dim(data)[[1]], "by", dim(data)[[2]], "\n")
  # Filter out cells with num zeros greater than a certain threshold
  thresh <- dim(data)[[2]] - g_thresh
  i <- rowSums(data == 0, na.rm=TRUE) <= thresh
  filtered <- data[i, ]
  
  # Filter out genes with num zeros greater than a certain threshold
  thresh <- dim(data)[[1]] - c_thresh
  i <- colSums(data == 0, na.rm=TRUE) <= thresh
  filtered <- filtered[, i]
  
  cat("Filtered dataframe dimension is", nrow(filtered), "by", ncol(filtered), "\n")
  return(filtered)
}

# Load datasets
datasets <- load_datasets(dataset_paths)

# Preprocess all datasets and save the filtered versions
for (name in names(datasets)) {
  filtered_data <- preprocess_data(datasets[[name]], zero_filtering_threshold)
  write.csv(filtered_data, file = paste0(output_folder, name, ".csv"), row.names = TRUE, col.names = TRUE)
}
