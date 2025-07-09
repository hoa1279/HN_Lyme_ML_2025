# ==============================================================================
# PROTEIN SEQUENCE ONE-HOT ENCODING PIPELINE
# ==============================================================================
# This script processes protein sequences from multiple proteins and converts 
# them to one-hot encoded format for machine learning applications.

# ==============================================================================
# REQUIRED LIBRARIES
# ==============================================================================
library(caret)          # For dummy variable creation
library(Peptides)       # For amino acid utilities
library(Biostrings)     # For reading FASTA files
library(readxl)         # For reading Excel files

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

#' Get list of amino acids plus gap character
#' @return Character vector of 21 amino acids plus "-" gap character
aaListN <- function() {
  c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
    "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "-")
}

#' Extract One-Hot Encoded features from protein sequences
#' 
#' This function converts protein sequences into one-hot encoded format where
#' each amino acid at each position is represented as a binary feature.
#' 
#' @param df Data frame containing sequences or vector of sequences
#' @param sequence_col Name of column containing sequences (default: "Sequence")
#' @param docking_col Optional name of column containing docking scores
#' @param n Optional integer to generate all possible sequences of length n
#' @return Data frame with one-hot encoded sequences
extract_OHE <- function(df = NULL, sequence_col = "Sequence", docking_col = NULL, n = NULL) {
  
  # Check required packages
  if (!requireNamespace("caret", quietly = TRUE)) {
    stop("Package 'caret' is needed for this function to work. Please install it.")
  }
  if (!requireNamespace("Peptides", quietly = TRUE)) {
    stop("Package 'Peptides' is needed for this function to work. Please install it.")
  }
  
  # Input validation and data preparation
  if (!is.null(df)) {
    if (is.vector(df)) {
      df <- data.frame(Sequence = df)
    } else if (!is.data.frame(df)) {
      stop("'df' must be a data frame or a vector.")
    }
  } else if (!is.null(n)) {
    # Generate all possible sequences of length n
    if (!is.numeric(n) || n <= 2 || (n %% 1 != 0)) {
      stop("n must be a positive integer greater than 2.")
    }
    
    # Create all possible peptide combinations
    PeList <- expand.grid(rep(list(aaListN()), n))
    sequences <- do.call(paste0, PeList)
    df <- data.frame(Sequence = sequences)
  } else {
    stop("Either 'df' or 'n' must be provided.")
  }
  
  # Validate column names
  col_names <- names(df)
  
  if (!(sequence_col %in% col_names)) {
    stop(paste("Column", sequence_col, "not found in the data frame."))
  }
  
  if (!is.null(docking_col) && !(docking_col %in% col_names)) {
    stop(paste("Column", docking_col, "not found in the data frame."))
  }
  
  # Sequence processing
  sequences <- toupper(df[[sequence_col]])
  max_length <- max(nchar(sequences))
  
  # Pad shorter sequences with "Z" to ensure uniform length
  if (any(nchar(sequences) < max_length)) {
    warning("Not all sequences are of the same length. Padding shorter sequences with 'Z'.")
  }
  
  padded_sequences <- sapply(sequences, function(sequence) {
    if (nchar(sequence) < max_length) {
      padding_needed <- max_length - nchar(sequence)
      return(paste(c(strsplit(sequence, "")[[1]], rep("Z", padding_needed)), collapse = ""))
    } else {
      return(sequence)
    }
  })
  
  # Convert sequences to matrix format (each column = position, each row = sequence)
  split_sequences <- strsplit(padded_sequences, "")
  sequence_df <- do.call("rbind", lapply(split_sequences, function(sequence) {
    data.frame(t(sequence))
  }))
  
  # Remove columns with only one unique value (no variation)
  check_single_levels <- function(df) {
    single_levels <- sapply(df, function(x) length(unique(x)) < 2)
    return(names(df)[single_levels])
  }
  
  single_level_cols <- check_single_levels(sequence_df)
  if (length(single_level_cols) > 0) {
    warning(paste("Removing columns with only one level:", paste(single_level_cols, collapse = ", ")))
    sequence_df <- sequence_df[, !(names(sequence_df) %in% single_level_cols)]
  }
  
  # Create one-hot encoded matrix
  dummies <- caret::dummyVars(~ ., data = sequence_df)
  one_hot_encoded_sequences <- data.frame(stats::predict(dummies, newdata = sequence_df))
  
  # Ensure consistent column structure across all datasets
  all_residues <- c(Peptides::aaList(), if (any(nchar(sequences) < max_length)) "Z")
  template_names <- unlist(sapply(1:max_length, function(x) paste0("X", x, all_residues)))
  
  # Add missing columns filled with zeros
  missing_cols <- setdiff(template_names, names(one_hot_encoded_sequences))
  if (length(missing_cols) > 0) {
    one_hot_encoded_sequences[, missing_cols] <- 0
  }
  
  # Organize final output columns
  one_hot_encoded_sequences[[sequence_col]] <- df[[sequence_col]]
  new_col_order <- c(sequence_col, setdiff(names(one_hot_encoded_sequences), sequence_col))
  
  if (!is.null(docking_col)) {
    one_hot_encoded_sequences[[docking_col]] <- df[[docking_col]]
    new_col_order <- c(sequence_col, docking_col, 
                       setdiff(names(one_hot_encoded_sequences), c(sequence_col, docking_col)))
  }
  
  one_hot_encoded_sequences <- one_hot_encoded_sequences[, new_col_order]
  
  return(one_hot_encoded_sequences)
}

# ==============================================================================
# MAIN PROCESSING PIPELINE
# ==============================================================================

# Get current working directory
cwd <- getwd()

# Define list of proteins to process
protein_list <- c("OspC", "OspA", "DbpA", "BBK32", "RevA", "P66", "BB_0406")


# Process each protein individually
for (i in seq_along(protein_list)) {
  
  prot <- protein_list[i]
  cat("Processing protein", i, "of", length(protein_list), ":", prot, "\n")

  # Step 1: Read protein sequences from FASTA file
  data_path <- file.path(cwd, "data", prot)
  setwd(data_path)
  fasta_file <- paste0(prot, "_cut_unique_seqs.fasta")
  cat("  Reading sequences from:", fasta_file, "\n")
  
  if (!file.exists(fasta_file)) {
    warning(paste("FASTA file not found:", fasta_file, "- Skipping", prot))
    next
  }
  
  sequences <- readAAStringSet(fasta_file)
  cat("  Found", length(sequences), "sequences\n")
  
  # Convert to data frame format
  sequence_df <- data.frame(
    name = names(sequences),
    sequence = as.character(sequences)
  )
  
  # Step 2: Read metadata from Excel file
  excel_file <- paste0(prot, "_unq_seqs.xlsx")
  cat("  Reading metadata from:", excel_file, "\n")
  
  if (!file.exists(excel_file)) {
    warning(paste("Excel file not found:", excel_file, "- Skipping", prot))
    next
  }
  
  metadata <- readxl::read_excel(excel_file, sheet = "use")
  metadata <- metadata[-1]  # Remove first column (assuming it's an index)
  
  # Step 3: Merge sequences with metadata
  cat("  Merging sequences with metadata...\n")
  combined_data <- merge(metadata, sequence_df, by.x = "seq_name", by.y = "name")
  
  # Step 4: Filter out sequences with missing dissemination data
  na_indices <- which(combined_data$Dissemination == "NA")
  if (length(na_indices) > 0) {
    cat("  Removing", length(na_indices), "sequences with missing dissemination data\n")
    combined_data <- combined_data[-na_indices, ]
  }
  
  cat("  Final dataset contains", nrow(combined_data), "sequences\n")
  
  # Step 5: Generate one-hot encoded features
  cat("  Generating one-hot encoded features...\n")
  
  # Note: The original code uses "x" as sequence column name, but the merged data
  # likely has "sequence" as the column name. Adjust as needed.
  sequence_column <- "sequence"  # Change this if your column name is different
  
  if (!sequence_column %in% names(combined_data)) {
    cat("  Available columns:", paste(names(combined_data), collapse = ", "), "\n")
    warning(paste("Sequence column", sequence_column, "not found in", prot, "data"))
    next
  }
  
  ohe_data <- extract_OHE(combined_data, sequence_col = sequence_column)
  rownames(ohe_data) <- combined_data$seq_name
  
  # Step 6: Save results
  output_file <- paste0(prot, "_OHE.csv")
  cat("  Saving results to:", output_file, "\n")
  
  # Remove the sequence column before saving (keep only numeric features)
  write.csv(ohe_data[, -1], output_file, row.names = TRUE)
  
  cat("  Successfully processed", prot, "\n\n")
}

cat("Pipeline completed!\n")

