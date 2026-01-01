#!/usr/bin/env Rscript
# Script to annotate Ensembl transcripts with biotype information
# 
# Purpose:
#   Takes a file containing Ensembl transcript IDs or gene IDs and annotates
#   them with biotype information from the Ensembl BioMart database.
#
# Requirements:
#   - R (>= 3.5.0)
#   - biomaRt package (available via BiocManager or apt: r-bioc-biomart)
#
# Installation:
#   Option 1 - Using BiocManager (requires internet access):
#     if (!requireNamespace("BiocManager", quietly = TRUE))
#         install.packages("BiocManager")
#     BiocManager::install("biomaRt")
#   
#   Option 2 - Using apt (Debian/Ubuntu):
#     sudo apt-get install r-bioc-biomart
#
# Usage:
#   Rscript annotate_transcripts.R <input_file> <output_file> [id_type]
#
# Arguments:
#   input_file  - Path to input file with IDs (one per line or tab/comma-separated)
#   output_file - Path to output file with annotations
#   id_type     - Optional. Type of ID in input: "transcript" (default), "gene", or "auto"
#
# Input Format:
#   The input file should contain Ensembl IDs, one per line or in the first column
#   of a tab-separated or comma-separated file. Supported ID formats:
#   - Transcript IDs: ENST00000000000
#   - Gene IDs: ENSG00000000000
#
# Output Format:
#   Tab-separated file with original data plus annotation columns:
#   - ensembl_gene_id
#   - ensembl_transcript_id
#   - transcript_biotype
#   - external_gene_name (gene symbol)
#
# Examples:
#   # Annotate transcript IDs
#   Rscript annotate_transcripts.R transcripts.txt annotated_transcripts.txt
#   
#   # Annotate gene IDs
#   Rscript annotate_transcripts.R genes.txt annotated_genes.txt gene
#   
#   # Auto-detect ID type
#   Rscript annotate_transcripts.R ids.txt annotated_ids.txt auto
#
# Note:
#   This script requires internet access to connect to the Ensembl BioMart database.

# Load required library
library(biomaRt)

# Configuration constants
BATCH_SIZE <- 500  # Number of IDs to query per batch
MAX_BIOTYPE_DISPLAY <- 10  # Maximum number of biotypes to display in summary

# Function to detect ID type based on prefix
detect_id_type <- function(ids) {
  # Remove NA and empty values
  ids <- ids[!is.na(ids) & ids != ""]
  
  if (length(ids) == 0) {
    stop("No valid IDs found in input file")
  }
  
  # Count ENST and ENSG prefixes
  enst_count <- sum(grepl("^ENST", ids))
  ensg_count <- sum(grepl("^ENSG", ids))
  
  if (enst_count > ensg_count) {
    return("transcript")
  } else if (ensg_count > enst_count) {
    return("gene")
  } else {
    stop("Cannot auto-detect ID type. Please specify 'transcript' or 'gene' as id_type parameter.")
  }
}

# Function to read input file
read_input_file <- function(input_file) {
  # Try to read as simple list first (one ID per line)
  tryCatch({
    data <- read.table(input_file, header = FALSE, stringsAsFactors = FALSE, 
                      sep = "\t", quote = "", comment.char = "")
    return(data)
  }, error = function(e) {
    # Try comma-separated
    tryCatch({
      data <- read.table(input_file, header = FALSE, stringsAsFactors = FALSE, 
                        sep = ",", quote = "", comment.char = "")
      return(data)
    }, error = function(e2) {
      stop("Failed to read input file. Please ensure it's a valid tab or comma-separated file.\n",
           "Error: ", conditionMessage(e))
    })
  })
}

# Main annotation function
annotate_transcripts <- function(input_file, output_file, id_type = "auto") {
  
  message("Reading input file: ", input_file)
  
  # Read input file
  input_data <- read_input_file(input_file)
  
  # Extract IDs from first column
  ids <- as.character(input_data[, 1])
  
  # Auto-detect ID type if needed
  if (id_type == "auto") {
    id_type <- detect_id_type(ids)
    message("Auto-detected ID type: ", id_type)
  }
  
  # Validate ID type
  if (!id_type %in% c("transcript", "gene")) {
    stop("Invalid id_type. Must be 'transcript', 'gene', or 'auto'")
  }
  
  message("Connecting to Ensembl BioMart database...")
  
  # Connect to Ensembl database (human genome)
  tryCatch({
    ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  }, error = function(e) {
    stop("Failed to connect to Ensembl BioMart database.\n",
         "Error: ", conditionMessage(e), "\n",
         "Please check your internet connection and try again.")
  })
  
  message("Querying BioMart for annotations...")
  
  # Define attributes and filters based on ID type
  attributes <- c(
    "ensembl_gene_id",
    "ensembl_transcript_id",
    "transcript_biotype",
    "external_gene_name"
  )
  
  if (id_type == "transcript") {
    filter_name <- "ensembl_transcript_id"
  } else {
    filter_name <- "ensembl_gene_id"
  }
  
  # Query BioMart in batches to handle large input
  all_annotations <- NULL
  
  # Split IDs into batches
  num_batches <- ceiling(length(ids) / BATCH_SIZE)
  
  for (i in 1:num_batches) {
    start_idx <- (i - 1) * BATCH_SIZE + 1
    end_idx <- min(i * BATCH_SIZE, length(ids))
    batch_ids <- ids[start_idx:end_idx]
    
    message(sprintf("Processing batch %d/%d (%d IDs)...", i, num_batches, length(batch_ids)))
    
    tryCatch({
      batch_annotations <- getBM(
        attributes = attributes,
        filters = filter_name,
        values = batch_ids,
        mart = ensembl
      )
      
      if (is.null(all_annotations)) {
        all_annotations <- batch_annotations
      } else {
        all_annotations <- rbind(all_annotations, batch_annotations)
      }
    }, error = function(e) {
      warning("Failed to retrieve data for batch ", i, ": ", conditionMessage(e))
    })
  }
  
  if (is.null(all_annotations) || nrow(all_annotations) == 0) {
    stop("No annotations retrieved from BioMart. Please check your IDs and try again.")
  }
  
  message(paste("Retrieved annotations for", nrow(all_annotations), "records"))
  
  # Merge annotations with input data
  # Create a data frame with the input IDs
  if (id_type == "transcript") {
    input_df <- data.frame(ensembl_transcript_id = ids, stringsAsFactors = FALSE)
    merge_by <- "ensembl_transcript_id"
  } else {
    input_df <- data.frame(ensembl_gene_id = ids, stringsAsFactors = FALSE)
    merge_by <- "ensembl_gene_id"
  }
  
  # Add original data if input has multiple columns
  if (ncol(input_data) > 1) {
    for (j in 2:ncol(input_data)) {
      input_df[[paste0("V", j)]] <- input_data[, j]
    }
  }
  
  # Merge with annotations
  annotated_data <- merge(input_df, all_annotations, by = merge_by, all.x = TRUE)
  
  # Write to output file
  tryCatch({
    write.table(
      annotated_data,
      file = output_file,
      sep = "\t",
      row.names = FALSE,
      col.names = TRUE,
      quote = FALSE,
      na = ""
    )
  }, error = function(e) {
    stop("Failed to write output file.\n",
         "Error: ", conditionMessage(e))
  })
  
  message(paste("Annotated data saved to:", output_file))
  
  # Print summary statistics
  message("\nSummary:")
  message(paste("  Input IDs:", length(ids)))
  message(paste("  IDs with annotations:", sum(!is.na(annotated_data$transcript_biotype))))
  message(paste("  IDs without annotations:", sum(is.na(annotated_data$transcript_biotype))))
  
  # Show biotype distribution
  if (sum(!is.na(annotated_data$transcript_biotype)) > 0) {
    message("\nBiotype distribution:")
    biotype_counts <- table(annotated_data$transcript_biotype, useNA = "no")
    print(head(sort(biotype_counts, decreasing = TRUE), MAX_BIOTYPE_DISPLAY))
  }
  
  return(annotated_data)
}

# Main execution
if (!interactive()) {
  # Parse command line arguments
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) < 2) {
    cat("Usage: Rscript annotate_transcripts.R <input_file> <output_file> [id_type]\n")
    cat("\n")
    cat("Arguments:\n")
    cat("  input_file  - Path to input file with Ensembl IDs\n")
    cat("  output_file - Path to output file for annotated data\n")
    cat("  id_type     - Optional. Type of ID: 'transcript' (default), 'gene', or 'auto'\n")
    cat("\n")
    cat("Examples:\n")
    cat("  Rscript annotate_transcripts.R transcripts.txt annotated.txt\n")
    cat("  Rscript annotate_transcripts.R genes.txt annotated.txt gene\n")
    cat("  Rscript annotate_transcripts.R ids.txt annotated.txt auto\n")
    quit(status = 1)
  }
  
  input_file <- args[1]
  output_file <- args[2]
  id_type <- ifelse(length(args) >= 3, args[3], "auto")
  
  # Check if input file exists
  if (!file.exists(input_file)) {
    stop("Input file does not exist: ", input_file)
  }
  
  # Run the annotation
  result <- annotate_transcripts(input_file, output_file, id_type)
  
  message("\nDone!")
}
