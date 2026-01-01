#!/usr/bin/env Rscript
# Script to extract gene attributes from Ensembl using ensembldb
# 
# Purpose:
#   Extracts the following attributes for all human genes from a local Ensembl database:
#   - Transcript stable ID
#   - Stable Ensembl gene ID
#   - Gene name
#   - Peptide stable ID
#   - CDS length
#   - Transcript type
#
# Requirements:
#   - R (>= 3.5.0)
#   - ensembldb package (available via BiocManager)
#   - An EnsDb annotation package (e.g., EnsDb.Hsapiens.v86, EnsDb.Hsapiens.v109)
#
# Installation:
#   # Install ensembldb and an annotation package
#   if (!requireNamespace("BiocManager", quietly = TRUE))
#       install.packages("BiocManager")
#   BiocManager::install("ensembldb")
#   
#   # Install a specific human genome annotation (choose one based on your needs):
#   # For Ensembl release 86 (GRCh38):
#   BiocManager::install("EnsDb.Hsapiens.v86")
#   
#   # For Ensembl release 109 (GRCh38):
#   BiocManager::install("EnsDb.Hsapiens.v109")
#   
#   # Or for the latest version, install from AnnotationHub:
#   BiocManager::install("AnnotationHub")
#
# Usage:
#   # Using a specific EnsDb package:
#   Rscript extract_gene_attributes_ensembldb.R [output_file] [ensdb_package]
#
#   # Using AnnotationHub (auto-downloads latest):
#   Rscript extract_gene_attributes_ensembldb.R [output_file] hub
#
# Arguments:
#   output_file    - Optional. Path to output file (default: gene_attributes_ensembldb.txt)
#   ensdb_package  - Optional. EnsDb package name (default: EnsDb.Hsapiens.v86)
#                    Use "hub" to auto-download latest from AnnotationHub
#
# Examples:
#   # Use default EnsDb.Hsapiens.v86:
#   Rscript extract_gene_attributes_ensembldb.R
#
#   # Use a specific version:
#   Rscript extract_gene_attributes_ensembldb.R output.txt EnsDb.Hsapiens.v109
#
#   # Use AnnotationHub for latest version:
#   Rscript extract_gene_attributes_ensembldb.R output.txt hub
#
# Output Format:
#   Pipe-separated (|) file with no header, containing:
#   ensembl_transcript_id|ensembl_gene_id|gene_name|ensembl_peptide_id|cds_length|transcript_biotype
#
# Example output:
#   ENST00000373020|ENSG00000000003|TSPAN6|ENSP00000362111|879|protein_coding
#   ENST00000373031|ENSG00000000005|TNMD|ENSP00000362122|1455|protein_coding
#
# Advantages over biomaRt:
#   - Works offline (no internet connection required after package installation)
#   - Faster queries (local database)
#   - More reproducible (fixed version)
#   - No dependency on remote server availability
#
# Note:
#   The output format is identical to extract_gene_attributes.R and matches
#   the expected input for frameshift.Rmd.

# Load required libraries
suppressPackageStartupMessages({
  library(ensembldb)
})

# Function to get EnsDb from AnnotationHub
get_ensdb_from_hub <- function() {
  message("Loading AnnotationHub to get latest EnsDb for human...")
  
  # Check if AnnotationHub is available
  if (!requireNamespace("AnnotationHub", quietly = TRUE)) {
    stop("AnnotationHub package is required when using 'hub' option.\n",
         "Please install it with: BiocManager::install('AnnotationHub')")
  }
  
  library(AnnotationHub)
  ah <- AnnotationHub()
  
  # Query for human EnsDb annotations
  query_result <- query(ah, c("EnsDb", "Homo sapiens"))
  
  if (length(query_result) == 0) {
    stop("No EnsDb databases found for Homo sapiens in AnnotationHub")
  }
  
  # Get the most recent version
  # EnsDb entries are typically named like "Ensembl 109 EnsDb for Homo sapiens"
  latest_idx <- length(query_result)
  message(paste("Using:", query_result[latest_idx]$title))
  message(paste("Version:", query_result[latest_idx]$description))
  
  edb <- query_result[[latest_idx]]
  return(edb)
}

# Function to load EnsDb package
load_ensdb <- function(ensdb_name) {
  if (ensdb_name == "hub") {
    return(get_ensdb_from_hub())
  }
  
  # Try to load the specified package
  if (!requireNamespace(ensdb_name, quietly = TRUE)) {
    stop(paste0("Package '", ensdb_name, "' is not installed.\n",
                "Please install it with: BiocManager::install('", ensdb_name, "')\n",
                "Or use 'hub' to auto-download from AnnotationHub."))
  }
  
  # Load the package and get the EnsDb object
  message(paste("Loading", ensdb_name, "..."))
  library(ensdb_name, character.only = TRUE)
  
  # Get the EnsDb object (usually has the same name as the package)
  edb <- get(ensdb_name)
  
  return(edb)
}

# Function to extract gene attributes using ensembldb
extract_gene_attributes_ensembldb <- function(output_file = "gene_attributes_ensembldb.txt",
                                               ensdb_name = "EnsDb.Hsapiens.v86") {
  
  # Load the EnsDb database
  edb <- load_ensdb(ensdb_name)
  
  message(paste("Database version:", ensembldb::ensemblVersion(edb)))
  message(paste("Genome build:", ensembldb::genome(edb)))
  message(paste("Organism:", ensembldb::organism(edb)))
  
  message("\nFetching gene attributes...")
  
  # Use the transcripts() function to get all transcripts with associated information
  # This is more efficient than multiple queries
  tryCatch({
    # Get all transcripts with their associated gene and protein information
    tx_data <- transcripts(
      edb,
      columns = c(
        "tx_id",              # Transcript stable ID
        "gene_id",            # Gene stable ID
        "gene_name",          # Gene symbol
        "protein_id",         # Protein stable ID
        "tx_cds_seq_start",   # CDS start (to calculate length)
        "tx_cds_seq_end",     # CDS end (to calculate length)
        "tx_biotype"          # Transcript biotype
      ),
      return.type = "DataFrame"
    )
  }, error = function(e) {
    stop("Failed to retrieve data from EnsDb.\n",
         "Error: ", conditionMessage(e))
  })
  
  # Convert to regular data frame for easier manipulation
  tx_df <- as.data.frame(tx_data)
  
  # Calculate CDS length from start and end positions
  # If either is NA, the CDS length will be NA
  tx_df$cds_length <- ifelse(
    !is.na(tx_df$tx_cds_seq_start) & !is.na(tx_df$tx_cds_seq_end),
    tx_df$tx_cds_seq_end - tx_df$tx_cds_seq_start + 1,
    NA
  )
  
  # Select and rename columns to match the expected output format
  # Order: tx_id, gene_id, gene_name, protein_id, cds_length, tx_biotype
  gene_data <- data.frame(
    ensembl_transcript_id = tx_df$tx_id,
    ensembl_gene_id = tx_df$gene_id,
    external_gene_name = tx_df$gene_name,
    ensembl_peptide_id = tx_df$protein_id,
    cds_length = tx_df$cds_length,
    transcript_biotype = tx_df$tx_biotype,
    stringsAsFactors = FALSE
  )
  
  # Remove any rows where transcript ID is missing (shouldn't happen, but just in case)
  gene_data <- gene_data[!is.na(gene_data$ensembl_transcript_id), ]
  
  if (nrow(gene_data) == 0) {
    stop("No data retrieved from EnsDb. Please check your database installation.")
  }
  
  message(paste("Retrieved", nrow(gene_data), "records"))
  
  # Calculate summary statistics
  unique_genes <- length(unique(gene_data$ensembl_gene_id))
  unique_transcripts <- length(unique(gene_data$ensembl_transcript_id))
  records_with_peptide <- sum(!is.na(gene_data$ensembl_peptide_id) & nzchar(as.character(gene_data$ensembl_peptide_id)))
  records_with_cds <- sum(!is.na(gene_data$cds_length))
  
  # Write to file with pipe separator (matching the format used in frameshift.Rmd)
  tryCatch({
    write.table(
      gene_data,
      file = output_file,
      sep = "|",
      row.names = FALSE,
      col.names = FALSE,
      quote = FALSE,
      na = ""
    )
  }, error = function(e) {
    stop("Failed to write output file.\n",
         "Error: ", conditionMessage(e))
  })
  
  message(paste("Gene attributes saved to:", output_file))
  
  # Print summary statistics
  message("\nSummary:")
  message(paste("  Total records:", nrow(gene_data)))
  message(paste("  Unique genes:", unique_genes))
  message(paste("  Unique transcripts:", unique_transcripts))
  message(paste("  Records with peptide IDs:", records_with_peptide))
  message(paste("  Records with CDS length:", records_with_cds))
  
  # Show transcript type distribution
  message("\nTranscript type distribution:")
  type_counts <- table(gene_data$transcript_biotype)
  print(head(sort(type_counts, decreasing = TRUE), 10))
  
  return(gene_data)
}

# Main execution
if (!interactive()) {
  # Parse command line arguments
  args <- commandArgs(trailingOnly = TRUE)
  
  output_file <- ifelse(length(args) > 0, args[1], "gene_attributes_ensembldb.txt")
  ensdb_name <- ifelse(length(args) > 1, args[2], "EnsDb.Hsapiens.v86")
  
  # Run the extraction
  result <- extract_gene_attributes_ensembldb(output_file, ensdb_name)
  
  message("\nDone!")
}
