#!/usr/bin/env Rscript
# Script to extract gene attributes from Ensembl using biomaRt
# 
# Purpose:
#   Extracts the following attributes for all human genes from Ensembl BioMart:
#   - Transcript stable ID
#   - Stable Ensembl gene ID
#   - Gene name
#   - Peptide stable ID
#   - CDS length
#   - Transcript type
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
#   Rscript extract_gene_attributes.R [output_file]
#
# Arguments:
#   output_file - Optional. Path to output file (default: gene_attributes.txt)
#
# Output Format:
#   Pipe-separated (|) file with no header, containing:
#   ensembl_transcript_id|ensembl_gene_id|gene_name|ensembl_peptide_id|cds_length|transcript_biotype
#
# Example output:
#   ENST00000373020|ENSG00000000003|TSPAN6|ENSP00000362111|879|protein_coding
#   ENST00000373031|ENSG00000000005|TNMD|ENSP00000362122|1455|protein_coding
#
# Note:
#   This script requires internet access to connect to the Ensembl BioMart database.
#   The query may take several minutes to complete depending on network speed.
#   The output format matches the expected input for frameshift.Rmd.
#   If the main Ensembl server times out, the script will automatically try the
#   US East mirror (useast.ensembl.org) as a fallback.

# Load required library
library(biomaRt)

# Function to extract gene attributes
extract_gene_attributes <- function(output_file = "gene_attributes.txt") {
  
  message("Connecting to Ensembl BioMart database...")
  
  # Connect to Ensembl database (human genome)
  # Try main server first, then fall back to mirror if timeout occurs
  ensembl <- NULL
  
  # List of hosts to try in order (main, then US East mirror)
  hosts <- c("https://www.ensembl.org", "https://useast.ensembl.org")
  host_names <- c("main Ensembl server", "US East mirror")
  
  # Pattern for detecting connection errors
  connection_error_pattern <- "timeout|timed out|cannot open|failed to connect"
  
  for (i in seq_along(hosts)) {
    tryCatch({
      message(sprintf("Attempting to connect to %s (%s)...", host_names[i], hosts[i]))
      ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = hosts[i])
      message(sprintf("Successfully connected to %s", host_names[i]))
      break  # Success, exit the loop
    }, error = function(e) {
      error_msg <- conditionMessage(e)
      # Check if it's a timeout or connection error
      if (grepl(connection_error_pattern, error_msg, ignore.case = TRUE)) {
        message(sprintf("Connection to %s failed with timeout or connection error.", host_names[i]))
        if (i < length(hosts)) {
          message("Trying alternative mirror...")
        }
      } else {
        # For non-timeout errors, still try the next host
        message(sprintf("Connection to %s failed: %s", host_names[i], error_msg))
        if (i < length(hosts)) {
          message("Trying alternative mirror...")
        }
      }
    })
  }
  
  # If all connection attempts failed
  if (is.null(ensembl)) {
    stop("Failed to connect to any Ensembl BioMart server (tried main and mirror sites).\n",
         "Please check your internet connection and try again later.\n",
         "The Ensembl servers may be temporarily unavailable.")
  }
  
  message("Fetching gene attributes...")
  
  # Define the attributes to retrieve
  # Order matches the expected format in frameshift.Rmd: ENST, ENSG, Symbol, ENSP, ProtLen, Type
  attributes <- c(
    "ensembl_transcript_id",    # Transcript stable ID (ENST)
    "ensembl_gene_id",          # Stable Ensembl gene ID (ENSG)
    "external_gene_name",       # Gene name (Symbol)
    "ensembl_peptide_id",       # Peptide stable ID (ENSP)
    "cds_length",               # CDS length (ProtLen)
    "transcript_biotype"        # Transcript type (Type)
  )
  
  # Query BioMart
  tryCatch({
    gene_data <- getBM(
      attributes = attributes,
      mart = ensembl
    )
  }, error = function(e) {
    stop("Failed to retrieve data from BioMart.\n",
         "Error: ", conditionMessage(e), "\n",
         "The Ensembl server may be temporarily unavailable. Please try again later.")
  })
  
  if (nrow(gene_data) == 0) {
    stop("No data retrieved from BioMart. Please check your connection and try again.")
  }
  
  message(paste("Retrieved", nrow(gene_data), "records"))
  
  # Calculate summary statistics once for reuse
  unique_genes <- length(unique(gene_data$ensembl_gene_id))
  unique_transcripts <- length(unique(gene_data$ensembl_transcript_id))
  records_with_peptide <- sum(!is.na(gene_data$ensembl_peptide_id) & gene_data$ensembl_peptide_id != "")
  records_with_cds <- sum(!is.na(gene_data$cds_length) & gene_data$cds_length != "")
  
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
  
  output_file <- ifelse(length(args) > 0, args[1], "gene_attributes.txt")
  
  # Run the extraction
  result <- extract_gene_attributes(output_file)
  
  message("\nDone!")
}
