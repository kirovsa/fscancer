#!/usr/bin/env Rscript
# Script to extract gene attributes from Ensembl using biomaRt
# 
# Purpose:
#   Extracts the following attributes for all human genes from Ensembl BioMart:
#   - Stable Ensembl gene ID
#   - Transcript stable ID
#   - Peptide stable ID
#   - Gene name
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
#   ensembl_gene_id|ensembl_transcript_id|ensembl_peptide_id|gene_name|cds_length|transcript_biotype
#
# Example output:
#   ENSG00000000003|ENST00000373020|ENSP00000362111|TSPAN6|879|protein_coding
#   ENSG00000000005|ENST00000373031|ENSP00000362122|TNMD|1455|protein_coding
#
# Note:
#   This script requires internet access to connect to the Ensembl BioMart database.
#   The query may take several minutes to complete depending on network speed.

# Load required library
library(biomaRt)

# Function to extract gene attributes
extract_gene_attributes <- function(output_file = "gene_attributes.txt") {
  
  message("Connecting to Ensembl BioMart database...")
  
  # Connect to Ensembl database (human genome)
  tryCatch({
    ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  }, error = function(e) {
    stop("Failed to connect to Ensembl BioMart database.\n",
         "Error: ", conditionMessage(e), "\n",
         "Please check your internet connection and try again.")
  })
  
  message("Fetching gene attributes...")
  
  # Define the attributes to retrieve
  attributes <- c(
    "ensembl_gene_id",          # Stable Ensembl gene ID
    "ensembl_transcript_id",    # Transcript stable ID
    "ensembl_peptide_id",       # Peptide stable ID
    "external_gene_name",       # Gene name
    "cds_length",               # CDS length
    "transcript_biotype"        # Transcript type
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
  message(paste("  Unique genes:", length(unique(gene_data$ensembl_gene_id))))
  message(paste("  Unique transcripts:", length(unique(gene_data$ensembl_transcript_id))))
  message(paste("  Records with peptide IDs:", sum(!is.na(gene_data$ensembl_peptide_id) & gene_data$ensembl_peptide_id != "")))
  message(paste("  Records with CDS length:", sum(!is.na(gene_data$cds_length) & gene_data$cds_length != "")))
  
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
