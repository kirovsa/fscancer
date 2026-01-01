#!/usr/bin/env Rscript
# Example demonstrating the annotate_transcripts.R script functionality
# This creates sample output showing what the script produces when run with internet access

# Example 1: Annotating transcript IDs
cat("Example 1: Annotating Ensembl Transcript IDs\n")
cat("=" , rep("=", 50), "\n", sep = "")

# Sample input data (transcript IDs)
transcript_ids <- c(
  "ENST00000373020",
  "ENST00000373031",
  "ENST00000371582",
  "ENST00000367770",
  "ENST00000286031"
)

# Sample annotated output
transcript_annotations <- data.frame(
  ensembl_transcript_id = transcript_ids,
  ensembl_gene_id = c(
    "ENSG00000000003",
    "ENSG00000000005",
    "ENSG00000000419",
    "ENSG00000000457",
    "ENSG00000000460"
  ),
  transcript_biotype = c(
    "protein_coding",
    "protein_coding",
    "protein_coding",
    "protein_coding",
    "protein_coding"
  ),
  external_gene_name = c(
    "TSPAN6",
    "TNMD",
    "DPM1",
    "SCYL3",
    "C1orf112"
  ),
  stringsAsFactors = FALSE
)

# Create example input file
input_file_transcripts <- "example_transcript_ids.txt"
writeLines(transcript_ids, input_file_transcripts)

# Create example output file
output_file_transcripts <- "example_transcript_annotations.txt"
write.table(
  transcript_annotations,
  file = output_file_transcripts,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE,
  na = ""
)

cat("\nInput file created:", input_file_transcripts, "\n")
cat("Output file created:", output_file_transcripts, "\n")
cat("\nSample output (first 3 lines):\n")
cat(paste(readLines(output_file_transcripts, n = 4), collapse = "\n"), "\n\n")

# Example 2: Annotating gene IDs
cat("\nExample 2: Annotating Ensembl Gene IDs\n")
cat("=" , rep("=", 50), "\n", sep = "")

# Sample input data (gene IDs)
gene_ids <- c(
  "ENSG00000000003",
  "ENSG00000000005",
  "ENSG00000000419"
)

# Sample annotated output (note: genes can have multiple transcripts)
gene_annotations <- data.frame(
  ensembl_gene_id = c(
    "ENSG00000000003",
    "ENSG00000000003",
    "ENSG00000000005",
    "ENSG00000000419",
    "ENSG00000000419"
  ),
  ensembl_transcript_id = c(
    "ENST00000373020",
    "ENST00000494424",
    "ENST00000373031",
    "ENST00000371582",
    "ENST00000371588"
  ),
  transcript_biotype = c(
    "protein_coding",
    "processed_transcript",
    "protein_coding",
    "protein_coding",
    "protein_coding"
  ),
  external_gene_name = c(
    "TSPAN6",
    "TSPAN6",
    "TNMD",
    "DPM1",
    "DPM1"
  ),
  stringsAsFactors = FALSE
)

# Create example input file
input_file_genes <- "example_gene_ids.txt"
writeLines(gene_ids, input_file_genes)

# Create example output file
output_file_genes <- "example_gene_annotations.txt"
write.table(
  gene_annotations,
  file = output_file_genes,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE,
  na = ""
)

cat("\nInput file created:", input_file_genes, "\n")
cat("Output file created:", output_file_genes, "\n")
cat("\nSample output (first 4 lines):\n")
cat(paste(readLines(output_file_genes, n = 5), collapse = "\n"), "\n\n")

# Summary
cat("\n")
cat("Usage Summary\n")
cat("=" , rep("=", 50), "\n", sep = "")
cat("\nThese examples demonstrate the format of the output from annotate_transcripts.R\n\n")
cat("Key Features:\n")
cat("  - Annotates Ensembl transcript IDs or gene IDs with biotype information\n")
cat("  - Queries the Ensembl BioMart database\n")
cat("  - Returns tab-separated output with annotation columns\n")
cat("  - Handles multiple transcripts per gene when annotating gene IDs\n\n")

cat("Command Examples:\n")
cat("  # Annotate transcript IDs (auto-detect):\n")
cat("  Rscript annotate_transcripts.R", input_file_transcripts, "output.txt\n\n")
cat("  # Annotate gene IDs (explicit type):\n")
cat("  Rscript annotate_transcripts.R", input_file_genes, "output.txt gene\n\n")

cat("Output Format:\n")
cat("  Tab-separated file with columns:\n")
cat("  - ensembl_gene_id\n")
cat("  - ensembl_transcript_id\n")
cat("  - transcript_biotype\n")
cat("  - external_gene_name\n\n")
