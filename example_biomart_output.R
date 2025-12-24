#!/usr/bin/env Rscript
# Example demonstrating the extract_gene_attributes.R script functionality
# This creates sample output showing what the script produces when run with internet access

# Create sample gene attribute data matching the expected format
# Column order: ENST, ENSG, Symbol, ENSP, ProtLen, Type
sample_data <- data.frame(
  ensembl_transcript_id = c(
    "ENST00000373020", "ENST00000373031", "ENST00000371582",
    "ENST00000367770", "ENST00000286031", "ENST00000374003",
    "ENST00000374080", "ENST00000371475", "ENST00000367429",
    "ENST00000358731"
  ),
  ensembl_gene_id = c(
    "ENSG00000000003", "ENSG00000000005", "ENSG00000000419",
    "ENSG00000000457", "ENSG00000000460", "ENSG00000000938",
    "ENSG00000001036", "ENSG00000001084", "ENSG00000001167",
    "ENSG00000001460"
  ),
  external_gene_name = c(
    "TSPAN6", "TNMD", "DPM1", "SCYL3", "C1orf112",
    "FGR", "FUCA2", "GCLC", "NFYA", "STPG1"
  ),
  ensembl_peptide_id = c(
    "ENSP00000362111", "ENSP00000362122", "ENSP00000360644",
    "ENSP00000356839", "ENSP00000275493", "ENSP00000363079",
    "ENSP00000363155", "ENSP00000360542", "ENSP00000356496",
    "ENSP00000347939"
  ),
  cds_length = c(
    879, 1455, 834, 2268, 888,
    1584, 1437, 1860, 1089, 1149
  ),
  transcript_biotype = c(
    "protein_coding", "protein_coding", "protein_coding",
    "protein_coding", "protein_coding", "protein_coding",
    "protein_coding", "protein_coding", "protein_coding",
    "protein_coding"
  ),
  stringsAsFactors = FALSE
)

# Write to example file
output_file <- "example_gene_attributes.txt"
write.table(
  sample_data,
  file = output_file,
  sep = "|",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE,
  na = ""
)

cat("Example output file created:", output_file, "\n")
cat("\nThis demonstrates the format of the output from extract_gene_attributes.R\n")
cat("\nSample output (first 3 lines):\n")
cat(paste(readLines(output_file, n = 3), collapse = "\n"), "\n")
cat("\nFormat: ensembl_transcript_id|ensembl_gene_id|gene_name|ensembl_peptide_id|cds_length|transcript_biotype\n")
cat("Matches: ENST|ENSG|Symbol|ENSP|ProtLen|Type (as expected by frameshift.Rmd)\n")
