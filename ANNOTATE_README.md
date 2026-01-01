# Transcript Annotation Script

## Overview

`annotate_transcripts.R` is an R script that annotates Ensembl transcript IDs or gene IDs with biotype information from the Ensembl BioMart database. This script is designed to take specific IDs as input and return them with added annotation columns, making it useful for downstream analysis workflows.

## Key Differences from extract_gene_attributes.R

| Feature | extract_gene_attributes.R | annotate_transcripts.R |
|---------|--------------------------|------------------------|
| Purpose | Extract all genes from Ensembl | Annotate specific input IDs |
| Input | None (queries all genes) | File with specific IDs |
| Output | All gene attributes | Input IDs + annotations |
| Use Case | Generate reference database | Annotate your own data |

## Extracted Attributes

The script adds the following annotation columns:

1. **Ensembl Gene ID** (`ensembl_gene_id`) - Unique identifier for genes
2. **Ensembl Transcript ID** (`ensembl_transcript_id`) - Unique identifier for transcripts
3. **Transcript Biotype** (`transcript_biotype`) - Classification of the transcript
4. **Gene Name** (`external_gene_name`) - Common gene symbol

## Transcript Biotypes

Common biotype values include:

- `protein_coding` - Protein-coding genes
- `lncRNA` - Long non-coding RNA
- `processed_pseudogene` - Processed pseudogenes
- `miRNA` - microRNA
- `snRNA` - Small nuclear RNA
- `snoRNA` - Small nucleolar RNA
- `rRNA` - Ribosomal RNA
- `misc_RNA` - Miscellaneous RNA
- `processed_transcript` - Transcripts that don't contain an open reading frame

## Requirements

- **R** (version 3.5.0 or higher)
- **biomaRt** package from Bioconductor

## Installation

### Option 1: Install biomaRt via BiocManager (requires internet access)

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("biomaRt")
```

### Option 2: Install via apt (Debian/Ubuntu systems)

```bash
sudo apt-get install r-bioc-biomart
```

## Usage

### Basic Syntax

```bash
Rscript annotate_transcripts.R <input_file> <output_file> [id_type]
```

### Arguments

- **input_file** (required) - Path to input file with Ensembl IDs
- **output_file** (required) - Path to output file for annotated data
- **id_type** (optional) - Type of ID in input:
  - `auto` (default) - Automatically detect ID type
  - `transcript` - Input contains transcript IDs (ENST...)
  - `gene` - Input contains gene IDs (ENSG...)

### Examples

#### Example 1: Annotate Transcript IDs (auto-detect)

```bash
# Input file: transcripts.txt
# ENST00000373020
# ENST00000373031
# ENST00000371582

Rscript annotate_transcripts.R transcripts.txt annotated_transcripts.txt
```

#### Example 2: Annotate Gene IDs (explicit type)

```bash
# Input file: genes.txt
# ENSG00000000003
# ENSG00000000005
# ENSG00000000419

Rscript annotate_transcripts.R genes.txt annotated_genes.txt gene
```

#### Example 3: Run from within R

```r
source("annotate_transcripts.R")
result <- annotate_transcripts("input.txt", "output.txt", "auto")
```

## Input Format

The input file should contain Ensembl IDs in one of these formats:

### Simple List (one ID per line)

```
ENST00000373020
ENST00000373031
ENST00000371582
```

### Tab-separated File

```
ENST00000373020	sample1	data1
ENST00000373031	sample2	data2
ENST00000371582	sample3	data3
```

### Comma-separated File

```
ENST00000373020,sample1,data1
ENST00000373031,sample2,data2
ENST00000371582,sample3,data3
```

**Note:** The script will use the first column as the ID column. Additional columns will be preserved in the output.

## Output Format

The script generates a **tab-separated** file with **headers**. The output includes:

1. Original ID column (ensembl_transcript_id or ensembl_gene_id)
2. Additional annotation columns
3. Any additional columns from the input file (if present)

### Output for Transcript IDs

```
ensembl_transcript_id	ensembl_gene_id	transcript_biotype	external_gene_name
ENST00000373020	ENSG00000000003	protein_coding	TSPAN6
ENST00000373031	ENSG00000000005	protein_coding	TNMD
ENST00000371582	ENSG00000000419	protein_coding	DPM1
```

### Output for Gene IDs (Multiple Transcripts per Gene)

```
ensembl_gene_id	ensembl_transcript_id	transcript_biotype	external_gene_name
ENSG00000000003	ENST00000373020	protein_coding	TSPAN6
ENSG00000000003	ENST00000494424	processed_transcript	TSPAN6
ENSG00000000005	ENST00000373031	protein_coding	TNMD
```

**Important:** When annotating gene IDs, a single gene may produce multiple output rows (one for each transcript of that gene).

## Output Statistics

When the script runs, it displays summary statistics:

```
Reading input file: transcripts.txt
Auto-detected ID type: transcript
Connecting to Ensembl BioMart database...
Querying BioMart for annotations...
Processing batch 1/1 (5 IDs)...
Retrieved annotations for 5 records
Annotated data saved to: annotated_transcripts.txt

Summary:
  Input IDs: 5
  IDs with annotations: 5
  IDs without annotations: 0

Biotype distribution:
protein_coding    5
```

## Example Output

An example script (`example_annotate_transcripts.R`) is provided to demonstrate the expected output format without requiring internet access or R installation:

```bash
Rscript example_annotate_transcripts.R
```

This creates sample files:
- `example_transcript_ids.txt` - Sample input (transcript IDs)
- `example_transcript_annotations.txt` - Sample output (annotated transcripts)
- `example_gene_ids.txt` - Sample input (gene IDs)
- `example_gene_annotations.txt` - Sample output (annotated genes)

## Features

### Batch Processing

The script automatically processes large input files in batches of 500 IDs to avoid overwhelming the BioMart server. Progress is displayed for each batch.

### Auto-Detection

When `id_type` is set to `auto` (default), the script analyzes the input IDs:
- Counts IDs starting with `ENST` (transcripts)
- Counts IDs starting with `ENSG` (genes)
- Selects the predominant type

### Preserving Original Data

If your input file has multiple columns, the script preserves all columns in the output, making it easy to integrate annotations into existing datasets.

### Handling Missing Annotations

IDs without annotations in the BioMart database will appear in the output with empty annotation fields, allowing you to identify which IDs could not be annotated.

## Comparison with Other Tools

### vs. extract_gene_attributes.R

- **extract_gene_attributes.R**: Downloads all human genes from Ensembl (200,000+ records)
  - Use when: You need a complete reference database
  - Output: All genes and transcripts

- **annotate_transcripts.R**: Annotates your specific list of IDs
  - Use when: You have a list of IDs to annotate
  - Output: Only your input IDs with annotations

### Workflow Integration

```bash
# Example workflow:
# 1. Extract variants from your analysis
# 2. Get list of transcript IDs
# 3. Annotate them with biotypes

# Step 1: Your analysis produces transcript IDs
cat variants.txt | cut -f1 > transcript_ids.txt

# Step 2: Annotate with biotypes
Rscript annotate_transcripts.R transcript_ids.txt annotated.txt

# Step 3: Filter for protein-coding transcripts
awk '$3 == "protein_coding"' annotated.txt > protein_coding_variants.txt
```

## Important Notes

1. **Internet Connection Required**: The script needs internet access to connect to the Ensembl BioMart database at www.ensembl.org

2. **Execution Time**: Annotation time depends on the number of IDs:
   - Small lists (< 100 IDs): seconds
   - Medium lists (100-1000 IDs): 1-2 minutes
   - Large lists (> 1000 IDs): several minutes

3. **Ensembl Version**: The script connects to the current release of Ensembl. Results may vary slightly between releases as gene annotations are updated.

4. **Multiple Transcripts per Gene**: When annotating gene IDs, remember that genes can have multiple transcripts with different biotypes. The output will include all transcripts for each gene.

5. **Supported ID Formats**: Only Ensembl IDs are supported:
   - Transcript IDs: ENST00000000000
   - Gene IDs: ENSG00000000000
   - Other formats (RefSeq, UCSC, etc.) are not supported

## Troubleshooting

### Connection Errors

If you see connection errors:
```
Error: Failed to connect to Ensembl BioMart database.
```

Verify:
- Your internet connection is active
- You can access www.ensembl.org in your browser
- No firewall is blocking the connection
- The Ensembl server is not temporarily down

### Package Errors

If biomaRt is not found:
```bash
# For Debian/Ubuntu
sudo apt-get install r-bioc-biomart

# Or in R
BiocManager::install("biomaRt")
```

### No Annotations Retrieved

If no annotations are found:
```
Error: No annotations retrieved from BioMart.
```

Possible causes:
- IDs are not valid Ensembl IDs
- IDs are from an old Ensembl version
- The Ensembl server is temporarily unavailable

Check your IDs format:
```bash
# Valid formats:
ENST00000373020  # Transcript ID
ENSG00000000003  # Gene ID

# Invalid formats:
NM_001234        # RefSeq ID
uc001aaa.1       # UCSC ID
```

### Auto-Detection Fails

If auto-detection cannot determine ID type:
```
Error: Cannot auto-detect ID type.
```

Specify the ID type explicitly:
```bash
# For transcript IDs
Rscript annotate_transcripts.R input.txt output.txt transcript

# For gene IDs
Rscript annotate_transcripts.R input.txt output.txt gene
```

## Performance Tips

1. **Batch Large Requests**: The script automatically batches requests, but for very large lists (>10,000 IDs), consider splitting into multiple smaller files.

2. **Use Caching**: If annotating the same IDs multiple times, save the output and reuse it rather than querying BioMart repeatedly.

3. **Off-Peak Hours**: For very large annotation tasks, consider running during off-peak hours when the Ensembl server is less busy.

## Additional Information

- **Ensembl BioMart**: https://www.ensembl.org/biomart/
- **biomaRt Package Documentation**: https://bioconductor.org/packages/biomaRt/
- **Ensembl Release Notes**: https://www.ensembl.org/info/website/news.html
- **BioMart User Guide**: https://www.ensembl.org/info/data/biomart/index.html

## Support

For issues specific to this script, please refer to the repository issues page.

For BioMart-related questions, see the Ensembl help desk: https://www.ensembl.org/Help/Contact
