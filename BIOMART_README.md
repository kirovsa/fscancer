# Gene Attributes Extraction Script

## Overview

`extract_gene_attributes.R` is an R script that uses the biomaRt package to extract gene attributes from the Ensembl BioMart database for all human genes.

## Extracted Attributes

The script retrieves the following attributes for each gene (in this order):

1. **Transcript Stable ID** (`ensembl_transcript_id`) - Unique identifier for transcripts
2. **Stable Ensembl Gene ID** (`ensembl_gene_id`) - Unique identifier for genes in Ensembl
3. **Gene Name** (`external_gene_name`) - Common gene symbol (e.g., TP53, BRCA1)
4. **Peptide Stable ID** (`ensembl_peptide_id`) - Unique identifier for protein sequences
5. **CDS Length** (`cds_length`) - Length of the coding sequence in base pairs
6. **Transcript Type** (`transcript_biotype`) - Classification of the transcript (e.g., protein_coding, lncRNA)

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

### Basic usage (default output file):

```bash
Rscript extract_gene_attributes.R
```

This will create a file named `gene_attributes.txt` in the current directory.

### Specify custom output file:

```bash
Rscript extract_gene_attributes.R my_output.txt
```

### Run from within R:

```r
source("extract_gene_attributes.R")
result <- extract_gene_attributes("output.txt")
```

## Output Format

The script generates a pipe-separated (`|`) file with **no header**. Each line contains the six attributes in this order:

```
ensembl_transcript_id|ensembl_gene_id|gene_name|ensembl_peptide_id|cds_length|transcript_biotype
```

### Example output:

```
ENST00000373020|ENSG00000000003|TSPAN6|ENSP00000362111|879|protein_coding
ENST00000373031|ENSG00000000005|TNMD|ENSP00000362122|1455|protein_coding
ENST00000371582|ENSG00000000419|DPM1|ENSP00000360644|834|protein_coding
```

### Note on empty fields:

- Some transcripts may not have peptide IDs (e.g., non-coding RNAs)
- Some transcripts may not have CDS length (e.g., pseudogenes)
- Empty fields are represented as blank in the output

## Output Statistics

When the script runs, it displays summary statistics including:

- Total number of records retrieved
- Number of unique genes
- Number of unique transcripts
- Number of records with peptide IDs
- Number of records with CDS length
- Distribution of transcript types

### Example output:

```
Connecting to Ensembl BioMart database...
Fetching gene attributes...
Retrieved 245031 records
Gene attributes saved to: gene_attributes.txt

Summary:
  Total records: 245031
  Unique genes: 71520
  Unique transcripts: 245031
  Records with peptide IDs: 115234
  Records with CDS length: 115234

Transcript type distribution:
protein_coding                 115234
processed_pseudogene            11234
lncRNA                          18345
...
```

## Example Output

An example script (`example_biomart_output.R`) is provided to demonstrate the expected output format. Run it with:

```bash
Rscript example_biomart_output.R
```

This creates a sample file `example_gene_attributes.txt` showing the format without requiring internet access.

## Compatibility with frameshift.Rmd

This script generates output in the same format as the `protlen.txt` file referenced in `frameshift.Rmd`. The pipe-separated format with column order matches the expected input:

```r
enslen<-read.table("gene_attributes.txt",sep="|",head=F)
colnames(enslen)<-c("ENST","ENSG","Symbol","ENSP","ProtLen","Type")
```

**Column mapping:**
1. ENST = Transcript ID (`ensembl_transcript_id`)
2. ENSG = Gene ID (`ensembl_gene_id`)
3. Symbol = Gene Name (`external_gene_name`)
4. ENSP = Peptide ID (`ensembl_peptide_id`)
5. ProtLen = CDS Length (`cds_length`)
6. Type = Transcript Type (`transcript_biotype`)

## Important Notes

1. **Internet Connection Required**: The script needs internet access to connect to the Ensembl BioMart database at www.ensembl.org

2. **Mirror Fallback**: If the main Ensembl server times out or is unavailable, the script automatically tries the US East mirror (useast.ensembl.org) as a fallback

3. **Execution Time**: Querying all human genes may take several minutes depending on network speed and server load

4. **Ensembl Version**: The script connects to the current release of Ensembl. Results may vary slightly between releases as gene annotations are updated

5. **Multiple Transcripts**: Genes can have multiple transcripts, so the same gene ID may appear multiple times with different transcript IDs

## Troubleshooting

### Connection errors:

If you see connection errors, the script will automatically try alternative mirrors. If all mirrors fail, verify:
- Your internet connection is active
- You can access www.ensembl.org or useast.ensembl.org in your browser
- No firewall is blocking the connection

### Package errors:

If biomaRt is not found:
```bash
# For Debian/Ubuntu
sudo apt-get install r-bioc-biomart

# Or in R
BiocManager::install("biomaRt")
```

### Empty results:

If no data is retrieved from any server:
- The Ensembl servers (both main and mirrors) may be temporarily unavailable
- Try again later or check the Ensembl status page

## Additional Information

- **Ensembl BioMart**: https://www.ensembl.org/biomart/
- **biomaRt Package Documentation**: https://bioconductor.org/packages/biomaRt/
- **Ensembl Release Notes**: https://www.ensembl.org/info/website/news.html
