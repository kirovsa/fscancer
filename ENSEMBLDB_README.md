# Gene Attributes Extraction Script (ensembldb version)

## Overview

`extract_gene_attributes_ensembldb.R` is an R script that uses the **ensembldb** package to extract gene attributes from local Ensembl databases for all human genes. This is an alternative to `extract_gene_attributes.R` which uses biomaRt.

## Key Advantages over biomaRt

1. **Works Offline** - No internet connection required after initial package installation
2. **Faster Queries** - Uses local databases instead of remote API calls
3. **More Reproducible** - Fixed Ensembl version for consistent results
4. **No Server Dependencies** - Independent of remote server availability
5. **Better for Automation** - More reliable in scripts and pipelines

## Extracted Attributes

The script retrieves the following attributes for each transcript (in this order):

1. **Transcript Stable ID** (`tx_id`) - Unique identifier for transcripts
2. **Stable Ensembl Gene ID** (`gene_id`) - Unique identifier for genes in Ensembl
3. **Gene Name** (`gene_name`) - Common gene symbol (e.g., TP53, BRCA1)
4. **Peptide Stable ID** (`protein_id`) - Unique identifier for protein sequences
5. **CDS Length** (`cds_length`) - Length of the coding sequence in base pairs
6. **Transcript Type** (`tx_biotype`) - Classification of the transcript (e.g., protein_coding, lncRNA)

## Requirements

- **R** (version 3.5.0 or higher)
- **ensembldb** package from Bioconductor
- An **EnsDb annotation package** (see Installation section)

## Installation

### Step 1: Install ensembldb package

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ensembldb")
```

### Step 2: Install an EnsDb annotation package

You have several options:

#### Option A: Install a specific Ensembl version (Recommended for reproducibility)

```r
# For Ensembl release 86 (GRCh38) - Good balance of stability and recency
BiocManager::install("EnsDb.Hsapiens.v86")

# For Ensembl release 109 (GRCh38) - More recent
BiocManager::install("EnsDb.Hsapiens.v109")

# For Ensembl release 75 (GRCh37/hg19) - If you need older genome build
BiocManager::install("EnsDb.Hsapiens.v75")
```

#### Option B: Use AnnotationHub to auto-download latest version

```r
BiocManager::install("AnnotationHub")
```

This allows the script to automatically download the latest EnsDb when you use the `hub` option.

### Available EnsDb Packages

| Package | Ensembl Version | Genome Build | Description |
|---------|----------------|--------------|-------------|
| EnsDb.Hsapiens.v75 | 75 | GRCh37 (hg19) | Older genome build |
| EnsDb.Hsapiens.v79 | 79 | GRCh38 | Early GRCh38 release |
| EnsDb.Hsapiens.v86 | 86 | GRCh38 | Stable, widely used |
| EnsDb.Hsapiens.v109 | 109 | GRCh38 | Recent release |

## Usage

### Basic usage with default settings:

```bash
Rscript extract_gene_attributes_ensembldb.R
```

This will:
- Use `EnsDb.Hsapiens.v86` as the database
- Create a file named `gene_attributes_ensembldb.txt`

### Specify custom output file:

```bash
Rscript extract_gene_attributes_ensembldb.R my_output.txt
```

### Use a specific EnsDb version:

```bash
Rscript extract_gene_attributes_ensembldb.R output.txt EnsDb.Hsapiens.v109
```

### Use AnnotationHub for latest version:

```bash
Rscript extract_gene_attributes_ensembldb.R output.txt hub
```

### Run from within R:

```r
source("extract_gene_attributes_ensembldb.R")

# Use default database
result <- extract_gene_attributes_ensembldb("output.txt")

# Use specific version
result <- extract_gene_attributes_ensembldb("output.txt", "EnsDb.Hsapiens.v109")

# Use AnnotationHub
result <- extract_gene_attributes_ensembldb("output.txt", "hub")
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

**This format is identical to the biomaRt version and fully compatible with frameshift.Rmd.**

## Output Statistics

When the script runs, it displays summary statistics including:

- Database version and genome build information
- Total number of records retrieved
- Number of unique genes
- Number of unique transcripts
- Number of records with peptide IDs
- Number of records with CDS length
- Distribution of transcript types

### Example output:

```
Loading EnsDb.Hsapiens.v86 ...
Database version: 86
Genome build: GRCh38
Organism: Homo sapiens

Fetching gene attributes...
Retrieved 216227 records
Gene attributes saved to: gene_attributes_ensembldb.txt

Summary:
  Total records: 216227
  Unique genes: 63970
  Unique transcripts: 216227
  Records with peptide IDs: 100836
  Records with CDS length: 100836

Transcript type distribution:
protein_coding           100836
processed_pseudogene      11168
lncRNA                    15778
...
```

## Compatibility with frameshift.Rmd

This script generates output in the exact same format as the biomaRt version (`extract_gene_attributes.R`). The output is compatible with `frameshift.Rmd`:

```r
enslen <- read.table("gene_attributes_ensembldb.txt", sep="|", head=F)
colnames(enslen) <- c("ENST","ENSG","Symbol","ENSP","ProtLen","Type")
```

**Column mapping:**
1. ENST = Transcript ID (`ensembl_transcript_id`)
2. ENSG = Gene ID (`ensembl_gene_id`)
3. Symbol = Gene Name (`external_gene_name`)
4. ENSP = Peptide ID (`ensembl_peptide_id`)
5. ProtLen = CDS Length (`cds_length`)
6. Type = Transcript Type (`transcript_biotype`)

## Comparison with biomaRt Version

| Feature | biomaRt | ensembldb |
|---------|---------|-----------|
| **Internet required** | Yes (always) | Only for AnnotationHub |
| **Speed** | Slower (API calls) | Faster (local queries) |
| **Reproducibility** | Version may change | Fixed version |
| **Server dependency** | Yes | No |
| **Offline usage** | No | Yes (after install) |
| **Database version** | Latest from server | User-specified |
| **Output format** | Pipe-separated | Pipe-separated (identical) |
| **Attributes** | Same 6 attributes | Same 6 attributes |

## Choosing Between biomaRt and ensembldb

**Use biomaRt (`extract_gene_attributes.R`) when:**
- You want the absolute latest Ensembl data
- You don't mind internet dependency
- You're doing a one-time extraction
- You want to avoid installing large annotation packages

**Use ensembldb (`extract_gene_attributes_ensembldb.R`) when:**
- You need offline capability
- You're running automated pipelines
- You need reproducible results with fixed Ensembl versions
- Speed is important
- You're doing multiple extractions

## Troubleshooting

### Package not found errors:

If you see errors about missing packages:

```r
# Install the required package
BiocManager::install("EnsDb.Hsapiens.v86")
```

### AnnotationHub errors:

If using `hub` option fails:

```r
# Clear cache and try again
library(AnnotationHub)
ah <- AnnotationHub()
cache(ah)  # Show cache location
# You may need to clear the cache directory if corrupted
```

### Empty results:

If no data is retrieved:
- Verify the EnsDb package is correctly installed
- Try a different EnsDb version
- Check that the package loaded successfully

### Version mismatch:

If you need a specific Ensembl version:
- Check available packages: https://bioconductor.org/packages/release/BiocViews.html#___EnsDb
- Install the exact version you need
- Use that version as the second argument to the script

## Additional Information

- **ensembldb Package**: https://bioconductor.org/packages/ensembldb/
- **EnsDb Packages**: https://bioconductor.org/packages/release/BiocViews.html#___EnsDb
- **Ensembl Versions**: https://www.ensembl.org/info/website/archives/index.html
- **AnnotationHub**: https://bioconductor.org/packages/AnnotationHub/

## Example Scripts

An example script (`example_ensembldb_output.R`) demonstrates the expected output format without requiring a database connection.
