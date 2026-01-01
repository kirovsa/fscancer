# Comparison: biomaRt vs ensembldb for Gene Attributes Extraction

This document helps you choose between the two gene attributes extraction scripts.

## Quick Comparison

| Feature | biomaRt (original) | ensembldb (alternative) |
|---------|-------------------|------------------------|
| **Script** | `extract_gene_attributes.R` | `extract_gene_attributes_ensembldb.R` |
| **Package** | biomaRt | ensembldb |
| **Internet Required** | ✅ Always | ❌ Only for AnnotationHub install |
| **Works Offline** | ❌ | ✅ |
| **Speed** | Slower (minutes) | Faster (seconds) |
| **Query Type** | Remote API calls | Local database queries |
| **Reproducibility** | Version may vary | Fixed version |
| **Server Dependency** | Yes (www.ensembl.org) | No |
| **Data Source** | Current Ensembl release | Specific EnsDb package version |
| **Output Format** | Pipe-separated (no header) | Pipe-separated (no header) |
| **Output Columns** | 6 attributes | 6 attributes (identical) |
| **Compatible with frameshift.Rmd** | ✅ | ✅ |

## Output Format

Both scripts produce **identical** output formats:

```
ENST00000373020|ENSG00000000003|TSPAN6|ENSP00000362111|879|protein_coding
ENST00000373031|ENSG00000000005|TNMD|ENSP00000362122|1455|protein_coding
```

Format: `transcript_id|gene_id|gene_name|peptide_id|cds_length|biotype`

## When to Use biomaRt

**Best for:**
- One-time data extraction
- Need the absolute latest Ensembl data
- Don't want to install large annotation packages (~100-500 MB)
- Quick analysis without local database setup
- Following the original workflow from the paper

**Command:**
```bash
Rscript extract_gene_attributes.R gene_attributes.txt
```

## When to Use ensembldb

**Best for:**
- Automated pipelines and workflows
- Offline or restricted network environments
- Need reproducible results with fixed Ensembl versions
- Multiple extractions (faster after initial setup)
- Integration with other ensembldb-based analyses
- Citation of specific Ensembl release in publications

**Commands:**
```bash
# Default (EnsDb.Hsapiens.v86):
Rscript extract_gene_attributes_ensembldb.R gene_attributes.txt

# Specific version:
Rscript extract_gene_attributes_ensembldb.R gene_attributes.txt EnsDb.Hsapiens.v109

# Latest from AnnotationHub:
Rscript extract_gene_attributes_ensembldb.R gene_attributes.txt hub
```

## Installation Comparison

### biomaRt Setup
```r
# Simple, one package
BiocManager::install("biomaRt")
```

**Size:** ~5 MB  
**Setup time:** < 1 minute  
**Network:** Required for installation and execution

### ensembldb Setup
```r
# Two packages needed
BiocManager::install("ensembldb")
BiocManager::install("EnsDb.Hsapiens.v86")  # or other version
```

**Size:** ~5 MB (ensembldb) + ~100-500 MB (EnsDb package)  
**Setup time:** 5-10 minutes (depending on EnsDb size)  
**Network:** Required only for installation

## Performance Comparison

Typical execution times for extracting all human genes:

| Aspect | biomaRt | ensembldb |
|--------|---------|-----------|
| **First run** | 5-15 minutes | 10-30 seconds |
| **Subsequent runs** | 5-15 minutes | 10-30 seconds |
| **Network issues** | May fail | Unaffected |
| **Server load** | Variable | Consistent |

## Data Version Control

### biomaRt
- Always queries current Ensembl release
- Version can change between runs
- Less reproducible for publications
- Good for latest annotations

### ensembldb
- Fixed version based on installed EnsDb package
- Consistent results over time
- Better for reproducible research
- Can cite specific Ensembl release

## Available Ensembl Versions

### Via ensembldb EnsDb packages:
- **v75** (Ensembl 75, GRCh37/hg19) - For older genome build
- **v79** (Ensembl 79, GRCh38) - Early GRCh38
- **v86** (Ensembl 86, GRCh38) - Stable, widely used (default)
- **v109** (Ensembl 109, GRCh38) - More recent

### Via AnnotationHub:
- Automatically downloads latest available version
- Cached locally for future use

## Code Comparison

### biomaRt approach:
```r
library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_data <- getBM(attributes = c(...), mart = ensembl)
```

### ensembldb approach:
```r
library(ensembldb)
library(EnsDb.Hsapiens.v86)
edb <- EnsDb.Hsapiens.v86
tx_data <- transcripts(edb, columns = c(...))
```

## Troubleshooting

### biomaRt Common Issues:
- "Connection timeout" → Server unavailable, try again later
- "SSL certificate problem" → Certificate verification issues
- "Empty result" → Server may be down

### ensembldb Common Issues:
- "Package not found" → Install EnsDb package first
- "Database version mismatch" → Install correct EnsDb version
- "No such column" → Update ensembldb package

## Recommendations

### For Paper Reproduction:
Use **biomaRt** (original) to match the publication workflow, or specify an exact Ensembl version with ensembldb for reproducibility.

### For New Analyses:
Use **ensembldb** for better reproducibility and offline capability. Document the specific EnsDb version used.

### For Production Pipelines:
Use **ensembldb** to avoid network dependencies and ensure consistent results.

### For Exploratory Analysis:
Either works fine. **biomaRt** has simpler setup; **ensembldb** is faster for repeated runs.

## Migration Between Approaches

The output formats are identical, so you can:
1. Switch between scripts without changing downstream code
2. Compare outputs from different Ensembl versions
3. Use both scripts in the same analysis pipeline

## Documentation

- **biomaRt:** See [BIOMART_README.md](BIOMART_README.md)
- **ensembldb:** See [ENSEMBLDB_README.md](ENSEMBLDB_README.md)

## Examples

Both scripts include example generators that work without databases:
- `example_biomart_output.R`
- `example_ensembldb_output.R`

These create sample output files showing the expected format.
