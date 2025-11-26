# fscancer
This code was made available to reproduce the frameshift analysis in the publication:
"Frameshifts may carry oncogenic potential beyond loss of function"

## Components

### 1. Original Script
**File:** `combineCBIO.fs.pl`

Script should be run from the github clone of cBioPortal one level above public directory. No params needed, just redirect output to `cbio.fs.nocelllinepdx.txt`.

This is read by the Rmd script that generates the report used to analyze the data that was published.

### 2. Enhanced Script with Per-Sample Filtering
**File:** `combineCBIO.fs.filter.pl`

Enhanced version that supports metadata-based per-sample filtering. Useful when mutation files contain a mix of patient and model samples.

**Usage:**
```bash
# Basic usage with auto-discovery of metadata
perl combineCBIO.fs.filter.pl > output.txt

# With explicit metadata file
perl combineCBIO.fs.filter.pl --metadata /path/to/metadata.tsv > output.txt

# Verbose mode (shows filtering statistics)
perl combineCBIO.fs.filter.pl --verbose > output.txt

# Include model samples (disable filtering)
perl combineCBIO.fs.filter.pl --include-model > output.txt
```

### 3. Sample Filter Module
**File:** `sample_filter.pl`

Perl module providing functions for metadata-based sample filtering. Used by `combineCBIO.fs.filter.pl`.

## Filtering Model Data

### Path-Based Filtering (Original Behavior)
Both scripts filter out model data based on study path names containing:
- `ccle` - Cancer Cell Line Encyclopedia
- `pdx` - Patient-Derived Xenograft models
- `cellline` or `cell_line` - Cell line studies
- `xenograft` - Xenograft models
- `test` - Test datasets

### Metadata-Based Per-Sample Filtering (New)
The enhanced script (`combineCBIO.fs.filter.pl`) supports per-sample filtering using metadata files. This is useful when a mutation file contains a mix of patient and model samples.

**How it works:**
1. The script auto-discovers metadata files in the current directory
2. Metadata files map sample barcodes to sample types
3. Mutation rows belonging to model samples are filtered out during merge

#### Supported Metadata Filenames
The script looks for these files (in order of priority):
- `metadata.tsv`
- `samples.tsv`
- `clinical_sample.tsv`
- `sample_metadata.csv`
- `sample_annotations.tsv`
- `clinical.tsv`
- Any file containing `sample`, `metadata`, or `clinical` in the name (`.tsv` or `.csv`)

#### Required Metadata Columns
The metadata file must have both:

**Sample ID Column** (one of, case-insensitive):
- `sample_id`
- `sample`
- `Tumor_Sample_Barcode`
- `sample_barcode`
- `sampleBarcode`

**Sample Type Column** (one of, case-insensitive):
- `sample_type`
- `sampleType`
- `sample_type_detail`
- `model`
- `is_model`
- `sample_class`

#### Model Sample Type Values
A sample is marked as a model if its type contains (case-insensitive):
- `pdx`
- `patient-derived xenograft`
- `xenograft`
- `cell line`, `cell_line`, or `cellline`
- `in vitro`
- `model`
- `ccle`

#### Example Metadata File
```tsv
sample_id	sample_type
SAMPLE001	Patient
SAMPLE002	PDX
SAMPLE003	Primary Tumor
SAMPLE004	Cell Line
SAMPLE005	Metastatic
```

## Running Tests
```bash
perl t/test_sample_filter.pl
```

## Backwards Compatibility
- The original `combineCBIO.fs.pl` script is unchanged
- The enhanced script maintains path-based filtering for compatibility
- If no metadata files are found, only path-based filtering is applied
- Mixed mutation files (with both patient and model samples) are now supported with per-sample filtering
