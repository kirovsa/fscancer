# fscancer
This code was made available to reproduce the frameshift analysis in the publication:
"Frameshifts may carry oncogenic potential beyond loss of function"

## Components

### 1. Perl Script (Original)
**File:** `combineCBIO.fs.pl`

Script should be run from the github clone of cBioPortal one level above public directory. No params needed, just redirect output to `cbio.fs.nocelllinepdx.txt`.

This is read by the Rmd script that generates the report used to analyze the data that was published.

### 2. Python Module (New)
**Files:** `src/sample_filter.py`, `src/combine_mutations.py`

A Python-based alternative that provides enhanced per-sample filtering using metadata files. This allows filtering out model samples even from mixed mutation files that contain both patient and model samples.

## Filtering Model Data

### Path-Based Filtering (Original Behavior)
Both scripts filter out model data based on study path names containing:
- `ccle` - Cancer Cell Line Encyclopedia
- `pdx` - Patient-Derived Xenograft models
- `cellline` or `cell_line` - Cell line studies
- `xenograft` - Xenograft models
- `test` - Test datasets

### Metadata-Based Per-Sample Filtering (New Behavior)
The Python module additionally supports per-sample filtering using metadata files. This is useful when a mutation file contains a mix of patient and model samples.

**How it works:**
1. The script auto-discovers metadata files in the input directory
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
- Any file containing `sample` or `metadata` in the name (`.tsv` or `.csv`)

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

## Usage

### Perl Script (Original)
```bash
cd /path/to/cbioportal/clone
perl combineCBIO.fs.pl > cbio.fs.nocelllinepdx.txt
```

### Python Script (New)
```bash
# Basic usage with auto-discovery of metadata
python src/combine_mutations.py --input-dir /path/to/cbioportal/data --output output.txt

# With explicit metadata file
python src/combine_mutations.py --input-dir /path/to/data --metadata /path/to/metadata.tsv --output output.txt

# Verbose mode (shows filtering statistics)
python src/combine_mutations.py --input-dir /path/to/data --output output.txt --verbose

# Include model samples (disable filtering)
python src/combine_mutations.py --input-dir /path/to/data --output output.txt --include-model
```

### Python API
```python
from src.sample_filter import load_sample_metadata, filter_mutation_rows, is_model_sample

# Load metadata and get model samples
metadata = load_sample_metadata('/path/to/metadata.tsv')
model_samples = {s for s, is_model in metadata.items() if is_model}

# Filter a mutation file
rows_kept, rows_filtered = filter_mutation_rows(
    '/path/to/mutations.tsv',
    model_samples,
    '/path/to/filtered_output.tsv'
)
```

## Running Tests
```bash
cd /path/to/fscancer
python -m unittest tests.test_sample_filter -v
```

## Backwards Compatibility
- The original Perl script behavior is unchanged
- The Python script maintains path-based filtering for compatibility
- If no metadata files are found, only path-based filtering is applied
- Mixed mutation files (with both patient and model samples) are now supported with per-sample filtering
