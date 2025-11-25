# fscancer
This code was made available to reproduce the frameshift analysis in the publication:
"Framshifts may carry oncogenic potential beyond loss of function"
There are 2 components:
Produce intermediate files from the perl script
combineCBIO.fs.pl
Script should be run from the github clone of CBIO one level above public directory. No params needed, just redirect output to 
cbio.fs.nocelllinepdx.txt
This is read by the Rmd script that generates the report used to analyze the data that was published

## Filtering Model Data
The script filters out all model data (cell lines and PDX models) from the mutation data by excluding studies containing:
- `ccle` - Cancer Cell Line Encyclopedia
- `pdx` - Patient-Derived Xenograft models
- `cellline` or `cell_line` - Cell line studies
- `xenograft` - Xenograft models
- `test` - Test datasets
