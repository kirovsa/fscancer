# fscancer
This code was made available to reproduce the frameshift analysis in the publication:
"Framshifts may carry oncogenic potential beyond loss of function"
There are 2 components:
Produce intermediate files from the perl script
combineCBIO.fs.pl
Script should be run from the github clone of CBIO one level above public directory. No params needed, jst redirect output to 
cbio.fs.nocelllinepdx.txt
This is read by the Rmd script that generates the report used to analyze the data that was published
