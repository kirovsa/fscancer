#!/usr/bin/perl
# sample_filter.pl - Functions for metadata-based per-sample filtering
# This module provides functionality to filter out model-derived samples
# (PDX, cell lines, etc.) from mutation files using metadata files.

use strict;
use warnings;

# Sample ID column names (case-insensitive matching)
my @SAMPLE_ID_COLUMNS = qw(
    sample_id sample tumor_sample_barcode sample_barcode
    samplebarcode sampleid tumor_sample_id
);

# Sample type column names (case-insensitive matching)
my @SAMPLE_TYPE_COLUMNS = qw(
    sample_type sampletype sample_type_detail model
    is_model sample_class sampleclass
);

# Model type patterns (case-insensitive matching)
my @MODEL_TYPE_PATTERNS = (
    'pdx',
    'patient.?derived.?xenograft',
    'xenograft',
    'cell.?line',
    'cellline',
    'in.?vitro',
    'model',
    'ccle'
);

# Metadata filenames to search for
my @METADATA_FILENAMES = qw(
    metadata.tsv samples.tsv clinical_sample.tsv
    sample_metadata.csv sample_annotations.tsv clinical.tsv
);

# Check if a sample type value indicates a model sample
# Returns 1 if model, 0 otherwise
sub is_model_sample {
    my ($sample_type) = @_;
    return 0 unless defined $sample_type && $sample_type ne '';
    
    my $value = lc($sample_type);
    $value =~ s/^\s+|\s+$//g;  # trim whitespace
    
    foreach my $pattern (@MODEL_TYPE_PATTERNS) {
        if ($value =~ /$pattern/) {
            return 1;
        }
    }
    return 0;
}

# Find column index matching any of the given patterns (case-insensitive)
sub find_column_index {
    my ($headers_ref, $patterns_ref) = @_;
    my @headers = @$headers_ref;
    my @patterns = @$patterns_ref;
    
    for my $i (0..$#headers) {
        my $header = lc($headers[$i]);
        $header =~ s/^\s+|\s+$//g;
        foreach my $pattern (@patterns) {
            if ($header eq lc($pattern)) {
                return $i;
            }
        }
    }
    return undef;
}

# Detect delimiter (tab or comma)
sub detect_delimiter {
    my ($first_line) = @_;
    my $tab_count = ($first_line =~ tr/\t//);
    my $comma_count = ($first_line =~ tr/,//);
    return $tab_count >= $comma_count ? "\t" : ",";
}

# Discover metadata files in a directory
sub discover_metadata_files {
    my ($directory) = @_;
    my @found_files;
    
    return @found_files unless -d $directory;
    
    # Check for exact filename matches
    foreach my $filename (@METADATA_FILENAMES) {
        my $path = "$directory/$filename";
        if (-f $path) {
            push @found_files, $path;
        }
    }
    
    # Check for files containing "sample" or "metadata" in name
    opendir(my $dh, $directory) or return @found_files;
    while (my $file = readdir($dh)) {
        next if $file =~ /^\./;
        next unless $file =~ /\.(tsv|csv)$/i;
        if ($file =~ /sample|metadata|clinical/i) {
            my $path = "$directory/$file";
            push @found_files, $path unless grep { $_ eq $path } @found_files;
        }
    }
    closedir($dh);
    
    return @found_files;
}

# Parse a metadata file and return hash of sample_id => is_model
sub parse_metadata_file {
    my ($file_path) = @_;
    my %sample_is_model;
    
    open(my $fh, '<', $file_path) or do {
        warn "Warning: Could not open metadata file $file_path: $!";
        return %sample_is_model;
    };
    
    # Skip comment lines and find header
    my $header_line;
    while (my $line = <$fh>) {
        chomp $line;
        next if $line =~ /^#/;
        $header_line = $line;
        last;
    }
    
    return %sample_is_model unless defined $header_line;
    
    my $delimiter = detect_delimiter($header_line);
    my @headers = split(/$delimiter/, $header_line);
    
    # Find sample ID and sample type columns
    my $sample_id_idx = find_column_index(\@headers, \@SAMPLE_ID_COLUMNS);
    my $sample_type_idx = find_column_index(\@headers, \@SAMPLE_TYPE_COLUMNS);
    
    return %sample_is_model unless defined $sample_id_idx && defined $sample_type_idx;
    
    # Parse data rows
    while (my $line = <$fh>) {
        chomp $line;
        next if $line =~ /^#/ || $line =~ /^\s*$/;
        
        my @fields = split(/$delimiter/, $line);
        
        my $max_idx = $sample_id_idx > $sample_type_idx ? $sample_id_idx : $sample_type_idx;
        next if scalar(@fields) <= $max_idx;
        
        my $sample_id = $fields[$sample_id_idx];
        my $sample_type = $fields[$sample_type_idx];
        
        $sample_id =~ s/^\s+|\s+$//g if defined $sample_id;
        
        if (defined $sample_id && $sample_id ne '') {
            $sample_is_model{$sample_id} = is_model_sample($sample_type);
        }
    }
    
    close($fh);
    return %sample_is_model;
}

# Load sample metadata from a file or directory
# Returns hash of sample_id => is_model (1 or 0)
sub load_sample_metadata {
    my ($path_or_dir, @additional_paths) = @_;
    my %sample_is_model;
    my @metadata_files;
    
    if (-f $path_or_dir) {
        push @metadata_files, $path_or_dir;
    } elsif (-d $path_or_dir) {
        push @metadata_files, discover_metadata_files($path_or_dir);
    }
    
    foreach my $path (@additional_paths) {
        if (-f $path && !grep { $_ eq $path } @metadata_files) {
            push @metadata_files, $path;
        }
    }
    
    foreach my $metadata_file (@metadata_files) {
        my %file_mapping = parse_metadata_file($metadata_file);
        while (my ($k, $v) = each %file_mapping) {
            $sample_is_model{$k} = $v;
        }
    }
    
    return %sample_is_model;
}

# Get sample ID column index from mutation file headers
sub get_sample_id_column_index {
    my ($headers_ref) = @_;
    return find_column_index($headers_ref, \@SAMPLE_ID_COLUMNS);
}

# Check if a sample barcode is a model sample
# Returns 1 if model, 0 if not, undef if not in metadata
sub is_sample_model {
    my ($sample_barcode, $model_samples_ref) = @_;
    return $model_samples_ref->{$sample_barcode} if exists $model_samples_ref->{$sample_barcode};
    return undef;
}

# Get set of model sample IDs from metadata hash
sub get_model_samples {
    my ($metadata_ref) = @_;
    my %model_samples;
    while (my ($sample, $is_model) = each %$metadata_ref) {
        $model_samples{$sample} = 1 if $is_model;
    }
    return %model_samples;
}

1;  # Return true for module loading
