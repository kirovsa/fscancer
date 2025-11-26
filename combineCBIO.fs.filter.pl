#!/usr/bin/perl
# combineCBIO.fs.filter.pl - Enhanced mutation file merger with per-sample filtering
#
# This script extends combineCBIO.fs.pl with metadata-based per-sample filtering.
# It filters out model-derived samples (PDX, cell lines, etc.) using metadata files.
#
# Usage:
#   perl combineCBIO.fs.filter.pl [--metadata FILE] [--verbose] > output.txt
#
# If metadata files are found in the current directory or specified via --metadata,
# individual mutation rows belonging to model samples will be filtered out.
# Otherwise, falls back to path-based filtering only (original behavior).

use strict;
use warnings;
use Getopt::Long;

# Include sample_filter functions
require "./sample_filter.pl";

# Parse command line arguments
my $metadata_file = '';
my $verbose = 0;
my $include_model = 0;

GetOptions(
    'metadata=s'    => \$metadata_file,
    'verbose'       => \$verbose,
    'include-model' => \$include_model,
) or die "Usage: $0 [--metadata FILE] [--verbose] [--include-model]\n";

# Load metadata for per-sample filtering
my %model_samples;
my $has_metadata = 0;

unless ($include_model) {
    # Auto-discover metadata files in current directory
    my %metadata = load_sample_metadata('.');
    
    # Also load from specified metadata file
    if ($metadata_file && -f $metadata_file) {
        my %extra = parse_metadata_file($metadata_file);
        while (my ($k, $v) = each %extra) {
            $metadata{$k} = $v;
        }
    }
    
    # Extract model samples
    %model_samples = get_model_samples(\%metadata);
    $has_metadata = scalar(keys %model_samples) > 0;
    
    if ($verbose) {
        if ($has_metadata) {
            print STDERR "Loaded " . scalar(keys %model_samples) . " model samples from metadata\n";
        } else {
            print STDERR "No metadata files found, using path-based filtering only\n";
        }
    }
}

# Find mutation files with path-based filtering (original behavior)
my $cbfs=`find . -name 'data_mutations*' |grep -vi ccle|grep -vi pdx|grep -vi cellline|grep -vi cell_line|grep -vi xenograft|grep -vi test`;
open(CNT,">cnt");
my @cbf=split(/\n/,$cbfs);

my %seen;
my %cnt;
my $total_filtered = 0;
my $total_kept = 0;

foreach my $proj (@cbf) {
    my ($dor,$fm,$projname,$file)=split(/\//,$proj);
    my ($study,$center,$r)=split("\_",$projname);
    my $uid=$study.$center;
    
    # Skip duplicates
    next if ($seen{$uid}++);
    
    open(F,$proj) or do {
        warn "Could not open $proj: $!";
        next;
    };
    
    # Skip comment lines to find header
    my $h;
    do {
        $h=<F>;
    } until (!defined($h) || $h!~/^#/);
    
    next unless defined $h;
    chomp $h;
    
    my @h=split(/\t/,$h);
    my $gene;
    my $vtype;
    my $classification;
    my $hgvsp;
    my $sample;
    
    for my $i (0..$#h) {
        if ($h[$i] eq "Hugo_Symbol") {
            $gene=$i;
        }
        if ($h[$i] eq "Variant_Type") {
            $vtype=$i;
        }
        if ($h[$i] eq "HGVSp") {
            $hgvsp=$i;
        }
        if ($h[$i] eq "Tumor_Sample_Barcode") {
            $sample=$i;
        }
        if ($h[$i] eq "Consequence") {
            $classification=$i;
        }
    }
    
    # Also check for other sample ID columns if Tumor_Sample_Barcode not found
    if (!defined $sample) {
        $sample = get_sample_id_column_index(\@h);
    }
    
    unless ($hgvsp) { 
        warn "Missing HGVSp column in $projname"; 
        next; 
    }
    
    my $file_filtered = 0;
    my $file_kept = 0;
    
    while (my $b=<F>) {
        next if ($b=~/^\#/);
        chomp $b;
        my @d=split(/\t/,$b);
        
        # Per-sample filtering
        if ($has_metadata && defined $sample && !$include_model) {
            my $sample_id = $d[$sample];
            $sample_id =~ s/^\s+|\s+$//g if defined $sample_id;
            
            if (defined $sample_id && exists $model_samples{$sample_id}) {
                $file_filtered++;
                $total_filtered++;
                next;  # Skip this row - it's a model sample
            }
        }
        
        my ($fsstart,$fslen);
        $d[$vtype]="SNP" if (defined $classification && $d[$classification]=~/inframe/i);
        
        if (!defined $classification || $d[$classification] !~/frameshift/i) {
            $fslen=0;
        }
        else {
            ($fsstart,$fslen)=$d[$hgvsp]=~/(\d+)/g;
        }
        unless ($fslen) { $fslen=0; }
        
        my $sample_val = defined $sample ? $d[$sample] : "";
        print "$projname\t".$d[$gene]."\t".$sample_val."\t".$d[$vtype]."\t".$d[$hgvsp]."\t".$fsstart."\t".$fslen."\n";
        $cnt{$d[$gene]}++;
        $file_kept++;
        $total_kept++;
    }
    
    close(F);
    
    if ($verbose && $file_filtered > 0) {
        print STDERR "Filtered $file_filtered model sample rows from $projname\n";
    }
}

foreach my $genek (keys %cnt) {
    print CNT $genek."\t".$cnt{$genek}."\n";
}
close(CNT);

if ($verbose) {
    print STDERR "Total: kept $total_kept rows, filtered $total_filtered model sample rows\n";
}
