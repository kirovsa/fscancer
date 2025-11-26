#!/usr/bin/perl
# combineCBIO.fs.filter.pl - Enhanced mutation file merger with per-sample filtering
#
# This script extends combineCBIO.fs.pl with metadata-based per-sample filtering.
# It filters out model-derived samples (PDX, cell lines, etc.) using metadata files.
#
# Usage:
#   perl combineCBIO.fs.filter.pl [--metadata FILE] [--verbose] [--stats FILE] > output.txt
#
# If metadata files are found in the current directory or specified via --metadata,
# individual mutation rows belonging to model samples will be filtered out.
# Otherwise, falls back to path-based filtering only (original behavior).
#
# Output files:
#   - Standard output: filtered mutation data
#   - cnt: gene counts
#   - filter_stats.tsv (or --stats FILE): per-project filtering statistics

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use File::Find;
use Cwd 'abs_path';

# Include sample_filter functions from same directory as this script
my $script_dir = dirname(abs_path($0));
require "$script_dir/sample_filter.pl";

# Parse command line arguments
my $metadata_file = '';
my $verbose = 0;
my $include_model = 0;
my $include_duplicates = 0;
my $stats_file = 'filter_stats.tsv';

GetOptions(
    'metadata=s'       => \$metadata_file,
    'verbose'          => \$verbose,
    'include-model'    => \$include_model,
    'include-duplicates' => \$include_duplicates,
    'stats=s'          => \$stats_file,
) or die "Usage: $0 [--metadata FILE] [--verbose] [--include-model] [--include-duplicates] [--stats FILE]\n";

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
# Using File::Find for safer file discovery
my @mutation_files;
find(sub {
    return unless -f $_;
    return unless /^data_mutations/;
    my $path = $File::Find::name;
    # Skip model-related paths (backwards compatibility)
    return if $path =~ /ccle/i;
    return if $path =~ /pdx/i;
    return if $path =~ /cellline/i;
    return if $path =~ /cell_line/i;
    return if $path =~ /xenograft/i;
    return if $path =~ /test/i;
    push @mutation_files, $path;
}, '.');

open(CNT,">cnt");

my %seen;
my %cnt;
my $total_filtered = 0;
my $total_kept = 0;
my $total_duplicates = 0;

# Global sample tracking for duplicate detection
my %seen_samples;      # normalized_barcode => 1
my %first_project;     # normalized_barcode => first_project_name

# Per-project statistics
my %project_stats;  # project_name => { accepted => N, rejected => N, duplicates => N, samples_accepted => {}, samples_rejected => {}, samples_duplicate => {} }

foreach my $proj (@mutation_files) {
    my ($dor,$projname,$file)=split(/\//,$proj);
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
    
    # Initialize project stats
    $project_stats{$projname} = {
        accepted => 0,
        rejected => 0,
        duplicates => 0,
        samples_accepted => {},
        samples_rejected => {},
        samples_duplicate => {}
    };
    
    my $file_filtered = 0;
    my $file_kept = 0;
    my $file_duplicates = 0;
    
    while (my $b=<F>) {
        next if ($b=~/^\#/);
        chomp $b;
        my @d=split(/\t/,$b);
        
        my $sample_id = defined $sample ? $d[$sample] : "";
        $sample_id =~ s/^\s+|\s+$//g if $sample_id;
        
        # Per-sample model filtering
        if ($has_metadata && defined $sample && !$include_model) {
            if (defined $sample_id && $sample_id ne '' && exists $model_samples{$sample_id}) {
                $file_filtered++;
                $total_filtered++;
                $project_stats{$projname}{rejected}++;
                $project_stats{$projname}{samples_rejected}{$sample_id} = 1;
                next;  # Skip this row - it's a model sample
            }
        }
        
        # Duplicate sample detection (across projects)
        if (defined $sample && !$include_duplicates && $sample_id ne '') {
            if (is_duplicate_sample($sample_id, \%seen_samples, \%first_project)) {
                $file_duplicates++;
                $total_duplicates++;
                $project_stats{$projname}{duplicates}++;
                $project_stats{$projname}{samples_duplicate}{$sample_id} = 1;
                next;  # Skip this row - duplicate sample
            }
            # Record sample as seen
            record_sample_seen($sample_id, $projname, \%seen_samples, \%first_project);
        }
        
        my $fsstart='',
	my $fslen='';
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
        $project_stats{$projname}{accepted}++;
        $project_stats{$projname}{samples_accepted}{$sample_val} = 1 if $sample_val;
    }
    
    close(F);
    
    if ($verbose && $file_filtered > 0) {
        print STDERR "Filtered $file_filtered model sample rows from $projname\n";
    }
    if ($verbose && $file_duplicates > 0) {
        print STDERR "Skipped $file_duplicates duplicate sample rows from $projname\n";
    }
}

foreach my $genek (keys %cnt) {
    print CNT $genek."\t".$cnt{$genek}."\n";
}
close(CNT);

# Write statistics file
open(my $stats_fh, '>', $stats_file) or warn "Could not open stats file $stats_file: $!";
if ($stats_fh) {
    print $stats_fh "# Filtering Statistics Report\n";
    print $stats_fh "# Generated by combineCBIO.fs.filter.pl\n";
    print $stats_fh "#\n";
    print $stats_fh "project\trows_accepted\trows_rejected\trows_duplicate\tsamples_accepted\tsamples_rejected\tsamples_duplicate\n";
    
    my $total_projects = 0;
    my $total_samples_accepted = 0;
    my $total_samples_rejected = 0;
    my $total_samples_duplicate = 0;
    
    foreach my $projname (sort keys %project_stats) {
        my $stats = $project_stats{$projname};
        my $samples_accepted_count = scalar(keys %{$stats->{samples_accepted}});
        my $samples_rejected_count = scalar(keys %{$stats->{samples_rejected}});
        my $samples_duplicate_count = scalar(keys %{$stats->{samples_duplicate}});
        
        print $stats_fh "$projname\t$stats->{accepted}\t$stats->{rejected}\t$stats->{duplicates}\t$samples_accepted_count\t$samples_rejected_count\t$samples_duplicate_count\n";
        
        $total_projects++;
        $total_samples_accepted += $samples_accepted_count;
        $total_samples_rejected += $samples_rejected_count;
        $total_samples_duplicate += $samples_duplicate_count;
    }
    
    print $stats_fh "#\n";
    print $stats_fh "# Summary\n";
    print $stats_fh "# Total projects processed: $total_projects\n";
    print $stats_fh "# Total rows accepted: $total_kept\n";
    print $stats_fh "# Total rows rejected (model): $total_filtered\n";
    print $stats_fh "# Total rows skipped (duplicate): $total_duplicates\n";
    print $stats_fh "# Total unique samples accepted: $total_samples_accepted\n";
    print $stats_fh "# Total unique samples rejected (model): $total_samples_rejected\n";
    print $stats_fh "# Total unique samples skipped (duplicate): $total_samples_duplicate\n";
    
    close($stats_fh);
    
    if ($verbose) {
        print STDERR "\nStatistics written to $stats_file\n";
    }
}

if ($verbose) {
    print STDERR "\n=== Summary ===\n";
    print STDERR "Total projects processed: " . scalar(keys %project_stats) . "\n";
    print STDERR "Total rows kept: $total_kept\n";
    print STDERR "Total rows filtered (model): $total_filtered\n";
    print STDERR "Total rows skipped (duplicate): $total_duplicates\n";
}
