#!/usr/bin/perl
# test_sample_filter.pl - Tests for sample_filter.pl functions

use strict;
use warnings;
use Test::More tests => 33;
use File::Basename;
use Cwd 'abs_path';

# Get the directory of this script and include the module
my $script_dir = dirname(abs_path($0));
my $repo_dir = dirname($script_dir);
require "$repo_dir/sample_filter.pl";

my $fixtures_dir = "$script_dir/fixtures";

# Test is_model_sample function
subtest 'is_model_sample - PDX samples' => sub {
    plan tests => 5;
    ok(is_model_sample("PDX"), "PDX detected as model");
    ok(is_model_sample("pdx"), "pdx (lowercase) detected as model");
    ok(is_model_sample("PDX Model"), "PDX Model detected as model");
    ok(is_model_sample("patient-derived xenograft"), "patient-derived xenograft detected as model");
    ok(is_model_sample("Patient Derived Xenograft"), "Patient Derived Xenograft detected as model");
};

subtest 'is_model_sample - Cell line samples' => sub {
    plan tests => 6;
    ok(is_model_sample("Cell Line"), "Cell Line detected as model");
    ok(is_model_sample("cell_line"), "cell_line detected as model");
    ok(is_model_sample("cellline"), "cellline detected as model");
    ok(is_model_sample("CCLE"), "CCLE detected as model");
    ok(is_model_sample("In Vitro"), "In Vitro detected as model");
    ok(is_model_sample("in-vitro"), "in-vitro detected as model");
};

subtest 'is_model_sample - Xenograft samples' => sub {
    plan tests => 2;
    ok(is_model_sample("Xenograft"), "Xenograft detected as model");
    ok(is_model_sample("xenograft model"), "xenograft model detected as model");
};

subtest 'is_model_sample - Model keyword' => sub {
    plan tests => 2;
    ok(is_model_sample("Model"), "Model detected as model");
    ok(is_model_sample("tumor model"), "tumor model detected as model");
};

subtest 'is_model_sample - Patient samples (not model)' => sub {
    plan tests => 5;
    ok(!is_model_sample("Patient"), "Patient not detected as model");
    ok(!is_model_sample("Primary Tumor"), "Primary Tumor not detected as model");
    ok(!is_model_sample("Metastatic"), "Metastatic not detected as model");
    ok(!is_model_sample("Normal"), "Normal not detected as model");
    ok(!is_model_sample("Tumor"), "Tumor not detected as model");
};

subtest 'is_model_sample - Empty values' => sub {
    plan tests => 2;
    ok(!is_model_sample(""), "Empty string not detected as model");
    ok(!is_model_sample(undef), "undef not detected as model");
};

# Test find_column_index
subtest 'find_column_index' => sub {
    plan tests => 3;
    my @headers = ('Hugo_Symbol', 'Variant_Type', 'Tumor_Sample_Barcode', 'HGVSp');
    my @sample_patterns = ('tumor_sample_barcode', 'sample_id', 'sample');
    
    my $idx = find_column_index(\@headers, \@sample_patterns);
    is($idx, 2, "Found Tumor_Sample_Barcode at index 2");
    
    my @headers2 = ('gene', 'sample_id', 'mutation');
    $idx = find_column_index(\@headers2, \@sample_patterns);
    is($idx, 1, "Found sample_id at index 1");
    
    my @headers3 = ('gene', 'mutation', 'consequence');
    $idx = find_column_index(\@headers3, \@sample_patterns);
    ok(!defined($idx), "Returns undef when no column found");
};

# Test detect_delimiter
subtest 'detect_delimiter' => sub {
    plan tests => 2;
    is(detect_delimiter("col1\tcol2\tcol3"), "\t", "Detected TSV delimiter");
    is(detect_delimiter("col1,col2,col3"), ",", "Detected CSV delimiter");
};

# Test discover_metadata_files
subtest 'discover_metadata_files' => sub {
    plan tests => 2;
    my @files = discover_metadata_files($fixtures_dir);
    ok(scalar(@files) > 0, "Found metadata files in fixtures dir");
    ok(grep { /metadata\.tsv$/ } @files, "Found metadata.tsv");
};

subtest 'discover_metadata_files - nonexistent dir' => sub {
    plan tests => 1;
    my @files = discover_metadata_files("/nonexistent/path");
    is(scalar(@files), 0, "Returns empty list for nonexistent directory");
};

# Test parse_metadata_file
subtest 'parse_metadata_file' => sub {
    plan tests => 8;
    my %result = parse_metadata_file("$fixtures_dir/metadata.tsv");
    
    ok(exists $result{'SAMPLE001'}, "SAMPLE001 exists in result");
    ok(exists $result{'SAMPLE002'}, "SAMPLE002 exists in result");
    
    is($result{'SAMPLE001'}, 0, "SAMPLE001 (Patient) is not model");
    is($result{'SAMPLE002'}, 1, "SAMPLE002 (PDX) is model");
    is($result{'SAMPLE003'}, 0, "SAMPLE003 (Patient) is not model");
    is($result{'SAMPLE004'}, 1, "SAMPLE004 (Cell Line) is model");
    is($result{'SAMPLE005'}, 0, "SAMPLE005 (Primary Tumor) is not model");
    is($result{'SAMPLE006'}, 1, "SAMPLE006 (xenograft) is model");
};

# Test load_sample_metadata from file
subtest 'load_sample_metadata - from file' => sub {
    plan tests => 2;
    my %result = load_sample_metadata("$fixtures_dir/metadata.tsv");
    ok(exists $result{'SAMPLE001'}, "SAMPLE001 loaded from file");
    ok(exists $result{'SAMPLE008'}, "SAMPLE008 loaded from file");
};

# Test load_sample_metadata from directory
subtest 'load_sample_metadata - from directory' => sub {
    plan tests => 2;
    my %result = load_sample_metadata($fixtures_dir);
    ok(scalar(keys %result) > 0, "Loaded samples from directory");
    ok(exists $result{'SAMPLE001'}, "Found SAMPLE001 via directory discovery");
};

# Test load_sample_metadata from nonexistent path
subtest 'load_sample_metadata - nonexistent path' => sub {
    plan tests => 1;
    my %result = load_sample_metadata("/nonexistent/path");
    is(scalar(keys %result), 0, "Returns empty hash for nonexistent path");
};

# Test get_model_samples
subtest 'get_model_samples' => sub {
    plan tests => 4;
    my %metadata = parse_metadata_file("$fixtures_dir/metadata.tsv");
    my %model_samples = get_model_samples(\%metadata);
    
    ok(exists $model_samples{'SAMPLE002'}, "SAMPLE002 (PDX) in model samples");
    ok(exists $model_samples{'SAMPLE004'}, "SAMPLE004 (Cell Line) in model samples");
    ok(!exists $model_samples{'SAMPLE001'}, "SAMPLE001 (Patient) not in model samples");
    ok(!exists $model_samples{'SAMPLE005'}, "SAMPLE005 (Primary Tumor) not in model samples");
};

# Test get_sample_id_column_index
subtest 'get_sample_id_column_index' => sub {
    plan tests => 3;
    my @headers1 = ('Hugo_Symbol', 'Variant_Type', 'Tumor_Sample_Barcode', 'HGVSp');
    is(get_sample_id_column_index(\@headers1), 2, "Found Tumor_Sample_Barcode");
    
    my @headers2 = ('gene', 'sample_id', 'mutation');
    is(get_sample_id_column_index(\@headers2), 1, "Found sample_id");
    
    my @headers3 = ('gene', 'mutation');
    ok(!defined(get_sample_id_column_index(\@headers3)), "Returns undef when not found");
};

# Integration tests
subtest 'Integration - Mixed file filtering' => sub {
    plan tests => 4;
    my %metadata = load_sample_metadata($fixtures_dir);
    my %model_samples = get_model_samples(\%metadata);
    
    ok(exists $model_samples{'SAMPLE002'}, "PDX sample detected");
    ok(exists $model_samples{'SAMPLE004'}, "Cell Line sample detected");
    ok(!exists $model_samples{'SAMPLE001'}, "Patient sample not in model set");
    ok(!exists $model_samples{'SAMPLE003'}, "Patient sample not in model set");
};

subtest 'Integration - Model only file detection' => sub {
    plan tests => 1;
    my %metadata = load_sample_metadata($fixtures_dir);
    my %model_samples = get_model_samples(\%metadata);
    
    # Check if all samples in model_only file are model samples
    my @model_file_samples = qw(SAMPLE002 SAMPLE004 SAMPLE006 SAMPLE008);
    my $all_model = 1;
    foreach my $sample (@model_file_samples) {
        $all_model = 0 unless exists $model_samples{$sample};
    }
    ok($all_model, "All samples in model-only file are detected as model samples");
};

subtest 'Integration - Patient only file' => sub {
    plan tests => 1;
    my %metadata = load_sample_metadata($fixtures_dir);
    my %model_samples = get_model_samples(\%metadata);
    
    # Check that no samples in patient_only file are model samples
    my @patient_file_samples = qw(SAMPLE001 SAMPLE003 SAMPLE005 SAMPLE007);
    my $none_model = 1;
    foreach my $sample (@patient_file_samples) {
        $none_model = 0 if exists $model_samples{$sample};
    }
    ok($none_model, "No samples in patient-only file are detected as model samples");
};

# Additional tests for edge cases
subtest 'Edge case - Column name case insensitivity' => sub {
    plan tests => 2;
    my @headers = ('Gene', 'SAMPLE', 'Mutation');
    my @patterns = ('sample');
    my $idx = find_column_index(\@headers, \@patterns);
    is($idx, 1, "Case-insensitive matching works");
    
    my @headers2 = ('Gene', 'TUMOR_SAMPLE_BARCODE', 'Mutation');
    my @sample_patterns = ('tumor_sample_barcode');
    $idx = find_column_index(\@headers2, \@sample_patterns);
    is($idx, 1, "Case-insensitive matching for sample barcode");
};

subtest 'Edge case - Whitespace handling' => sub {
    plan tests => 2;
    ok(is_model_sample("  PDX  "), "PDX with whitespace detected");
    ok(is_model_sample("\tCell Line\n"), "Cell Line with tabs/newlines detected");
};

subtest 'Edge case - Model type variations' => sub {
    plan tests => 4;
    ok(is_model_sample("cell-line"), "cell-line detected");
    ok(is_model_sample("invitro"), "invitro detected");
    ok(is_model_sample("CCLE cell line"), "CCLE cell line detected");
    ok(is_model_sample("PDX xenograft"), "PDX xenograft detected");
};

# Test discover_subdirectories
subtest 'discover_subdirectories' => sub {
    plan tests => 3;
    my $projects_dir = "$fixtures_dir/projects";
    my @subdirs = discover_subdirectories($projects_dir);
    
    ok(scalar(@subdirs) >= 2, "Found at least 2 subdirectories");
    
    # Check that we get [name, path] pairs
    my %found_dirs = map { $_->[0] => $_->[1] } @subdirs;
    ok(exists $found_dirs{'project_alpha'}, "Found project_alpha subdirectory");
    ok(exists $found_dirs{'project_beta'}, "Found project_beta subdirectory");
};

subtest 'discover_subdirectories - nonexistent dir' => sub {
    plan tests => 1;
    my @subdirs = discover_subdirectories("/nonexistent/path");
    is(scalar(@subdirs), 0, "Returns empty list for nonexistent directory");
};

# Test load_sample_metadata_from_subdirs
subtest 'load_sample_metadata_from_subdirs' => sub {
    plan tests => 6;
    my $projects_dir = "$fixtures_dir/projects";
    my %project_metadata = load_sample_metadata_from_subdirs($projects_dir);
    
    ok(exists $project_metadata{'project_alpha'}, "project_alpha loaded as project ID");
    ok(exists $project_metadata{'project_beta'}, "project_beta loaded as project ID");
    
    # Verify project_alpha samples
    my $alpha_samples = $project_metadata{'project_alpha'};
    is($alpha_samples->{'ALPHA001'}, 0, "ALPHA001 (Patient) is not model");
    is($alpha_samples->{'ALPHA002'}, 1, "ALPHA002 (PDX) is model");
    
    # Verify project_beta samples
    my $beta_samples = $project_metadata{'project_beta'};
    is($beta_samples->{'BETA001'}, 1, "BETA001 (Cell Line) is model");
    is($beta_samples->{'BETA002'}, 0, "BETA002 (Patient) is not model");
};

subtest 'load_sample_metadata_from_subdirs - nonexistent dir' => sub {
    plan tests => 1;
    my %project_metadata = load_sample_metadata_from_subdirs("/nonexistent/path");
    is(scalar(keys %project_metadata), 0, "Returns empty hash for nonexistent directory");
};

subtest 'load_sample_metadata_from_subdirs - empty projects skipped' => sub {
    plan tests => 1;
    my $projects_dir = "$fixtures_dir/projects";
    my %project_metadata = load_sample_metadata_from_subdirs($projects_dir);
    
    # empty_project should not be in results (no metadata files)
    ok(!exists $project_metadata{'empty_project'}, "Empty project not included in results");
};

# Test normalize_sample_barcode
subtest 'normalize_sample_barcode - Basic normalization' => sub {
    plan tests => 5;
    is(normalize_sample_barcode("SAMPLE001"), "SAMPLE001", "Simple barcode unchanged");
    is(normalize_sample_barcode("  SAMPLE001  "), "SAMPLE001", "Whitespace trimmed");
    is(normalize_sample_barcode(""), "", "Empty string returns empty");
    is(normalize_sample_barcode(undef), "", "undef returns empty");
    is(normalize_sample_barcode("sample-123"), "sample-123", "Regular barcode unchanged");
};

subtest 'normalize_sample_barcode - TCGA samples' => sub {
    plan tests => 4;
    # TCGA format: TCGA-XX-XXXX-XXX-XXX-XXXX-XX
    # First 15 chars (TCGA-XX-XXXX-XX) identify sample
    is(normalize_sample_barcode("TCGA-AA-1234-01A-11D-1234-01"), "TCGA-AA-1234-01", "Full TCGA barcode normalized");
    is(normalize_sample_barcode("TCGA-AA-1234-01"), "TCGA-AA-1234-01", "Short TCGA barcode unchanged");
    is(normalize_sample_barcode("TCGA-BB-5678-02A-22R-5678-02"), "TCGA-BB-5678-02", "Different TCGA barcode normalized");
    # Same participant, different portions should match
    my $bc1 = normalize_sample_barcode("TCGA-AA-1234-01A-11D-1234-01");
    my $bc2 = normalize_sample_barcode("TCGA-AA-1234-01B-22R-5678-02");
    is($bc1, $bc2, "Same sample different portions normalize to same value");
};

# Test is_duplicate_sample and record_sample_seen
subtest 'is_duplicate_sample - Basic detection' => sub {
    plan tests => 4;
    my %seen_samples;
    my %first_project;
    
    ok(!is_duplicate_sample("SAMPLE001", \%seen_samples, \%first_project), "First sample not duplicate");
    record_sample_seen("SAMPLE001", "project_a", \%seen_samples, \%first_project);
    ok(is_duplicate_sample("SAMPLE001", \%seen_samples, \%first_project), "Same sample is duplicate");
    ok(!is_duplicate_sample("SAMPLE002", \%seen_samples, \%first_project), "Different sample not duplicate");
    ok(!is_duplicate_sample("", \%seen_samples, \%first_project), "Empty barcode not duplicate");
};

subtest 'is_duplicate_sample - TCGA normalization' => sub {
    plan tests => 3;
    my %seen_samples;
    my %first_project;
    
    # Record first TCGA sample
    record_sample_seen("TCGA-AA-1234-01A-11D-1234-01", "project_a", \%seen_samples, \%first_project);
    
    # Same sample with different portion/analyte should be duplicate
    ok(is_duplicate_sample("TCGA-AA-1234-01B-22R-5678-02", \%seen_samples, \%first_project), 
       "Same TCGA sample different portion is duplicate");
    
    # Different sample should not be duplicate
    ok(!is_duplicate_sample("TCGA-BB-5678-02A-11D-1234-01", \%seen_samples, \%first_project),
       "Different TCGA sample not duplicate");
    
    # Different sample type (e.g., 01 vs 02) should not be duplicate
    ok(!is_duplicate_sample("TCGA-AA-1234-02A-11D-1234-01", \%seen_samples, \%first_project),
       "Same patient different sample type not duplicate");
};

# Test get_sample_first_project
subtest 'get_sample_first_project' => sub {
    plan tests => 4;
    my %seen_samples;
    my %first_project;
    
    record_sample_seen("SAMPLE001", "project_a", \%seen_samples, \%first_project);
    record_sample_seen("SAMPLE002", "project_b", \%seen_samples, \%first_project);
    
    is(get_sample_first_project("SAMPLE001", \%first_project), "project_a", "Returns correct first project for SAMPLE001");
    is(get_sample_first_project("SAMPLE002", \%first_project), "project_b", "Returns correct first project for SAMPLE002");
    ok(!defined(get_sample_first_project("SAMPLE003", \%first_project)), "Returns undef for unknown sample");
    ok(!defined(get_sample_first_project("", \%first_project)), "Returns undef for empty barcode");
};

# Integration test for duplicate detection across projects
subtest 'Integration - Duplicate detection across projects' => sub {
    plan tests => 6;
    my %seen_samples;
    my %first_project;
    
    # Simulate processing samples from multiple projects
    my @project_a_samples = qw(SAMPLE001 SAMPLE002 SAMPLE003);
    my @project_b_samples = qw(SAMPLE002 SAMPLE004 SAMPLE005);  # SAMPLE002 is duplicate
    my @project_c_samples = qw(SAMPLE001 SAMPLE005 SAMPLE006);  # SAMPLE001 and SAMPLE005 are duplicates
    
    # Process project_a
    foreach my $sample (@project_a_samples) {
        record_sample_seen($sample, "project_a", \%seen_samples, \%first_project) 
            unless is_duplicate_sample($sample, \%seen_samples, \%first_project);
    }
    
    # Check project_b duplicates
    ok(!is_duplicate_sample("SAMPLE004", \%seen_samples, \%first_project), "SAMPLE004 is new in project_b");
    ok(is_duplicate_sample("SAMPLE002", \%seen_samples, \%first_project), "SAMPLE002 is duplicate from project_a");
    
    # Record project_b samples
    foreach my $sample (@project_b_samples) {
        record_sample_seen($sample, "project_b", \%seen_samples, \%first_project)
            unless is_duplicate_sample($sample, \%seen_samples, \%first_project);
    }
    
    # Check project_c duplicates
    ok(is_duplicate_sample("SAMPLE001", \%seen_samples, \%first_project), "SAMPLE001 duplicate in project_c");
    ok(is_duplicate_sample("SAMPLE005", \%seen_samples, \%first_project), "SAMPLE005 duplicate in project_c");
    ok(!is_duplicate_sample("SAMPLE006", \%seen_samples, \%first_project), "SAMPLE006 is new in project_c");
    
    # Verify first project tracking
    is(get_sample_first_project("SAMPLE002", \%first_project), "project_a", "SAMPLE002 first seen in project_a");
};

print "\nAll tests completed!\n";
