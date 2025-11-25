#!/usr/bin/env python3
"""
Merge cBioPortal mutation files with per-sample filtering.

This script merges mutation data files from cBioPortal studies while filtering
out model-derived samples (PDX, cell lines, etc.) using metadata files.

Usage:
    python combine_mutations.py [options]
    
Options:
    --input-dir DIR       Directory containing mutation files (default: current dir)
    --output FILE         Output file path (default: stdout)
    --metadata FILE       Optional metadata file for sample filtering
    --include-model       Include model samples (disable filtering)
    --verbose             Print verbose output
"""

import argparse
import os
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from sample_filter import (
    is_model_sample,
    load_sample_metadata,
    get_sample_id_column_index,
    is_file_entirely_model,
    _detect_delimiter,
    _discover_metadata_files,
)

# Patterns for study names that indicate model data (for whole-file filtering)
MODEL_STUDY_PATTERNS = [
    r'\bccle\b',
    r'\bpdx\b',
    r'\bcellline\b',
    r'\bcell_line\b',
    r'\bxenograft\b',
    r'\btest\b',
]


def is_model_study_path(path: str) -> bool:
    """
    Check if a file path indicates a model study based on directory/file names.
    This preserves backwards compatibility with the original Perl script behavior.
    
    Args:
        path: Path to check.
    
    Returns:
        True if the path indicates a model study.
    """
    path_lower = path.lower()
    for pattern in MODEL_STUDY_PATTERNS:
        if re.search(pattern, path_lower):
            return True
    return False


def find_mutation_files(directory: str) -> List[str]:
    """
    Find all mutation data files in a directory.
    
    Args:
        directory: Directory to search.
    
    Returns:
        List of paths to mutation files.
    """
    mutation_files = []
    
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.startswith('data_mutations'):
                mutation_files.append(os.path.join(root, file))
    
    return mutation_files


def extract_study_info(file_path: str) -> Tuple[str, str, str]:
    """
    Extract study information from a file path.
    
    Args:
        file_path: Path to the mutation file.
    
    Returns:
        Tuple of (project_name, study, center).
    """
    parts = file_path.split('/')
    
    if len(parts) >= 3:
        proj_name = parts[-2]
        study_parts = proj_name.split('_')
        if len(study_parts) >= 2:
            study = study_parts[0]
            center = study_parts[1]
        else:
            study = proj_name
            center = ''
    else:
        proj_name = os.path.basename(os.path.dirname(file_path))
        study = proj_name
        center = ''
    
    return proj_name, study, center


def process_mutation_file(
    file_path: str,
    model_samples: Set[str],
    seen_uids: Set[str],
    include_model: bool = False,
    verbose: bool = False
) -> List[str]:
    """
    Process a single mutation file and return filtered output lines.
    
    Args:
        file_path: Path to the mutation file.
        model_samples: Set of sample IDs that are model samples.
        seen_uids: Set of already-processed study UIDs (modified in place).
        include_model: If True, include model samples.
        verbose: Print verbose output.
    
    Returns:
        List of output lines.
    """
    output_lines = []
    
    proj_name, study, center = extract_study_info(file_path)
    uid = study + center
    
    # Skip if already seen
    if uid in seen_uids:
        if verbose:
            print(f"Skipping duplicate study: {proj_name}", file=sys.stderr)
        return []
    seen_uids.add(uid)
    
    # Check for model study path (backwards compatibility)
    if not include_model and is_model_study_path(file_path):
        if verbose:
            print(f"Skipping model study (path): {proj_name}", file=sys.stderr)
        return []
    
    # Check if entire file is model samples
    if not include_model and model_samples:
        if is_file_entirely_model(file_path, model_samples):
            if verbose:
                print(f"Skipping model study (all samples): {proj_name}", file=sys.stderr)
            return []
    
    delimiter = _detect_delimiter(file_path)
    
    try:
        with open(file_path, 'r', encoding='utf-8', errors='replace') as f:
            # Skip comment lines and find header
            header_line = None
            for line in f:
                if not line.startswith('#'):
                    header_line = line.strip()
                    break
            
            if not header_line:
                if verbose:
                    print(f"Warning: Empty file {file_path}", file=sys.stderr)
                return []
            
            headers = header_line.split(delimiter)
            
            # Find required column indices
            gene_idx = None
            vtype_idx = None
            hgvsp_idx = None
            sample_idx = None
            classification_idx = None
            
            for i, header in enumerate(headers):
                header_clean = header.strip()
                if header_clean == "Hugo_Symbol":
                    gene_idx = i
                elif header_clean == "Variant_Type":
                    vtype_idx = i
                elif header_clean == "HGVSp":
                    hgvsp_idx = i
                elif header_clean == "Tumor_Sample_Barcode":
                    sample_idx = i
                elif header_clean == "Consequence":
                    classification_idx = i
            
            # Also check for other sample ID columns if Tumor_Sample_Barcode not found
            if sample_idx is None:
                sample_idx = get_sample_id_column_index(headers)
            
            if hgvsp_idx is None:
                if verbose:
                    print(f"Warning: Missing HGVSp column in {proj_name}", file=sys.stderr)
                return []
            
            # Process data rows
            rows_kept = 0
            rows_filtered = 0
            
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                
                fields = line.strip().split(delimiter)
                
                # Per-sample filtering
                if not include_model and sample_idx is not None and model_samples:
                    if len(fields) > sample_idx:
                        sample_id = fields[sample_idx].strip()
                        if sample_id in model_samples:
                            rows_filtered += 1
                            continue
                
                # Extract values
                gene = fields[gene_idx] if gene_idx is not None and len(fields) > gene_idx else ""
                sample = fields[sample_idx] if sample_idx is not None and len(fields) > sample_idx else ""
                vtype = fields[vtype_idx] if vtype_idx is not None and len(fields) > vtype_idx else ""
                hgvsp = fields[hgvsp_idx] if len(fields) > hgvsp_idx else ""
                classification = fields[classification_idx] if classification_idx is not None and len(fields) > classification_idx else ""
                
                # Variant type adjustment (matching Perl logic)
                if 'inframe' in classification.lower():
                    vtype = "SNP"
                
                # Calculate frameshift info
                fsstart = ""
                fslen = "0"
                
                if 'frameshift' not in classification.lower():
                    fslen = "0"
                else:
                    # Extract numbers from HGVSp
                    numbers = re.findall(r'\d+', hgvsp)
                    if len(numbers) >= 1:
                        fsstart = numbers[0]
                    if len(numbers) >= 2:
                        fslen = numbers[1]
                    else:
                        fslen = "0"
                
                if not fslen:
                    fslen = "0"
                
                output_line = f"{proj_name}\t{gene}\t{sample}\t{vtype}\t{hgvsp}\t{fsstart}\t{fslen}"
                output_lines.append(output_line)
                rows_kept += 1
            
            if verbose and rows_filtered > 0:
                print(f"Filtered {rows_filtered} model sample rows from {proj_name}", file=sys.stderr)
    
    except Exception as e:
        print(f"Error processing {file_path}: {e}", file=sys.stderr)
        return []
    
    return output_lines


def main():
    parser = argparse.ArgumentParser(
        description='Merge cBioPortal mutation files with per-sample filtering'
    )
    parser.add_argument(
        '--input-dir', '-i',
        default='.',
        help='Directory containing mutation files (default: current dir)'
    )
    parser.add_argument(
        '--output', '-o',
        default=None,
        help='Output file path (default: stdout)'
    )
    parser.add_argument(
        '--metadata', '-m',
        action='append',
        help='Metadata file for sample filtering (can be specified multiple times)'
    )
    parser.add_argument(
        '--include-model',
        action='store_true',
        help='Include model samples (disable filtering)'
    )
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Print verbose output'
    )
    
    args = parser.parse_args()
    
    # Load metadata
    model_samples: Set[str] = set()
    
    if not args.include_model:
        # Auto-discover metadata files in input directory
        sample_metadata = load_sample_metadata(args.input_dir, args.metadata)
        
        # Extract model samples
        for sample_id, is_model in sample_metadata.items():
            if is_model:
                model_samples.add(sample_id)
        
        if args.verbose:
            if model_samples:
                print(f"Loaded {len(model_samples)} model samples from metadata", file=sys.stderr)
            else:
                print("No metadata files found, using path-based filtering only", file=sys.stderr)
    
    # Find mutation files
    mutation_files = find_mutation_files(args.input_dir)
    
    if args.verbose:
        print(f"Found {len(mutation_files)} mutation files", file=sys.stderr)
    
    # Process files
    seen_uids: Set[str] = set()
    all_output_lines: List[str] = []
    gene_counts: Dict[str, int] = {}
    
    for mutation_file in sorted(mutation_files):
        output_lines = process_mutation_file(
            mutation_file,
            model_samples,
            seen_uids,
            args.include_model,
            args.verbose
        )
        
        for line in output_lines:
            all_output_lines.append(line)
            # Count genes
            parts = line.split('\t')
            if len(parts) >= 2:
                gene = parts[1]
                gene_counts[gene] = gene_counts.get(gene, 0) + 1
    
    # Write output
    if args.output:
        with open(args.output, 'w', encoding='utf-8') as out:
            for line in all_output_lines:
                out.write(line + '\n')
        
        # Write gene counts
        cnt_file = os.path.splitext(args.output)[0] + '.cnt'
        with open(cnt_file, 'w', encoding='utf-8') as cnt:
            for gene, count in gene_counts.items():
                cnt.write(f"{gene}\t{count}\n")
    else:
        for line in all_output_lines:
            print(line)
    
    if args.verbose:
        print(f"Processed {len(all_output_lines)} mutation records", file=sys.stderr)


if __name__ == '__main__':
    main()
