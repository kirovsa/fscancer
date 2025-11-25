"""
Sample filtering module for mutation data.

This module provides functionality to filter out model-derived samples
(PDX, cell lines, etc.) from mutation files using metadata files.
"""

import csv
import os
import re
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple, Union

# Patterns for detecting metadata filenames
METADATA_FILENAME_PATTERNS = [
    'metadata.tsv',
    'samples.tsv',
    'clinical_sample.tsv',
    'sample_metadata.csv',
    'sample_annotations.tsv',
    'clinical.tsv',
]

# Patterns for files containing "sample" or "metadata" in name
METADATA_GLOB_PATTERNS = [
    '*sample*.tsv',
    '*sample*.csv',
    '*metadata*.tsv',
    '*metadata*.csv',
    '*clinical*.tsv',
    '*clinical*.csv',
]

# Common sample ID column names (case-insensitive)
SAMPLE_ID_COLUMNS = [
    'sample_id',
    'sample',
    'tumor_sample_barcode',
    'sample_barcode',
    'samplebarcode',
    'sampleid',
    'tumor_sample_id',
]

# Common sample type column names (case-insensitive)
SAMPLE_TYPE_COLUMNS = [
    'sample_type',
    'sampletype',
    'sample_type_detail',
    'model',
    'is_model',
    'sample_class',
    'sampleclass',
]

# Model type patterns (case-insensitive)
MODEL_TYPE_PATTERNS = [
    r'\bpdx\b',
    r'patient[- ]?derived[- ]?xenograft',
    r'\bxenograft\b',
    r'cell[- _]?line',
    r'\bcellline\b',
    r'in[- ]?vitro',
    r'\bmodel\b',
    r'\bccle\b',
]


def is_model_sample(sample_type_value: str) -> bool:
    """
    Determine if a sample type value indicates a model sample.
    
    Args:
        sample_type_value: The sample type string to check.
    
    Returns:
        True if the sample is a model (PDX, cell line, etc.), False otherwise.
    """
    if not sample_type_value:
        return False
    
    value_lower = str(sample_type_value).lower().strip()
    
    for pattern in MODEL_TYPE_PATTERNS:
        if re.search(pattern, value_lower, re.IGNORECASE):
            return True
    
    return False


def _find_column_index(headers: List[str], column_patterns: List[str]) -> Optional[int]:
    """
    Find the index of a column matching any of the given patterns.
    
    Args:
        headers: List of column headers.
        column_patterns: List of column name patterns to match (case-insensitive).
    
    Returns:
        The index of the matching column, or None if not found.
    """
    headers_lower = [h.lower().strip() for h in headers]
    
    for pattern in column_patterns:
        pattern_lower = pattern.lower()
        for i, header in enumerate(headers_lower):
            if header == pattern_lower:
                return i
    
    return None


def _detect_delimiter(file_path: str) -> str:
    """
    Detect the delimiter used in a file (TSV or CSV).
    
    Args:
        file_path: Path to the file.
    
    Returns:
        The detected delimiter ('\\t' for TSV, ',' for CSV).
    """
    with open(file_path, 'r', encoding='utf-8', errors='replace') as f:
        first_line = ''
        for line in f:
            if not line.startswith('#'):
                first_line = line
                break
        
        if not first_line:
            return '\t'
        
        tab_count = first_line.count('\t')
        comma_count = first_line.count(',')
        
        return '\t' if tab_count >= comma_count else ','


def _discover_metadata_files(directory: str) -> List[str]:
    """
    Discover metadata files in a directory.
    
    Args:
        directory: Directory to search for metadata files.
    
    Returns:
        List of paths to discovered metadata files.
    """
    found_files = []
    dir_path = Path(directory)
    
    if not dir_path.exists():
        return found_files
    
    # Check for exact filename matches
    for filename in METADATA_FILENAME_PATTERNS:
        file_path = dir_path / filename
        if file_path.exists():
            found_files.append(str(file_path))
    
    # Check for glob patterns
    for pattern in METADATA_GLOB_PATTERNS:
        for file_path in dir_path.glob(pattern):
            if str(file_path) not in found_files:
                found_files.append(str(file_path))
    
    return found_files


def load_sample_metadata(
    path_or_dir: str,
    additional_paths: Optional[List[str]] = None
) -> Dict[str, bool]:
    """
    Load sample metadata from a file or directory and return a mapping
    of sample barcodes to whether they are model samples.
    
    Args:
        path_or_dir: Path to a metadata file or directory containing metadata files.
        additional_paths: Optional list of additional metadata file paths to check.
    
    Returns:
        Dictionary mapping sample barcode -> is_model (bool).
        Returns empty dict if no metadata is found.
    """
    sample_is_model: Dict[str, bool] = {}
    metadata_files: List[str] = []
    
    path = Path(path_or_dir)
    
    if path.is_file():
        metadata_files.append(str(path))
    elif path.is_dir():
        metadata_files.extend(_discover_metadata_files(str(path)))
    
    if additional_paths:
        for p in additional_paths:
            if Path(p).is_file() and p not in metadata_files:
                metadata_files.append(p)
    
    for metadata_file in metadata_files:
        try:
            file_mapping = _parse_metadata_file(metadata_file)
            sample_is_model.update(file_mapping)
        except Exception as e:
            # Log warning but continue with other files
            print(f"Warning: Could not parse metadata file {metadata_file}: {e}")
            continue
    
    return sample_is_model


def _parse_metadata_file(file_path: str) -> Dict[str, bool]:
    """
    Parse a single metadata file and extract sample ID to model status mapping.
    
    Args:
        file_path: Path to the metadata file.
    
    Returns:
        Dictionary mapping sample barcode -> is_model (bool).
    """
    sample_is_model: Dict[str, bool] = {}
    delimiter = _detect_delimiter(file_path)
    
    with open(file_path, 'r', encoding='utf-8', errors='replace') as f:
        # Skip comment lines
        header_line = None
        for line in f:
            if not line.startswith('#'):
                header_line = line.strip()
                break
        
        if not header_line:
            return sample_is_model
        
        headers = header_line.split(delimiter)
        
        # Find sample ID column
        sample_id_idx = _find_column_index(headers, SAMPLE_ID_COLUMNS)
        if sample_id_idx is None:
            return sample_is_model
        
        # Find sample type column
        sample_type_idx = _find_column_index(headers, SAMPLE_TYPE_COLUMNS)
        if sample_type_idx is None:
            return sample_is_model
        
        # Parse data rows
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            
            fields = line.strip().split(delimiter)
            
            if len(fields) <= max(sample_id_idx, sample_type_idx):
                continue
            
            sample_id = fields[sample_id_idx].strip()
            sample_type = fields[sample_type_idx].strip()
            
            if sample_id:
                sample_is_model[sample_id] = is_model_sample(sample_type)
    
    return sample_is_model


def get_sample_id_column_index(headers: List[str]) -> Optional[int]:
    """
    Find the sample ID column index in mutation file headers.
    
    Args:
        headers: List of column headers from a mutation file.
    
    Returns:
        The index of the sample ID column, or None if not found.
    """
    return _find_column_index(headers, SAMPLE_ID_COLUMNS)


def filter_mutation_rows(
    mutation_file: str,
    model_samples: Set[str],
    output_file: Optional[str] = None
) -> Tuple[int, int]:
    """
    Filter mutation rows, removing those belonging to model samples.
    
    Args:
        mutation_file: Path to the mutation file.
        model_samples: Set of sample IDs that are model samples.
        output_file: Optional path for filtered output. If None, returns counts only.
    
    Returns:
        Tuple of (rows_kept, rows_filtered).
    """
    delimiter = _detect_delimiter(mutation_file)
    rows_kept = 0
    rows_filtered = 0
    
    output_lines = []
    
    with open(mutation_file, 'r', encoding='utf-8', errors='replace') as f:
        # Collect and preserve comment lines and find header
        comment_lines = []
        header_line = None
        
        for line in f:
            if line.startswith('#'):
                comment_lines.append(line)
            else:
                header_line = line
                break
        
        if not header_line:
            return (0, 0)
        
        headers = header_line.strip().split(delimiter)
        sample_id_idx = get_sample_id_column_index(headers)
        
        if sample_id_idx is None:
            # Can't filter without sample ID column - keep all rows
            if output_file:
                with open(output_file, 'w', encoding='utf-8') as out:
                    for comment in comment_lines:
                        out.write(comment)
                    out.write(header_line)
                    for line in f:
                        out.write(line)
                        rows_kept += 1
            return (rows_kept, 0)
        
        # Add comment lines and header to output
        output_lines.extend(comment_lines)
        output_lines.append(header_line)
        
        # Process data rows
        for line in f:
            if line.startswith('#') or not line.strip():
                output_lines.append(line)
                continue
            
            fields = line.strip().split(delimiter)
            
            if len(fields) <= sample_id_idx:
                output_lines.append(line)
                rows_kept += 1
                continue
            
            sample_id = fields[sample_id_idx].strip()
            
            if sample_id in model_samples:
                rows_filtered += 1
            else:
                output_lines.append(line)
                rows_kept += 1
    
    if output_file:
        with open(output_file, 'w', encoding='utf-8') as out:
            for line in output_lines:
                out.write(line if line.endswith('\n') else line + '\n')
    
    return (rows_kept, rows_filtered)


def is_file_entirely_model(
    mutation_file: str,
    model_samples: Set[str]
) -> bool:
    """
    Check if a mutation file contains only model samples.
    
    Args:
        mutation_file: Path to the mutation file.
        model_samples: Set of sample IDs that are model samples.
    
    Returns:
        True if all samples in the file are model samples, False otherwise.
    """
    delimiter = _detect_delimiter(mutation_file)
    
    with open(mutation_file, 'r', encoding='utf-8', errors='replace') as f:
        # Skip comment lines and find header
        header_line = None
        for line in f:
            if not line.startswith('#'):
                header_line = line
                break
        
        if not header_line:
            return False
        
        headers = header_line.strip().split(delimiter)
        sample_id_idx = get_sample_id_column_index(headers)
        
        if sample_id_idx is None:
            return False
        
        samples_in_file: Set[str] = set()
        
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            
            fields = line.strip().split(delimiter)
            
            if len(fields) <= sample_id_idx:
                continue
            
            sample_id = fields[sample_id_idx].strip()
            if sample_id:
                samples_in_file.add(sample_id)
        
        if not samples_in_file:
            return False
        
        # File is entirely model if all samples are in model_samples
        return samples_in_file.issubset(model_samples)
