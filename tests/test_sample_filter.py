#!/usr/bin/env python3
"""
Tests for sample_filter module and metadata-based per-sample filtering.
"""

import os
import sys
import tempfile
import unittest
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))

from sample_filter import (
    is_model_sample,
    load_sample_metadata,
    filter_mutation_rows,
    is_file_entirely_model,
    get_sample_id_column_index,
    _detect_delimiter,
    _discover_metadata_files,
)


class TestIsModelSample(unittest.TestCase):
    """Tests for the is_model_sample function."""
    
    def test_pdx_samples(self):
        """PDX samples should be identified as model samples."""
        self.assertTrue(is_model_sample("PDX"))
        self.assertTrue(is_model_sample("pdx"))
        self.assertTrue(is_model_sample("PDX Model"))
        self.assertTrue(is_model_sample("patient-derived xenograft"))
        self.assertTrue(is_model_sample("Patient Derived Xenograft"))
    
    def test_cell_line_samples(self):
        """Cell line samples should be identified as model samples."""
        self.assertTrue(is_model_sample("Cell Line"))
        self.assertTrue(is_model_sample("cell_line"))
        self.assertTrue(is_model_sample("cellline"))
        self.assertTrue(is_model_sample("CCLE"))
        self.assertTrue(is_model_sample("In Vitro"))
        self.assertTrue(is_model_sample("in-vitro"))
    
    def test_xenograft_samples(self):
        """Xenograft samples should be identified as model samples."""
        self.assertTrue(is_model_sample("Xenograft"))
        self.assertTrue(is_model_sample("xenograft model"))
    
    def test_model_keyword(self):
        """Samples with 'model' keyword should be identified."""
        self.assertTrue(is_model_sample("Model"))
        self.assertTrue(is_model_sample("tumor model"))
    
    def test_patient_samples(self):
        """Patient samples should not be identified as model samples."""
        self.assertFalse(is_model_sample("Patient"))
        self.assertFalse(is_model_sample("Primary Tumor"))
        self.assertFalse(is_model_sample("Metastatic"))
        self.assertFalse(is_model_sample("Normal"))
        self.assertFalse(is_model_sample("Tumor"))
    
    def test_empty_values(self):
        """Empty values should not be identified as model samples."""
        self.assertFalse(is_model_sample(""))
        self.assertFalse(is_model_sample(None))


class TestLoadSampleMetadata(unittest.TestCase):
    """Tests for the load_sample_metadata function."""
    
    @classmethod
    def setUpClass(cls):
        """Set up test fixtures path."""
        cls.fixtures_dir = Path(__file__).parent / 'fixtures'
    
    def test_load_metadata_from_file(self):
        """Test loading metadata from a specific file."""
        metadata_file = self.fixtures_dir / 'metadata.tsv'
        result = load_sample_metadata(str(metadata_file))
        
        # Check that we got the expected samples
        self.assertIn('SAMPLE001', result)
        self.assertIn('SAMPLE002', result)
        
        # Check model status
        self.assertFalse(result['SAMPLE001'])  # Patient
        self.assertTrue(result['SAMPLE002'])   # PDX
        self.assertFalse(result['SAMPLE003'])  # Patient
        self.assertTrue(result['SAMPLE004'])   # Cell Line
        self.assertFalse(result['SAMPLE005'])  # Primary Tumor
        self.assertTrue(result['SAMPLE006'])   # xenograft
        self.assertFalse(result['SAMPLE007'])  # Metastatic
        self.assertTrue(result['SAMPLE008'])   # patient-derived xenograft
    
    def test_load_metadata_from_directory(self):
        """Test auto-discovery of metadata files in a directory."""
        result = load_sample_metadata(str(self.fixtures_dir))
        
        # Should have found the metadata.tsv file
        self.assertGreater(len(result), 0)
        self.assertIn('SAMPLE001', result)
    
    def test_load_nonexistent_path(self):
        """Test loading from a nonexistent path returns empty dict."""
        result = load_sample_metadata('/nonexistent/path')
        self.assertEqual(result, {})


class TestFilterMutationRows(unittest.TestCase):
    """Tests for filtering mutation rows based on model samples."""
    
    @classmethod
    def setUpClass(cls):
        """Set up test fixtures path and load metadata."""
        cls.fixtures_dir = Path(__file__).parent / 'fixtures'
        metadata = load_sample_metadata(str(cls.fixtures_dir / 'metadata.tsv'))
        cls.model_samples = {s for s, is_model in metadata.items() if is_model}
    
    def test_filter_mixed_file(self):
        """Test filtering a mixed mutation file with patient and model samples."""
        mutation_file = self.fixtures_dir / 'data_mutations_mixed.tsv'
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as tmp:
            tmp_path = tmp.name
        
        try:
            rows_kept, rows_filtered = filter_mutation_rows(
                str(mutation_file),
                self.model_samples,
                tmp_path
            )
            
            # Mixed file has 5 rows: SAMPLE001, SAMPLE002, SAMPLE003, SAMPLE004, SAMPLE005
            # SAMPLE002 (PDX) and SAMPLE004 (Cell Line) are model samples
            self.assertEqual(rows_filtered, 2)
            self.assertEqual(rows_kept, 3)
            
            # Verify output file content
            with open(tmp_path, 'r') as f:
                lines = f.readlines()
            
            # Should have header + 3 data rows
            data_lines = [l for l in lines if not l.startswith('#')]
            self.assertEqual(len(data_lines), 4)  # 1 header + 3 data
            
            # Verify model samples are not in output
            content = ''.join(lines)
            self.assertNotIn('SAMPLE002', content)
            self.assertNotIn('SAMPLE004', content)
            
            # Verify patient samples are in output
            self.assertIn('SAMPLE001', content)
            self.assertIn('SAMPLE003', content)
            self.assertIn('SAMPLE005', content)
        
        finally:
            if os.path.exists(tmp_path):
                os.unlink(tmp_path)
    
    def test_filter_model_only_file(self):
        """Test filtering a file with only model samples."""
        mutation_file = self.fixtures_dir / 'data_mutations_model_only.tsv'
        
        rows_kept, rows_filtered = filter_mutation_rows(
            str(mutation_file),
            self.model_samples
        )
        
        # All 4 samples are model samples
        self.assertEqual(rows_filtered, 4)
        self.assertEqual(rows_kept, 0)
    
    def test_filter_patient_only_file(self):
        """Test filtering a file with only patient samples."""
        mutation_file = self.fixtures_dir / 'data_mutations_patient_only.tsv'
        
        rows_kept, rows_filtered = filter_mutation_rows(
            str(mutation_file),
            self.model_samples
        )
        
        # All 4 samples are patient samples
        self.assertEqual(rows_filtered, 0)
        self.assertEqual(rows_kept, 4)


class TestIsFileEntirelyModel(unittest.TestCase):
    """Tests for detecting files containing only model samples."""
    
    @classmethod
    def setUpClass(cls):
        """Set up test fixtures path and load metadata."""
        cls.fixtures_dir = Path(__file__).parent / 'fixtures'
        metadata = load_sample_metadata(str(cls.fixtures_dir / 'metadata.tsv'))
        cls.model_samples = {s for s, is_model in metadata.items() if is_model}
    
    def test_model_only_file(self):
        """File with only model samples should be detected."""
        mutation_file = self.fixtures_dir / 'data_mutations_model_only.tsv'
        result = is_file_entirely_model(str(mutation_file), self.model_samples)
        self.assertTrue(result)
    
    def test_mixed_file(self):
        """Mixed file should not be detected as entirely model."""
        mutation_file = self.fixtures_dir / 'data_mutations_mixed.tsv'
        result = is_file_entirely_model(str(mutation_file), self.model_samples)
        self.assertFalse(result)
    
    def test_patient_only_file(self):
        """Patient-only file should not be detected as entirely model."""
        mutation_file = self.fixtures_dir / 'data_mutations_patient_only.tsv'
        result = is_file_entirely_model(str(mutation_file), self.model_samples)
        self.assertFalse(result)


class TestColumnDetection(unittest.TestCase):
    """Tests for sample ID column detection."""
    
    def test_detect_tumor_sample_barcode(self):
        """Should detect Tumor_Sample_Barcode column."""
        headers = ['Hugo_Symbol', 'Variant_Type', 'Tumor_Sample_Barcode', 'HGVSp']
        idx = get_sample_id_column_index(headers)
        self.assertEqual(idx, 2)
    
    def test_detect_sample_id(self):
        """Should detect sample_id column."""
        headers = ['gene', 'sample_id', 'mutation']
        idx = get_sample_id_column_index(headers)
        self.assertEqual(idx, 1)
    
    def test_detect_sample_case_insensitive(self):
        """Should detect sample columns case-insensitively."""
        headers = ['Gene', 'SAMPLE', 'Mutation']
        idx = get_sample_id_column_index(headers)
        self.assertEqual(idx, 1)
    
    def test_no_sample_column(self):
        """Should return None when no sample column found."""
        headers = ['gene', 'mutation', 'consequence']
        idx = get_sample_id_column_index(headers)
        self.assertIsNone(idx)


class TestDelimiterDetection(unittest.TestCase):
    """Tests for file delimiter detection."""
    
    def test_detect_tsv(self):
        """Should detect TSV delimiter."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as tmp:
            tmp.write("col1\tcol2\tcol3\n")
            tmp.write("val1\tval2\tval3\n")
            tmp_path = tmp.name
        
        try:
            delimiter = _detect_delimiter(tmp_path)
            self.assertEqual(delimiter, '\t')
        finally:
            os.unlink(tmp_path)
    
    def test_detect_csv(self):
        """Should detect CSV delimiter."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as tmp:
            tmp.write("col1,col2,col3\n")
            tmp.write("val1,val2,val3\n")
            tmp_path = tmp.name
        
        try:
            delimiter = _detect_delimiter(tmp_path)
            self.assertEqual(delimiter, ',')
        finally:
            os.unlink(tmp_path)


class TestMetadataDiscovery(unittest.TestCase):
    """Tests for metadata file discovery."""
    
    def test_discover_metadata_tsv(self):
        """Should discover metadata.tsv file."""
        fixtures_dir = Path(__file__).parent / 'fixtures'
        files = _discover_metadata_files(str(fixtures_dir))
        
        # Should find metadata.tsv
        filenames = [os.path.basename(f) for f in files]
        self.assertIn('metadata.tsv', filenames)
    
    def test_discover_nonexistent_directory(self):
        """Should return empty list for nonexistent directory."""
        files = _discover_metadata_files('/nonexistent/path')
        self.assertEqual(files, [])


class TestIntegration(unittest.TestCase):
    """Integration tests for the complete filtering workflow."""
    
    @classmethod
    def setUpClass(cls):
        """Set up test fixtures."""
        cls.fixtures_dir = Path(__file__).parent / 'fixtures'
    
    def test_end_to_end_mixed_file_filtering(self):
        """Test complete workflow: load metadata, filter mixed file."""
        # Load metadata
        metadata = load_sample_metadata(str(self.fixtures_dir))
        model_samples = {s for s, is_model in metadata.items() if is_model}
        
        # Verify we loaded metadata correctly
        self.assertIn('SAMPLE002', model_samples)  # PDX
        self.assertIn('SAMPLE004', model_samples)  # Cell Line
        self.assertNotIn('SAMPLE001', model_samples)  # Patient
        
        # Filter mutation file
        mutation_file = self.fixtures_dir / 'data_mutations_mixed.tsv'
        rows_kept, rows_filtered = filter_mutation_rows(
            str(mutation_file),
            model_samples
        )
        
        # Verify filtering worked
        self.assertEqual(rows_filtered, 2)
        self.assertEqual(rows_kept, 3)
    
    def test_end_to_end_model_file_detection(self):
        """Test complete workflow: load metadata, detect model-only file."""
        # Load metadata
        metadata = load_sample_metadata(str(self.fixtures_dir))
        model_samples = {s for s, is_model in metadata.items() if is_model}
        
        # Check model-only file
        model_file = self.fixtures_dir / 'data_mutations_model_only.tsv'
        self.assertTrue(is_file_entirely_model(str(model_file), model_samples))
        
        # Check mixed file
        mixed_file = self.fixtures_dir / 'data_mutations_mixed.tsv'
        self.assertFalse(is_file_entirely_model(str(mixed_file), model_samples))


if __name__ == '__main__':
    unittest.main()
