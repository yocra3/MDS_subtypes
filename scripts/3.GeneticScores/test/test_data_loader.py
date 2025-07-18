"""
Tests unitarios para el módulo DataLoader.

Autor: Sistema automatizado
Fecha: 2025-07-04
"""

import unittest
import pandas as pd
import numpy as np
from unittest.mock import Mock, patch, MagicMock
import sys
from pathlib import Path

# Añadir utils al path
sys.path.append(str(Path(__file__).parent.parent / "utils"))

from data_loader import DataLoader


class TestDataLoader(unittest.TestCase):
    """Tests para la clase DataLoader."""
    
    def setUp(self):
        """Configurar tests."""
        self.test_data_paths = {
            'base_path': '/test/path',
            'clinical_data': 'test_clinical.tsv',
            'mutations_data': 'test_mutations.tsv',
            'cna_data': 'test_cna.tsv',
            'gene_embeddings': 'test_embeddings.json',
            'cv_folds': 'test_folds.rds'
        }
        self.data_loader = DataLoader(self.test_data_paths)
    
    def test_init(self):
        """Test inicialización."""
        self.assertEqual(str(self.data_loader.base_path), '/test/path')
        self.assertEqual(self.data_loader.data_paths, self.test_data_paths)
    
    def test_load_clinical_data_not_implemented(self):
        """Test que load_clinical_data lanza NotImplementedError."""
        with self.assertRaises(NotImplementedError):
            self.data_loader.load_clinical_data()
    
    def test_load_mutations_data_not_implemented(self):
        """Test que load_mutations_data lanza NotImplementedError."""
        with self.assertRaises(NotImplementedError):
            self.data_loader.load_mutations_data()
    
    def test_load_cna_data_not_implemented(self):
        """Test que load_cna_data lanza NotImplementedError."""
        with self.assertRaises(NotImplementedError):
            self.data_loader.load_cna_data()
    
    def test_load_gene_embeddings_not_implemented(self):
        """Test que load_gene_embeddings lanza NotImplementedError."""
        with self.assertRaises(NotImplementedError):
            self.data_loader.load_gene_embeddings()
    
    def test_load_cv_folds_not_implemented(self):
        """Test que load_cv_folds lanza NotImplementedError."""
        with self.assertRaises(NotImplementedError):
            self.data_loader.load_cv_folds()
    
    def test_validate_data_consistency(self):
        """Test validación de consistencia de datos."""
        test_data = {'clinical': pd.DataFrame(), 'mutations': pd.DataFrame()}
        result = self.data_loader.validate_data_consistency(test_data)
        self.assertTrue(result)


if __name__ == '__main__':
    unittest.main()
