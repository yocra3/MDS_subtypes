"""
Tests de integración para el pipeline unificado.

Autor: Sistema automatizado
Fecha: 2025-07-04
"""

import unittest
import tempfile
import shutil
from pathlib import Path
import yaml
import json
import sys
import pandas as pd

# Añadir el directorio padre al path para importar el módulo principal
sys.path.append(str(Path(__file__).parent.parent))

from prepare_data_unified import UnifiedDataPipeline


class TestUnifiedDataPipeline(unittest.TestCase):
    """Tests de integración para el pipeline unificado."""
    
    def setUp(self):
        """Configurar tests."""
        # Crear directorio temporal
        self.temp_dir = tempfile.mkdtemp()
        self.temp_path = Path(self.temp_dir)
        
        # Crear configuración de test
        self.test_config = {
            'data_paths': {
                'base_path': str(self.temp_path),
                'clinical_data': 'test_clinical.tsv',
                'mutations_data': 'test_mutations.tsv',
                'cv_folds': 'test_folds.rds'
            },
            'clinical_processing': {
                'variables_to_include': ['Age', 'Sex', 'OS', 'OS_STATUS'],
                'categorical_variables': ['Sex'],
                'continuous_variables': ['Age'],
                'scaling_method': 'standard',
                'encoding_method': 'onehot'
            },
            'data_generation': {
                'output_formats': ['table']
            }
        }
        
        # Guardar configuración
        self.config_path = self.temp_path / 'test_config.yaml'
        with open(self.config_path, 'w') as f:
            yaml.dump(self.test_config, f)
        
        # Crear directorio de salida
        self.output_dir = self.temp_path / 'output'
    
    def tearDown(self):
        """Limpiar después de tests."""
        shutil.rmtree(self.temp_dir)
    
    def test_pipeline_init(self):
        """Test inicialización del pipeline."""
        pipeline = UnifiedDataPipeline(str(self.config_path), str(self.output_dir))
        
        self.assertEqual(pipeline.config_path, self.config_path)
        self.assertEqual(pipeline.output_dir, self.output_dir)
        self.assertEqual(pipeline.config, self.test_config)
        
        # Verificar que el directorio de salida se creó
        self.assertTrue(self.output_dir.exists())
    
    def test_load_config_yaml(self):
        """Test carga de configuración YAML."""
        pipeline = UnifiedDataPipeline(str(self.config_path), str(self.output_dir))
        
        loaded_config = pipeline._load_config()
        self.assertEqual(loaded_config, self.test_config)
    
    def test_load_config_json(self):
        """Test carga de configuración JSON."""
        # Crear configuración JSON
        json_config_path = self.temp_path / 'test_config.json'
        with open(json_config_path, 'w') as f:
            json.dump(self.test_config, f)
        
        pipeline = UnifiedDataPipeline(str(json_config_path), str(self.output_dir))
        
        loaded_config = pipeline._load_config()
        self.assertEqual(loaded_config, self.test_config)
    
    def test_load_config_invalid_format(self):
        """Test error con formato de configuración inválido."""
        # Crear archivo con extensión inválida
        invalid_config_path = self.temp_path / 'test_config.txt'
        invalid_config_path.write_text('invalid config')
        
        with self.assertRaises(ValueError):
            UnifiedDataPipeline(str(invalid_config_path), str(self.output_dir))
    
    def test_run_pipeline_not_implemented(self):
        """Test que run_pipeline falla porque los módulos no están implementados."""
        pipeline = UnifiedDataPipeline(str(self.config_path), str(self.output_dir))
        
        # Debería fallar porque los módulos aún no están implementados
        with self.assertRaises(NotImplementedError):
            pipeline.run_pipeline()


class TestConfigurationFiles(unittest.TestCase):
    """Tests para validar archivos de configuración."""
    
    def setUp(self):
        """Configurar tests."""
        self.config_dir = Path(__file__).parent.parent / 'configs'
    
    def test_default_config_exists(self):
        """Test que existe configuración por defecto."""
        default_config = self.config_dir / 'default.yaml'
        self.assertTrue(default_config.exists())
    
    def test_gnn_config_exists(self):
        """Test que existe configuración para GNN."""
        gnn_config = self.config_dir / 'gnn_config.yaml'
        self.assertTrue(gnn_config.exists())
    
    def test_table_config_exists(self):
        """Test que existe configuración para tablas."""
        table_config = self.config_dir / 'table_config.yaml'
        self.assertTrue(table_config.exists())
    
    def test_default_config_valid_yaml(self):
        """Test que la configuración por defecto es YAML válido."""
        default_config = self.config_dir / 'default.yaml'
        
        with open(default_config, 'r') as f:
            config = yaml.safe_load(f)
        
        # Verificar secciones principales (actualizadas para configuración simplificada)
        self.assertIn('data_paths', config)
        self.assertIn('variable_processing', config)  # Cambió de 'clinical_processing'
        self.assertIn('data_generation', config)
    
    def test_gnn_config_valid_yaml(self):
        """Test que la configuración GNN es YAML válido."""
        gnn_config = self.config_dir / 'gnn_config.yaml'
        
        with open(gnn_config, 'r') as f:
            config = yaml.safe_load(f)
        
        # Verificar que tiene configuración específica para grafos
        self.assertIn('graph_format', config['data_generation'])
        self.assertEqual(config['data_generation']['output_formats'], ['graph'])
    
    def test_table_config_valid_yaml(self):
        """Test que la configuración de tablas es YAML válido."""
        table_config = self.config_dir / 'table_config.yaml'
        
        with open(table_config, 'r') as f:
            config = yaml.safe_load(f)
        
        # Verificar que tiene configuración específica para tablas
        self.assertIn('table_format', config['data_generation'])
        self.assertEqual(config['data_generation']['output_formats'], ['table'])


if __name__ == '__main__':
    unittest.main()
