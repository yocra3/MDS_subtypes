"""
Tests unitarios para gnn_cv_functions.py
=======================================

Tests para validar las funciones de cross-validation para modelos GNN.

Ejecutar:
python -m pytest scripts/3.GeneticScores/tests/test_gnn_cv_functions.py -v

O desde dentro del contenedor Docker:
cd /workspace && python scripts/3.GeneticScores/tests/test_gnn_cv_functions.py
"""

import os
import sys
import tempfile
import shutil
import yaml
from pathlib import Path
import numpy as np
import torch
import pytest
from torch_geometric.data import HeteroData

# Añadir path para imports
sys.path.append('scripts/3.GeneticScores/utils')
from gnn_cv_functions import GNNCrossValidator, run_gnn_cv_pipeline


class TestGNNCrossValidator:
    """Tests para la clase GNNCrossValidator."""
    
    @pytest.fixture
    def temp_dirs(self):
        """Fixture para crear directorios temporales."""
        temp_dir = tempfile.mkdtemp()
        yield temp_dir
        shutil.rmtree(temp_dir)
    
    @pytest.fixture
    def sample_graph_data(self, temp_dirs):
        """Fixture para crear datos de grafos de prueba."""
        def create_sample_hetero_graph(patient_idx: int) -> HeteroData:
            """Crea un grafo HeteroData de muestra."""
            graph = HeteroData()
            
            # Patient node
            graph['patient'].x = torch.randn(1, 13)  # 13 features
            graph['patient'].y = torch.tensor([[np.random.exponential(2), np.random.binomial(1, 0.3)]], dtype=torch.float)
            graph['patient'].ID = patient_idx
            
            # Gene nodes (4 genes)
            n_genes = 4
            graph['gene'].x = torch.randn(n_genes, 16)  # 16-dim embeddings
            graph['gene'].weights = torch.rand(n_genes)  # VAF weights
            
            # Edges (patient -> genes)
            graph['gene', 'patient'].edge_index = torch.tensor([
                list(range(n_genes)),  # src: genes
                [0] * n_genes          # dst: patient
            ], dtype=torch.long)
            
            return graph
        
        # Crear datos de grafos para 3 folds con 10 pacientes cada uno
        graph_data = {}
        for fold in range(3):
            fold_graphs = []
            for patient_idx in range(10):
                graph = create_sample_hetero_graph(patient_idx + fold * 10)
                fold_graphs.append(graph)
            graph_data[f'fold_{fold}'] = fold_graphs
        
        # Guardar archivo de prueba
        graph_file = Path(temp_dirs) / "test_graphs.pt"
        torch.save(graph_data, graph_file)
        
        return str(graph_file), graph_data
    
    def test_initialization(self, temp_dirs):
        """Test inicialización del validador."""
        graph_path = Path(temp_dirs) / "dummy_graphs.pt"
        results_dir = Path(temp_dirs) / "results"
        
        cv_validator = GNNCrossValidator(
            graph_data_path=str(graph_path),
            results_dir=str(results_dir)
        )
        
        assert cv_validator.graph_data_path == str(graph_path)
        assert cv_validator.results_dir == results_dir
        assert results_dir.exists()
    
    def test_load_graph_folds(self, sample_graph_data):
        """Test carga de folds de grafos."""
        graph_file, expected_data = sample_graph_data
        
        cv_validator = GNNCrossValidator(graph_data_path=graph_file)
        loaded_data = cv_validator.load_graph_folds()
        
        assert cv_validator.n_folds == 3
        assert len(loaded_data) == 3
        assert all(f'fold_{i}' in loaded_data for i in range(3))
        assert all(len(loaded_data[f'fold_{i}']) == 10 for i in range(3))
    
    def test_create_fold_dataloaders(self, sample_graph_data):
        """Test creación de DataLoaders por fold."""
        graph_file, _ = sample_graph_data
        
        cv_validator = GNNCrossValidator(graph_data_path=graph_file)
        cv_validator.load_graph_folds()
        
        # Test fold 0 como test
        train_loader, test_loader = cv_validator.create_fold_dataloaders(
            test_fold=0,
            batch_size_train=4,
            batch_size_test=8
        )
        
        # Verificar tamaños
        assert test_loader.batch_size == 8
        assert train_loader.batch_size == 4
        
        # Verificar que test tiene datos del fold 0 y train del resto
        test_batch = next(iter(test_loader))
        train_batch = next(iter(train_loader))
        
        # Test batch debería tener 10 pacientes (del fold 0)
        assert test_batch['patient'].x.shape[0] == 8
        # Train batch debería tener datos de los otros folds
        assert train_batch['patient'].x.shape[0] == 4  # Máximo 20 (2 folds)

        # Verificar que el total de pacientes del test_loader es 10
        assert sum(b['patient'].x.shape[0] for b in test_loader) == 10
        assert sum(b['patient'].x.shape[0] for b in train_loader) == 20

    def test_create_fold_dataloaders_invalid_fold(self, sample_graph_data):
        """Test manejo de fold inválido."""
        graph_file, _ = sample_graph_data
        
        cv_validator = GNNCrossValidator(graph_data_path=graph_file)
        cv_validator.load_graph_folds()
        
        with pytest.raises(ValueError):
            cv_validator.create_fold_dataloaders(test_fold=5)  # Solo hay 3 folds
    
    def test_create_model_instance(self, temp_dirs):
        """Test creación de instancias de modelos."""
        # Crear archivo de configuración temporal
        temp_config = {
            'models': {
                'test_v1': {
                    'model_type': 'v1',
                    'use_vaf': False,
                    'hyperparameters': {
                        'hidden_gene_dim': 32,
                        'hidden_dim': 16,
                        'out_dim': 1,
                        'learning_rate': 0.001,
                        'hidden_vaf': 8
                    }
                }
            }
        }
        
        config_file = Path(temp_dirs) / "test_config.yaml"
        with open(config_file, 'w') as f:
            yaml.dump(temp_config, f)
        
        cv_validator = GNNCrossValidator(
            graph_data_path="dummy_path",
            config_path=str(config_file),
            results_dir=temp_dirs
        )
        
        # Simular dimensiones cargadas
        cv_validator.patient_feat_dim = 13
        cv_validator.gene_feat_dim = 16
        # Test creación de modelo
        model_v1 = cv_validator.create_model_instance('test_v1')
        assert hasattr(model_v1, 'training_step')
        print("Test 1")

        # Test con configuración manual
        manual_config = {
            'model_type': 'v1',
            'use_vaf': False,
            'hidden_gene_dim': 32,
            'hidden_dim': 16,
            'out_dim': 1,
            'learning_rate': 0.001,
            'hidden_vaf': 8
        }
        model_manual = cv_validator.create_model_instance('manual_model', manual_config)
        assert hasattr(model_manual, 'training_step')
        
        # Test tipo inválido
        with pytest.raises(ValueError):
            cv_validator.create_model_instance('non_existent_model')
    
    def test_define_checkpoint_callback(self, temp_dirs):
        """Test definición de callback de checkpoint."""
        cv_validator = GNNCrossValidator(
            graph_data_path="dummy_path",
            results_dir=temp_dirs
        )
        
        callback = cv_validator.define_checkpoint_callback("test_model", 0)
        
        assert callback.monitor == "val_loss"
        assert callback.mode == "min"
        assert callback.save_top_k == 1
        assert "fold_0" in callback.dirpath


class TestGNNCVPipeline:
    """Tests para la función de pipeline completo."""
    
    @pytest.fixture
    def temp_setup(self):
        """Fixture para setup temporal completo."""
        temp_dir = tempfile.mkdtemp()
        
        # Crear datos de grafos mínimos
        def create_minimal_graph(patient_idx: int) -> HeteroData:
            graph = HeteroData()
            graph['patient'].x = torch.randn(1, 13)
            graph['patient'].y = torch.tensor([[1.0, 1.0]], dtype=torch.float)
            graph['patient'].ID = patient_idx
            graph['gene'].x = torch.randn(2, 16)
            graph['gene'].weights = torch.rand(2)
            graph['gene', 'patient'].edge_index = torch.tensor([[0, 1], [0, 0]], dtype=torch.long)
            return graph
        
        # 2 folds con 5 pacientes cada uno (datos mínimos)
        graph_data = {}
        for fold in range(2):
            graphs = [create_minimal_graph(i + fold * 5) for i in range(5)]
            graph_data[f'fold_{fold}'] = graphs
        
        graph_file = Path(temp_dir) / "minimal_graphs.pt"
        torch.save(graph_data, graph_file)
        
        yield str(graph_file), temp_dir
        shutil.rmtree(temp_dir)
    
    def test_run_gnn_cv_pipeline_minimal(self, temp_setup):
        """Test del pipeline CV con configuración mínima."""
        graph_file, temp_dir = temp_setup
        
        # Crear configuración de prueba
        test_config = {
            'models': {
                'test_v1': {
                    'model_type': 'v1',
                    'use_vaf': False,
                    'hyperparameters': {
                        'hidden_gene_dim': 16,
                        'hidden_dim': 8,
                        'out_dim': 1,
                        'learning_rate': 0.01,
                        'hidden_vaf': 4
                    },
                    'training': {
                        'max_epochs': 2,  # Muy pocas épocas para test rápido
                        'batch_size_train': 4,
                        'batch_size_test': 8
                    }
                }
            }
        }
        
        config_file = Path(temp_dir) / "test_config.yaml"
        with open(config_file, 'w') as f:
            yaml.dump(test_config, f)
        
        # Ejecutar pipeline
        results = run_gnn_cv_pipeline(
            graph_data_path=graph_file,
            config_path=str(config_file),
            results_dir=Path(temp_dir) / "results",
            models_to_run=['test_v1']
        )

        assert 'test_v1' in results

        result = results['test_v1']
        assert result['model_name'] == 'test_v1'
        assert result['n_folds'] == 2
        assert len(result['fold_results']) == 2
        assert 'mean_val_c_index' in result
        assert 'data_dimensions' in result
        
        # Verificar que se guardaron archivos
        results_dir = Path(temp_dir) / "results"
        assert (results_dir / "test_v1_cv_results.pkl").exists()
        assert (results_dir / "test_v1_cv_summary.csv").exists()
        assert (results_dir / "comparative_cv_summary.csv").exists()


def run_tests():
    """Ejecuta todos los tests."""
    print("Ejecutando tests para gnn_cv_functions.py...")
    
    # Configurar paths
    current_dir = Path.cwd()
    if not (current_dir / "scripts" / "3.GeneticScores").exists():
        print("Error: Ejecutar desde el directorio raíz del proyecto")
        return False
    
    try:
        # Importar y ejecutar tests básicos
        test_class = TestGNNCrossValidator()
        
        # Test de inicialización
        with tempfile.TemporaryDirectory() as temp_dir:
            test_class.test_initialization(temp_dir)
            print("✓ Test de inicialización pasado")
        
        # Test de creación de modelos
        with tempfile.TemporaryDirectory() as temp_dir:
            test_class.test_create_model_instance(temp_dir)
            print("✓ Test de creación de modelos pasado")
        
        print("\nTodos los tests básicos pasaron!")
        print("\nNOTA: Las nuevas funcionalidades incluyen:")
        print("- Configuración YAML para hiperparámetros")
        print("- Detección automática de dimensiones de entrada")
        print("- Configuración flexible por modelo")
        print("\nPara ejecutar tests completos con datos reales:")
        print("python -m pytest scripts/3.GeneticScores/tests/test_gnn_cv_functions.py -v")
        
        return True
        
    except Exception as e:
        print(f"Error en tests: {str(e)}")
        return False


if __name__ == "__main__":
    success = run_tests()
    sys.exit(0 if success else 1)
