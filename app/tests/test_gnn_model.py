"""
Tests mínimos y completos para el módulo gnn_model.py

Tests que cubren inicialización, validación, predicción y manejo de errores.
Requiere archivos reales de transformer y gene embeddings para funcionar.
"""

import pytest
import tempfile
import os
import yaml
from unittest.mock import Mock, patch

from app.models.gnn_model import MDSGNNModel, load_mds_gnn_model
from app.models.base_model import PatientData, RiskPrediction, ValidationResult


# Fixture de configuración con paths reales
@pytest.fixture
def real_config():
    """Configuración con paths reales a archivos existentes."""
    return {
        'metadata': {'model_version': '1.0.0'},
        'models': {
            'gnn': {
                'checkpoint_path': 'results/gnn/full_input_model/model_checkpoints/fold_0/v1_boolean_fold_0-v1.ckpt',
                'architecture': {
                    'patient_feat_dim': 17,
                    'gene_feat_dim': 16,
                    'hidden_gene_dim': 32,
                    'hidden_dim': 16
                },
                'hyperparameters': {
                    'learning_rate': 0.001
                }
            }
        },
        'data': {
            'clinical_transformer': {
                'file_path': 'results/gnn/preprocess/full_input/fitted_transformers/preprocessing_pipeline.pkl'
            },
            'gene_embeddings': {
                'file_path': 'results/preprocess/gene_embedding_reduced.tsv',
                'supported_genes': ['TP53', 'SF3B1', 'ASXL1', 'SRSF2', 'TET2', 'DNMT3A', 'IDH1'],
                'embedding_dim': 16
            }
        },
        'validation': {
            'clinical_ranges': {
                'AGE': {'min': 18, 'max': 100, 'unit': 'años'},
                'BM_BLAST': {'min': 0, 'max': 100, 'unit': '%'},
                'WBC': {'min': 0.1, 'max': 500, 'unit': '10^9/L'},
                'ANC': {'min': 0, 'max': 100, 'unit': '10^9/L'},
                'HB': {'min': 3, 'max': 20, 'unit': 'g/dL'},
                'PLT': {'min': 1, 'max': 2000, 'unit': '10^9/L'},
                'MONOCYTES': {'min': 0, 'max': 50, 'unit': '10^9/L'}
            },
            'categorical_values': {
                'SEX': ['M', 'F'],
                'CYTO_IPSSR': ['Very Good', 'Good', 'Intermediate', 'Poor', 'Very Poor'],
                'complex': ['complex', 'non-complex']
            },
            'warnings': {
                'missing_mutations': {
                    'message': 'No se detectaron mutaciones'
                }
            }
        },
        'risk_categories': {
            'thresholds': {
                'very_low': -1.5,
                'low': -0.5,
                'moderate_low': 0.0,
                'moderate_high': 0.5,
                'high': 1.0
            },
            'labels': {
                'very_low': 'Muy Bajo',
                'low': 'Bajo',
                'moderate_low': 'Moderado-Bajo',
                'moderate_high': 'Moderado-Alto',
                'high': 'Alto',
                'very_high': 'Muy Alto'
            },
            'colors': {
                'very_low': '#28a745',
                'low': '#6c757d',
                'moderate_low': '#ffc107',
                'moderate_high': '#fd7e14',
                'high': '#dc3545',
                'very_high': '#6f42c1'
            }
        }
    }


@pytest.fixture
def config_file(real_config):
    """Crea archivo de configuración temporal."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
        yaml.dump(real_config, f)
        config_path = f.name
    
    yield config_path
    
    # Cleanup
    os.unlink(config_path)


@pytest.fixture
def sample_patient_data():
    """Datos de paciente válidos para tests."""
    return PatientData(
        AGE=68,
        SEX='M',
        BM_BLAST=13,
        WBC=1,
        ANC=0,
        HB=12,
        PLT=110,
        MONOCYTES=0.5,
        CYTO_IPSSR='Good',
        complex="non-complex",
        plus8=False,
        del7=False,
        del20q=False,
        del7q=False,
        delY=False,
        del5q=False,
        mutations=['IDH1', 'ASXL1', 'SRSF2']
    )


class TestMDSGNNModel:
    """Tests principales para la clase MDSGNNModel."""
    
    def test_init_with_real_files(self, config_file):
        """Test inicialización con archivos reales."""
        # Solo mockear el checkpoint del modelo GNN si no existe
        with patch('app.models.gnn_model.Path.exists') as mock_exists:
            # def side_effect(path):
            #     # Mockear solo el checkpoint del GNN, permitir otros archivos
            #     if 'v1_boolean_fold_0-v1.ckpt' in str(path):
            #         return False
            #     return True  # Permite que otros archivos existan
            
            # mock_exists.side_effect = side_effect
            
            model = MDSGNNModel(config_file)
            
            # Verificar que todos los componentes se inicializaron
            assert model.config is not None
            assert model.model_config is not None
            assert model.data_config is not None
            assert model.validation_config is not None
            assert model.risk_config is not None
            assert model.model is not None
            assert model.transformer is not None
            assert model.gene_embeddings is not None
            assert model.graph_builder is not None
    
    def test_validate_input_valid_data(self, config_file, sample_patient_data):
        """Test validación con datos válidos."""
        with patch('app.models.gnn_model.Path.exists') as mock_exists:
            # def side_effect(path):
            #     if 'v1_boolean_fold_0-v1.ckpt' in str(path):
            #         return False
            #     return True
            # mock_exists.side_effect = side_effect
            
            model = MDSGNNModel(config_file)
            validation = model.validate_input(sample_patient_data)
            
            #assert isinstance(validation, ValidationResult)
            assert validation.is_valid is True
            assert len(validation.errors) == 0
            assert len(validation.validated_fields) > 0
    
    def test_validate_input_invalid_data(self, config_file):
        """Test validación con datos inválidos."""
        with patch('app.models.gnn_model.Path.exists') as mock_exists:
            # def side_effect(path):
            #     if 'v1_boolean_fold_0-v1.ckpt' in str(path):
            #         return False
            #     return True
            # mock_exists.side_effect = side_effect
            
            model = MDSGNNModel(config_file)
            
            invalid_data = PatientData(
                AGE=150,  # Fuera de rango
                SEX='X',  # Sexo inválido
                BM_BLAST=8.5,
                WBC=3.2,
                ANC=1.8,
                HB=1.0,  # Fuera de rango
                PLT=120,
                MONOCYTES=0.5,
                CYTO_IPSSR='Invalid',  # Categoría inválida
                complex="complex",
                plus8=False,
                del7=True,
                del20q=False,
                del7q=False,
                delY=False,
                del5q=False,
                mutations=['TP53']
            )
            
            validation = model.validate_input(invalid_data)
            
            assert validation.is_valid is False
            assert len(validation.errors) > 0
            assert any('AGE' in error for error in validation.errors)
            assert any('HB' in error for error in validation.errors)
            assert any('Sexo inválido' in error for error in validation.errors)
            assert any('citogenético inválido' in error for error in validation.errors)
    
    def test_predict_valid_data(self, config_file, sample_patient_data):
        """Test predicción con datos válidos."""
        
        model = MDSGNNModel(config_file)
        prediction = model.predict(sample_patient_data)
            
        #assert isinstance(prediction, RiskPrediction)
        assert isinstance(prediction.raw_score, float)
        assert prediction.raw_score == 0.6024637391136364
        assert isinstance(prediction.risk_category, str)
        assert prediction.risk_category in ['very_low', 'low', 'moderate_low', 
                                              'moderate_high', 'high', 'very_high']
            #assert isinstance(prediction.risk_probability, float)
            #assert 0.0 <= prediction.risk_probability <= 1.0
        assert isinstance(prediction.hazard_ratio, float)
        assert prediction.hazard_ratio > 0.0
            #assert isinstance(prediction.confidence, float)
            #assert 0.5 <= prediction.confidence <= 1.0
        assert prediction.processing_time > 0
        assert prediction.model_version == '1.0.0'
    
    def test_predict_invalid_data(self, config_file):
        """Test predicción con datos inválidos (debe fallar)."""
        with patch('app.models.gnn_model.Path.exists') as mock_exists:
            # def side_effect(path):
            #     if 'v1_boolean_fold_0-v1.ckpt' in str(path):
            #         return False
            #     return True
            # mock_exists.side_effect = side_effect
            
            model = MDSGNNModel(config_file)
            
            invalid_data = PatientData(
                AGE=150,  # Edad inválida
                SEX='X',  # Sexo inválido
                BM_BLAST=8.5,
                WBC=3.2,
                ANC=1.8,
                HB=9.5,
                PLT=120,
                MONOCYTES=0.5,
                CYTO_IPSSR='Invalid',  # Categoría inválida
                complex="complex",
                plus8=False,
                del7=True,
                del20q=False,
                del7q=False,
                delY=False,
                del5q=False,
                mutations=['TP53']
            )
            
            with pytest.raises(RuntimeError):
                model.predict(invalid_data)
    
    def test_get_risk_category(self, config_file):
        """Test categorización de riesgo."""
        with patch('app.models.gnn_model.Path.exists') as mock_exists:
            # def side_effect(path):
            #     if 'v1_boolean_fold_0-v1.ckpt' in str(path):
            #         return False
            #     return True
            # mock_exists.side_effect = side_effect
            
            model = MDSGNNModel(config_file)
            
            # Test diferentes scores
            assert model.get_risk_category(-2.0) == 'very_low'
            assert model.get_risk_category(-1.0) == 'low'
            assert model.get_risk_category(-0.25) == 'moderate_low'
            assert model.get_risk_category(0.25) == 'moderate_high'
            assert model.get_risk_category(0.75) == 'high'
            assert model.get_risk_category(1.5) == 'very_high'
    
    def test_get_risk_category_labels_and_colors(self, config_file):
        """Test obtención de etiquetas y colores de riesgo."""
        with patch('app.models.gnn_model.Path.exists') as mock_exists:
            # def side_effect(path):
            #     if 'v1_boolean_fold_0-v1.ckpt' in str(path):
            #         return False
            #     return True
            # mock_exists.side_effect = side_effect
            
            model = MDSGNNModel(config_file)
            
            # Test etiquetas
            assert model.get_risk_category_label('very_low') == 'Muy Bajo'
            assert model.get_risk_category_label('high') == 'Alto'
            assert model.get_risk_category_label('unknown') == 'unknown'
            
            # Test colores
            assert model.get_risk_category_color('very_low') == '#28a745'
            assert model.get_risk_category_color('high') == '#dc3545'
            assert model.get_risk_category_color('unknown') == '#666666'
    
    def test_get_model_info(self, config_file):
        """Test información del modelo."""
        with patch('app.models.gnn_model.Path.exists') as mock_exists:
            # def side_effect(path):
            #     if 'v1_boolean_fold_0-v1.ckpt' in str(path):
            #         return False
            #     return True
            # mock_exists.side_effect = side_effect
            
            model = MDSGNNModel(config_file)
            info = model.get_model_info()
            
            assert isinstance(info, dict)
            assert 'name' in info
            assert 'type' in info
            assert 'version' in info
            assert 'framework' in info
            assert 'description' in info
            assert info['version'] == '1.0.0'


class TestLoadMDSGNNModel:
    """Tests para la función de conveniencia."""
    
    def test_load_with_default_config(self):
        """Test carga con configuración por defecto."""
        with patch('app.models.gnn_model.MDSGNNModel') as mock_model_class:
            mock_instance = Mock()
            mock_model_class.return_value = mock_instance
            
            result = load_mds_gnn_model()
            
            mock_model_class.assert_called_once_with("app/config.yaml")
            assert result == mock_instance
    
    def test_load_with_custom_config(self):
        """Test carga con configuración personalizada."""
        custom_path = "/custom/config.yaml"
        
        with patch('app.models.gnn_model.MDSGNNModel') as mock_model_class:
            mock_instance = Mock()
            mock_model_class.return_value = mock_instance
            
            result = load_mds_gnn_model(custom_path)
            
            mock_model_class.assert_called_once_with(custom_path)
            assert result == mock_instance


class TestErrorHandling:
    """Tests para manejo de errores."""
    
    def test_invalid_config_path(self):
        """Test con path de configuración inexistente."""
        with pytest.raises(FileNotFoundError):
            MDSGNNModel('/path/to/nonexistent/config.yaml')
    
    def test_invalid_yaml_format(self):
        """Test con YAML inválido."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            f.write("invalid: yaml: content: [")
            config_path = f.name
        
        try:
            with pytest.raises(yaml.YAMLError):
                MDSGNNModel(config_path)
        finally:
            os.unlink(config_path)
    
    def test_missing_required_files(self, real_config):
        """Test cuando faltan archivos requeridos."""
        # Crear config con paths inexistentes
        real_config['data']['clinical_transformer']['file_path'] = '/nonexistent/transformer.pkl'
        real_config['data']['gene_embeddings']['file_path'] = '/nonexistent/embeddings.tsv'
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            yaml.dump(real_config, f)
            config_path = f.name
        
        try:
            with pytest.raises(RuntimeError):
                MDSGNNModel(config_path)
        finally:
            os.unlink(config_path)
