#!/usr/bin/env python3
"""
Script de prueba simple para el pipeline unificado.
Prueba básica de funcionamiento de los módulos principales.
"""
import sys
import os
from pathlib import Path
import logging

# Configurar logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Añadir el directorio utils al path
script_dir = Path(__file__).parent
utils_dir = script_dir.parent / "utils"
sys.path.insert(0, str(utils_dir))

logger.info(f"Buscando módulos en: {utils_dir}")

try:
    from data_loader import DataLoader
    from clinical_processor import ClinicalProcessor
    from data_generators import DataGenerators
    logger.info("✓ Módulos importados correctamente")
except ImportError as e:
    logger.error(f"✗ Error importando módulos: {e}")
    sys.exit(1)

def test_basic_functionality():
    """Test básico de funcionalidad."""
    logger.info("Iniciando test básico...")
    
    # Configuración de prueba mínima simplificada
    test_config = {
        'data_paths': {
            'base_path': '/tmp',
            'clinical_mutations_file': 'test.tsv',
            'mutations_vaf_file': None,
            'gene_embeddings_file': None
        },
        'variable_processing': {
            'continuous_variables': ['Age', 'HB'],
            'categorical_variables': ['Sex'],
            'binary_variables': ['Complex_Karyotype'],
            'survival_variables': {
                'time_variable': 'OS',
                'status_variable': 'OS_STATUS'
            },
            'gene_selection': {
                'selected_genes': ['TP53', 'ASXL1'],
                'auto_selection_criteria': {
                    'min_mutation_frequency': 0.05,
                    'max_genes': 10
                }
            }
        },
        'data_generation': {
            'output_formats': ['table'],
            'table_format': {
                'include_all_clinical': True,
                'output_filename': 'test_table'
            },
            'graph_format': {
                'graph_type': 'patient_gene',
                'patient_node_features': 'clinical_only',
                'gene_node_features': {
                    'embedding_source': 'scgpt',
                    'embedding_dimension': 512
                },
                'edge_criteria': {
                    'use_vaf_filter': False,
                    'min_vaf_threshold': 0.05
                },
                'output_filename': 'test_graph'
            }
        },
        'cross_validation': {
            'n_folds': 5,
            'fold_column': 'fold'
        },
        'output': {
            'file_formats': {
                'tables': 'pickle',
                'graphs': 'pickle',
                'metadata': 'json'
            },
            'compression': None,
            'validate_output': True,
            'save_config_copy': True
        }
    }
    
    try:
        # Probar inicialización de módulos
        logger.info("Probando inicialización de DataLoader...")
        data_loader = DataLoader(test_config)
        logger.info("✓ DataLoader inicializado")
        
        logger.info("Probando inicialización de ClinicalProcessor...")
        clinical_processor = ClinicalProcessor(test_config)
        logger.info("✓ ClinicalProcessor inicializado")
        
        logger.info("Probando inicialización de DataGenerators...")
        data_generators = DataGenerators(test_config)
        logger.info("✓ DataGenerators inicializado")
        
        # Probar que los módulos tienen los métodos esperados
        logger.info("Verificando métodos disponibles...")
        
        # DataLoader
        assert hasattr(data_loader, 'load_clinical_data'), "DataLoader missing load_clinical_data"
        assert hasattr(data_loader, 'load_mutation_vaf_data'), "DataLoader missing load_mutation_vaf_data"
        assert hasattr(data_loader, 'load_gene_embeddings'), "DataLoader missing load_gene_embeddings"
        logger.info("✓ DataLoader tiene todos los métodos requeridos")
        
        # ClinicalProcessor
        assert hasattr(clinical_processor, 'prepare_clinical_features'), "ClinicalProcessor missing prepare_clinical_features"
        assert hasattr(clinical_processor, 'extract_survival_data'), "ClinicalProcessor missing extract_survival_data"
        assert hasattr(clinical_processor, 'extract_cv_folds'), "ClinicalProcessor missing extract_cv_folds"
        assert hasattr(clinical_processor, 'get_feature_matrix'), "ClinicalProcessor missing get_feature_matrix"
        logger.info("✓ ClinicalProcessor tiene todos los métodos requeridos")
        
        # DataGenerators
        assert hasattr(data_generators, 'generate_table_data'), "DataGenerators missing generate_table_data"
        assert hasattr(data_generators, 'generate_graph_data'), "DataGenerators missing generate_graph_data"
        logger.info("✓ DataGenerators tiene todos los métodos requeridos")
        
        logger.info("🎉 Test básico completado exitosamente")
        logger.info("Todos los módulos están correctamente inicializados y tienen los métodos esperados")
        return True
        
    except Exception as e:
        logger.error(f"✗ Error en test básico: {e}")
        import traceback
        logger.error(f"Traceback: {traceback.format_exc()}")
        return False

if __name__ == "__main__":
    logger.info("="*50)
    logger.info("TEST BÁSICO DEL PIPELINE UNIFICADO")
    logger.info("="*50)
    
    success = test_basic_functionality()
    
    if success:
        logger.info("\n✅ Test básico EXITOSO - El pipeline está listo para usar")
    else:
        logger.error("\n❌ Test básico FALLÓ - Revisar errores arriba")
    
    sys.exit(0 if success else 1)
