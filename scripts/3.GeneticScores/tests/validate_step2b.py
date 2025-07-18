#!/usr/bin/env python3
"""
Script de validación completa para el Paso 2B del pipeline unificado.
Valida la funcionalidad de todos los módulos implementados.
"""
import sys
import os
from pathlib import Path
import logging
import pandas as pd
import numpy as np
import tempfile
import yaml
import json

# Configurar logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Añadir el directorio utils al path para importar módulos
script_dir = Path(__file__).parent
utils_dir = script_dir.parent / "utils"
sys.path.insert(0, str(utils_dir))

logger.info(f"Script ejecutándose desde: {script_dir}")
logger.info(f"Buscando módulos en: {utils_dir}")

# Importar módulos del pipeline
try:
    from data_loader import DataLoader
    from clinical_processor import ClinicalProcessor
    from data_generators import DataGenerators
    logger.info("✓ Módulos del pipeline importados correctamente")
except ImportError as e:
    logger.error(f"✗ Error importando módulos: {e}")
    logger.error(f"  Verificar que estés en el directorio correcto")
    logger.error(f"  Path actual: {os.getcwd()}")
    logger.error(f"  Buscando en: {utils_dir}")
    sys.exit(1)

# Verificar dependencias adicionales
try:
    import torch
    from torch_geometric.data import HeteroData
    logger.info("✓ Dependencias PyTorch Geometric disponibles")
except ImportError as e:
    logger.warning(f"⚠ PyTorch Geometric no disponible: {e}")
    logger.warning("  Los tests de grafos podrían fallar")


class Step2BValidator:
    """Validador completo para el Paso 2B del pipeline."""
    
    def __init__(self):
        """Inicializa el validador."""
        self.temp_dir = None
        self.test_config = None
        self.test_data = None
        
        logger.info("Inicializando validador Step 2B...")
    
    def setup_test_environment(self):
        """Configura el entorno de prueba con datos simulados."""
        logger.info("Configurando entorno de prueba...")
        
        # Crear directorio temporal
        self.temp_dir = Path(tempfile.mkdtemp())
        logger.info(f"Directorio temporal creado: {self.temp_dir}")
        
        # Crear datos de prueba
        self._create_test_data()
        
        # Crear configuración de prueba
        self._create_test_config()
        
        logger.info("✓ Entorno de prueba configurado")
    
    def _create_test_data(self):
        """Crea datos de prueba simulados."""
        logger.info("Creando datos de prueba...")
        
        # Datos clínicos con mutaciones combinadas
        np.random.seed(42)
        n_patients = 100
        
        clinical_data = {
            'ID': [f'Patient_{i:03d}' for i in range(n_patients)],
            'Age': np.random.normal(65, 10, n_patients),
            'Sex': np.random.choice(['M', 'F'], n_patients),
            'HB': np.random.normal(9.5, 2, n_patients),
            'ANC': np.random.lognormal(0, 1, n_patients),
            'PLT': np.random.normal(150, 50, n_patients),
            'BM_BLAST': np.random.uniform(0, 20, n_patients),
            'PB_BLAST': np.random.uniform(0, 10, n_patients),
            'Cytogenetics_IPSS': np.random.choice(['Good', 'Intermediate', 'Poor'], n_patients),
            'WHO_Classification': np.random.choice(['RCUD', 'RARS', 'RCMD', 'RAEB-1', 'RAEB-2'], n_patients),
            'Complex_Karyotype': np.random.choice([0.0, 1.0], n_patients, p=[0.8, 0.2]),
            'OS': np.random.exponential(24, n_patients),  # Meses
            'OS_STATUS': np.random.choice([0.0, 1.0], n_patients, p=[0.3, 0.7]),
            'fold': np.random.randint(1, 11, n_patients),
            # Mutaciones binarizadas (genes frecuentes en MDS)
            'TP53': np.random.choice([0.0, 1.0], n_patients, p=[0.85, 0.15]),
            'ASXL1': np.random.choice([0.0, 1.0], n_patients, p=[0.75, 0.25]),
            'SF3B1': np.random.choice([0.0, 1.0], n_patients, p=[0.8, 0.2]),
            'DNMT3A': np.random.choice([0.0, 1.0], n_patients, p=[0.9, 0.1]),
            'SRSF2': np.random.choice([0.0, 1.0], n_patients, p=[0.85, 0.15]),
            'TET2': np.random.choice([0.0, 1.0], n_patients, p=[0.88, 0.12]),
            'RUNX1': np.random.choice([0.0, 1.0], n_patients, p=[0.92, 0.08]),
            'U2AF1': np.random.choice([0.0, 1.0], n_patients, p=[0.95, 0.05])
        }
        
        df_clinical = pd.DataFrame(clinical_data)
        clinical_file = self.temp_dir / "clinical_mutations_combined.tsv"
        df_clinical.to_csv(clinical_file, sep='\t', index=False)
        
        # Datos de mutaciones con VAF
        mutation_data = []
        genes = ['TP53', 'ASXL1', 'SF3B1', 'DNMT3A', 'SRSF2', 'TET2', 'RUNX1', 'U2AF1']
        
        for _, row in df_clinical.iterrows():
            patient_id = row['ID']
            for gene in genes:
                if row[gene] == 1:  # Si tiene la mutación
                    vaf = np.random.uniform(0.05, 0.8)*100  # VAF entre 5% y 80%
                    mutation_data.append({
                        'ID': patient_id,
                        'GENE': gene,
                        'VAF': vaf
                    })
        
        df_mutations = pd.DataFrame(mutation_data)
        mutations_file = self.temp_dir / "mutations_prioritized.tsv"
        df_mutations.to_csv(mutations_file, sep='\t', index=False)
        
        # Embeddings de genes simulados
        # Crear embeddings de genes simulados y guardarlos en un DataFrame CSV
        embedding_dim = 4
        gene_embeddings = []
        for gene in genes:
            embedding = np.random.normal(0, 1, embedding_dim)
            gene_embeddings.append([gene] + embedding.tolist())
        
        df_embeddings = pd.DataFrame(
            gene_embeddings,
            columns=["gene"] + [f"emb_{i}" for i in range(embedding_dim)]
        )
        embeddings_file = self.temp_dir / "gene_embeddings.csv"
        df_embeddings.to_csv(embeddings_file, index=False, sep = "\t")
        
        logger.info(f"✓ Datos de prueba creados: {n_patients} pacientes, {len(df_mutations)} mutaciones")
    
    def _create_test_config(self):
        """Crea configuración de prueba."""
        self.test_config = {
            'data_paths': {
                'base_path': str(self.temp_dir),
                'clinical_mutations_file': 'clinical_mutations_combined.tsv',
                'mutations_vaf_file': 'mutations_prioritized.tsv',
                'gene_embeddings_file': 'gene_embeddings.csv'
            },
            'variable_processing': {
                'continuous_variables': ['Age', 'HB', 'ANC', 'PLT', 'BM_BLAST', 'PB_BLAST'],
                'categorical_variables': ['Sex', 'Cytogenetics_IPSS', 'WHO_Classification'],
                'binary_variables': ['Complex_Karyotype'],
                'survival_variables': {
                    'time_variable': 'OS',
                    'status_variable': 'OS_STATUS'
                },
                'gene_selection': {
                    'selected_genes': ['TP53', 'ASXL1', 'SF3B1', 'DNMT3A', 'SRSF2'],
                    'auto_selection_criteria': {
                        'min_mutation_frequency': 0.05,
                        'max_genes': 10
                    }
                }
            },
            'data_generation': {
                'output_formats': ['table', 'graph'],
                'table_format': {
                    'include_all_clinical': True,
                    'output_filename': 'test_table_data'
                },
                'graph_format': {
                    'graph_type': 'patient_gene',
                    'patient_node_features': 'clinical_only',
                    'gene_node_features': {
                        'embedding_source': 'scgpt',
                        'embedding_dimension': 4
                    },
                    'edge_criteria': {
                        'use_vaf_filter': True,
                        'min_vaf_threshold': 0.05
                    },
                    'output_filename': 'test_graph_data'
                }
            },
            'cross_validation': {
                'n_folds': 10,
                'fold_column': 'fold'
            },
            'output': {
                'file_formats': {
                    'tables': 'csv',
                    'graphs': 'pickle',
                    'metadata': 'json'
                },
                'compression': None,
                'validate_output': True,
                'save_config_copy': True
            }
        }
    
    def test_data_loader(self):
        """Prueba el módulo DataLoader."""
        logger.info("🧪 Probando DataLoader...")
        
        try:
            data_loader = DataLoader(self.test_config)
            
            # Probar carga de datos clínicos
            df_clinical = data_loader.load_clinical_data()
            assert df_clinical is not None, "Error: datos clínicos no cargados"
            assert df_clinical.shape[0] == 100, f"Error: esperados 100 pacientes, obtenidos {df_clinical.shape[0]}"
            logger.info(f"✓ Datos clínicos cargados: {df_clinical.shape}")
            
            # Probar carga de mutaciones VAF
            df_mutations = data_loader.load_mutation_vaf_data()
            assert df_mutations is not None, "Error: datos de mutaciones no cargados"
            assert 'VAF' in df_mutations.columns, "Error: columna VAF no encontrada"
            logger.info(f"✓ Datos de mutaciones VAF cargados: {df_mutations.shape}")
            
            # Probar carga de embeddings
            embeddings = data_loader.load_gene_embeddings()
            assert embeddings is not None, "Error: embeddings no cargados"
            assert len(embeddings) > 0, "Error: embeddings vacíos"
            logger.info(f"✓ Embeddings de genes cargados: {len(embeddings)} genes")
            
            logger.info("✅ DataLoader: PASSED")
            return True
            
        except Exception as e:
            logger.error(f"❌ DataLoader: FAILED - {e}")
            return False
    
    def test_clinical_processor(self):
        """Prueba el módulo ClinicalProcessor."""
        logger.info("🧪 Probando ClinicalProcessor...")
        
        try:
            # Cargar datos
            data_loader = DataLoader(self.test_config)
            df_clinical = data_loader.load_clinical_data()
            
            # Procesar datos
            processor = ClinicalProcessor(self.test_config)
            selected_genes = self.test_config['variable_processing']['gene_selection']['selected_genes']
            
            df_processed, metadata = processor.prepare_clinical_features(df_clinical, selected_genes)
            assert df_processed is not None, "Error: datos procesados nulos"
            assert 'feature_columns' in metadata, "Error: metadatos incompletos"
            logger.info(f"✓ Features procesadas: {len(metadata['feature_columns'])} columnas")
            
            # Probar extracción de datos de supervivencia
            times, events = processor.extract_survival_data(df_clinical)
            assert len(times) == 100, "Error: longitud incorrecta de tiempos"
            assert len(events) == 100, "Error: longitud incorrecta de eventos"
            logger.info(f"✓ Datos de supervivencia extraídos: {len(times)} observaciones")
            
            # Probar extracción de folds CV
            folds = processor.extract_cv_folds(df_clinical)
            assert len(folds) == 100, "Error: longitud incorrecta de folds"
            assert len(np.unique(folds)) == 10, "Error: número incorrecto de folds únicos"
            logger.info(f"✓ Folds CV extraídos: {len(np.unique(folds))} folds únicos")
            
            # Probar extracción de matriz de features
            X = processor.get_feature_matrix(df_processed, metadata['feature_columns'])
            assert X.shape[0] == 100, "Error: número incorrecto de muestras"
            assert X.shape[1] == len(metadata['feature_columns']), "Error: número incorrecto de features"
            logger.info(f"✓ Matriz de features extraída: {X.shape}")
            
            logger.info("✅ ClinicalProcessor: PASSED")
            return True
            
        except Exception as e:
            logger.error(f"❌ ClinicalProcessor: FAILED - {e}")
            return False
    
    def test_data_generators(self):
        """Prueba el módulo DataGenerators."""
        logger.info("🧪 Probando DataGenerators...")
        
        try:
            # Preparar datos
            data_loader = DataLoader(self.test_config)
            df_clinical = data_loader.load_clinical_data()
            df_mutations = data_loader.load_mutation_vaf_data()
            embeddings = data_loader.load_gene_embeddings()
            
            processor = ClinicalProcessor(self.test_config)
            selected_genes = self.test_config['variable_processing']['gene_selection']['selected_genes']
            df_processed, metadata = processor.prepare_clinical_features(df_clinical, selected_genes)
            
            # Generar datos
            generators = DataGenerators(self.test_config)
            
            # Crear directorio de salida
            output_dir = self.temp_dir / "output"
            output_dir.mkdir(exist_ok=True)
            
            # Probar generación de tabla
            logger.info("  Generando datos de tabla...")
            table_data = generators.generate_table_data(
                df_processed, metadata, output_dir
            )
            assert table_data is not None, "Error: datos de tabla no generados"
            assert 'file_path' in table_data, "Error: ruta de archivo no encontrada"
            logger.info("✓ Datos de tabla generados")
            
            # Probar generación de grafo
            logger.info("  Generando datos de grafo...")
            try:
                graph_data = generators.generate_graph_data(
                    df_processed, metadata, df_mutations, embeddings, output_dir
                )
                assert graph_data is not None, "Error: datos de grafo no generados"
                assert 'file_path' in graph_data, "Error: ruta de archivo no encontrada"
                logger.info("✓ Datos de grafo generados")
            except ImportError as e:
                logger.warning(f"⚠ Test de grafo omitido (dependencias faltantes): {e}")
            
            logger.info("✅ DataGenerators: PASSED")
            return True
            
        except Exception as e:
            logger.error(f"❌ DataGenerators: FAILED - {e}")
            import traceback
            logger.error(f"Traceback: {traceback.format_exc()}")
            return False
    
    def test_pipeline_integration(self):
        """Prueba la integración completa del pipeline."""
        logger.info("🧪 Probando integración completa del pipeline...")
        
        try:
            # Simular ejecución del pipeline completo
            data_loader = DataLoader(self.test_config)
            processor = ClinicalProcessor(self.test_config)
            generators = DataGenerators(self.test_config)
            
            # Cargar datos
            logger.info("  Cargando datos...")
            df_clinical = data_loader.load_clinical_data()
            df_mutations = data_loader.load_mutation_vaf_data()
            embeddings = data_loader.load_gene_embeddings()
            
            # Procesar datos
            logger.info("  Procesando datos...")
            selected_genes = self.test_config['variable_processing']['gene_selection']['selected_genes']
            df_processed, metadata = processor.prepare_clinical_features(df_clinical, selected_genes)
            
            # Generar datos de salida
            output_dir = self.temp_dir / "integration_test"
            output_dir.mkdir(exist_ok=True)
            
            # Generar ambos formatos
            logger.info("  Generando formatos de salida...")
            for output_format in self.test_config['data_generation']['output_formats']:
                logger.info(f"    - Formato: {output_format}")
                
                if output_format == 'table':
                    table_data = generators.generate_table_data(
                        df_processed, metadata, output_dir
                    )
                    assert table_data is not None, "Error en generación de tabla"
                    logger.info("    ✓ Tabla generada exitosamente")
                    
                elif output_format == 'graph':
                    try:
                        graph_data = generators.generate_graph_data(
                            df_processed, metadata, df_mutations, embeddings, output_dir
                        )
                        assert graph_data is not None, "Error en generación de grafo"
                        logger.info("    ✓ Grafo generado exitosamente")
                    except ImportError as e:
                        logger.warning(f"    ⚠ Grafo omitido (dependencias faltantes): {e}")
            
            logger.info("✅ Integración completa: PASSED")
            return True
            
        except Exception as e:
            logger.error(f"❌ Integración completa: FAILED - {e}")
            import traceback
            logger.error(f"Traceback: {traceback.format_exc()}")
            return False
    
    def test_error_handling(self):
        """Prueba el manejo de errores."""
        logger.info("🧪 Probando manejo de errores...")
        
        try:
            # Configuración con archivos inexistentes
            bad_config = self.test_config.copy()
            bad_config['data_paths']['clinical_mutations_file'] = 'nonexistent.tsv'
            
            data_loader = DataLoader(bad_config)
            
            # Esto debería lanzar una excepción
            try:
                df_clinical = data_loader.load_clinical_data()
                logger.error("❌ Error: debería haber fallado con archivo inexistente")
                return False
            except FileNotFoundError:
                logger.info("✓ Error correctamente manejado: archivo inexistente")
            
            # Configuración con columnas faltantes
            df_test = pd.DataFrame({'col1': [1, 2, 3]})
            processor = ClinicalProcessor(self.test_config)
            
            try:
                times, events = processor.extract_survival_data(df_test)
                logger.error("❌ Error: debería haber fallado con columnas faltantes")
                return False
            except ValueError:
                logger.info("✓ Error correctamente manejado: columnas faltantes")
            
            logger.info("✅ Manejo de errores: PASSED")
            return True
            
        except Exception as e:
            logger.error(f"❌ Manejo de errores: FAILED - {e}")
            return False
    
    def cleanup(self):
        """Limpia el entorno de prueba."""
        if self.temp_dir and self.temp_dir.exists():
            import shutil
            shutil.rmtree(self.temp_dir)
            logger.info(f"✓ Directorio temporal eliminado: {self.temp_dir}")
    
    def run_all_tests(self):
        """Ejecuta todas las pruebas."""
        logger.info("🚀 Iniciando validación completa del Paso 2B...")
        
        tests = [
            ("DataLoader", self.test_data_loader),
            ("ClinicalProcessor", self.test_clinical_processor),
            ("DataGenerators", self.test_data_generators),
            ("Integración Completa", self.test_pipeline_integration),
            ("Manejo de Errores", self.test_error_handling)
        ]
        
        results = []
        
        try:
            self.setup_test_environment()
            
            for test_name, test_func in tests:
                logger.info(f"\n{'='*50}")
                logger.info(f"Ejecutando: {test_name}")
                logger.info(f"{'='*50}")
                
                success = test_func()
                results.append((test_name, success))
        
        finally:
            print("\nLimpiando entorno de prueba...")
            #self.cleanup()
        
        # Resumen de resultados
        logger.info(f"\n{'='*60}")
        logger.info("📊 RESUMEN DE VALIDACIÓN PASO 2B")
        logger.info(f"{'='*60}")
        
        passed = 0
        failed = 0
        
        for test_name, success in results:
            status = "✅ PASSED" if success else "❌ FAILED"
            logger.info(f"{test_name:30} : {status}")
            
            if success:
                passed += 1
            else:
                failed += 1
        
        logger.info(f"\n{'='*60}")
        logger.info(f"Total: {len(results)} | Passed: {passed} | Failed: {failed}")
        
        if failed == 0:
            logger.info("🎉 VALIDACIÓN PASO 2B: COMPLETAMENTE EXITOSA")
            return True
        else:
            logger.error(f"💥 VALIDACIÓN PASO 2B: {failed} PRUEBAS FALLARON")
            return False


def main():
    """Función principal."""
    logger.info("="*60)
    logger.info("VALIDACIÓN DEL PASO 2B - PIPELINE UNIFICADO")
    logger.info("="*60)
    logger.info(f"Directorio de trabajo: {os.getcwd()}")
    logger.info(f"Python version: {sys.version}")
    
    validator = Step2BValidator()
    success = validator.run_all_tests()
    
    if success:
        logger.info("\n🎉 ¡Todos los tests pasaron! El pipeline está listo para usar.")
    else:
        logger.error("\n💥 Algunos tests fallaron. Revisar los errores arriba.")
    
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
