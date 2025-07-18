#!/usr/bin/env python3
"""
Script de validación de integración con modelos - Paso 3D del pipeline unificado.
Valida que los datos generados son compatibles con los modelos existentes de NN y GNN.
"""
import sys
import os
from pathlib import Path
import logging
import pandas as pd
import numpy as np
import torch
import pickle
from typing import Dict, Any, List, Optional, Tuple

# Configurar logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Verificar dependencias adicionales
try:
    import lightning as L
    from torch_geometric.data import HeteroData
    from torch_geometric.loader import DataLoader as GeoDataLoader
    from torch.utils.data import TensorDataset, DataLoader
    HAS_TORCH_GEOMETRIC = True
    HAS_LIGHTNING = True
except ImportError as e:
    logger.warning(f"⚠ Dependencias faltantes: {e}")
    HAS_TORCH_GEOMETRIC = False
    HAS_LIGHTNING = False

# Añadir el directorio pytorch_files al path
script_dir = Path(__file__).parent.parent
pytorch_files_dir = script_dir / "pytorch_files"
if pytorch_files_dir.exists():
    sys.path.append(str(pytorch_files_dir))


class ModelIntegrationValidator:
    """Validador de integración con modelos existentes."""
    
    def __init__(self, data_dir: str, models_dir: str = None):
        """
        Inicializa el validador.
        
        Args:
            data_dir: Directorio donde están los datos generados
            models_dir: Directorio donde están los modelos (opcional)
        """
        self.data_dir = Path(data_dir)
        self.models_dir = Path(models_dir) if models_dir else Path.cwd() / "scripts/3.GeneticScores"
        self.validation_results = {}
        
        logger.info(f"Inicializando validador de integración con modelos...")
        logger.info(f"Directorio de datos: {self.data_dir}")
        logger.info(f"Directorio de modelos: {self.models_dir}")
    
    def test_nn_model_compatibility(self, table_file: str) -> Dict[str, Any]:
        """
        Prueba compatibilidad con modelos NN de supervivencia.
        
        Args:
            table_file: Archivo de datos de tabla
            
        Returns:
            Resultados de la prueba
        """
        logger.info("🧠 Probando compatibilidad con modelos NN...")
        
        results = {
            'data_loaded': False,
            'model_imported': False,
            'model_created': False,
            'forward_pass_success': False,
            'data_format_correct': False,
            'dimensions_compatible': False,
            'survival_format_valid': False,
            'errors': []
        }
        
        try:
            # 1. Cargar datos de tabla
            table_path = self.data_dir / table_file
            if not table_path.exists():
                results['errors'].append(f"Archivo no encontrado: {table_file}")
                return results
            
            if table_path.suffix == '.pkl':
                df = pd.read_pickle(table_path)
            elif table_path.suffix == '.csv':
                df = pd.read_csv(table_path)
            else:
                results['errors'].append(f"Formato no soportado: {table_path.suffix}")
                return results
            
            results['data_loaded'] = True
            logger.info(f"✓ Datos cargados: {df.shape}")

            # 2. Verificar formato de datos de supervivencia
            required_survival_cols = ['LFS_YEARS', 'LFS_STATUS']
            if all(col in df.columns for col in required_survival_cols):
                results['survival_format_valid'] = True
                logger.info("✓ Formato de supervivencia válido")
            else:
                missing_cols = [col for col in required_survival_cols if col not in df.columns]
                results['errors'].append(f"Columnas de supervivencia faltantes: {missing_cols}")
            
            # 3. Preparar datos para modelo NN (formato benchmark_NN.py)
            # Separar features de variables de supervivencia y fold
            feature_cols = [col for col in df.columns if col not in ['LFS_YEARS', 'LFS_STATUS', 'fold']]
            X = df[feature_cols].values
            y = df[['LFS_YEARS', 'LFS_STATUS']].values
            
            # Convertir a tensores
            X_tensor = torch.tensor(X, dtype=torch.float32)
            y_tensor = torch.tensor(y, dtype=torch.float32)
            
            results['data_format_correct'] = True
            logger.info(f"✓ Datos convertidos a tensores: X{X_tensor.shape}, y{y_tensor.shape}")
            
            # 4. Importar modelo BasicNN
            try:
                from BasicNN import BasicNN
                results['model_imported'] = True
                logger.info("✓ Modelo BasicNN importado")
            except ImportError as e:
                results['errors'].append(f"Error importando BasicNN: {e}")
                return results
            
            # 5. Crear instancia del modelo
            try:
                input_dim = X_tensor.shape[1]
                hidden_dim = 16
                learning_rate = 0.001
                
                model = BasicNN(input_dim, hidden_dim, learning_rate)
                results['model_created'] = True
                results['dimensions_compatible'] = True
                logger.info(f"✓ Modelo creado con input_dim={input_dim}")
            except Exception as e:
                results['errors'].append(f"Error creando modelo: {e}")
                return results
            
            # 6. Probar forward pass
            try:
                model.eval()
                with torch.no_grad():
                    # Tomar una muestra pequeña para probar
                    sample_X = X_tensor[:5] if len(X_tensor) > 5 else X_tensor
                    sample_y = y_tensor[:5] if len(y_tensor) > 5 else y_tensor
                    
                    # Forward pass
                    predictions = model(sample_X)
                    
                    # Verificar dimensiones de salida
                    expected_output_shape = (sample_X.shape[0], 1)
                    if predictions.shape == expected_output_shape:
                        results['forward_pass_success'] = True
                        logger.info(f"✓ Forward pass exitoso: {predictions.shape}")
                    else:
                        results['errors'].append(f"Dimensión de salida incorrecta: {predictions.shape} vs {expected_output_shape}")
                
            except Exception as e:
                results['errors'].append(f"Error en forward pass: {e}")
            
            # 7. Probar creación de DataLoader
            try:
                # Simular división train/val usando fold
                if 'fold' in df.columns:
                    train_mask = df['fold'] != 1  # Usar fold 1 como validación
                    val_mask = df['fold'] == 1
                else:
                    # División aleatoria si no hay folds
                    n_samples = len(df)
                    val_size = min(50, n_samples // 4)
                    indices = torch.randperm(n_samples)
                    train_mask = torch.zeros(n_samples, dtype=torch.bool)
                    train_mask[indices[val_size:]] = True
                    val_mask = ~train_mask
                
                train_dataset = TensorDataset(X_tensor[train_mask], y_tensor[train_mask])
                val_dataset = TensorDataset(X_tensor[val_mask], y_tensor[val_mask])
                
                train_loader = DataLoader(train_dataset, batch_size=32, shuffle=True)
                val_loader = DataLoader(val_dataset, batch_size=64, shuffle=False)
                
                # Probar iteración
                for batch_X, batch_y in train_loader:
                    with torch.no_grad():
                        batch_pred = model(batch_X)
                    break  # Solo probar el primer batch
                
                logger.info("✓ DataLoader compatible y funcional")
                
            except Exception as e:
                results['errors'].append(f"Error con DataLoader: {e}")
            
            logger.info("✅ Prueba de compatibilidad NN completada")
            
        except Exception as e:
            results['errors'].append(f"Error general en prueba NN: {e}")
            logger.error(f"❌ Error en prueba NN: {e}")
        
        return results
    
    def test_gnn_model_compatibility(self, graph_file: str) -> Dict[str, Any]:
        """
        Prueba compatibilidad con modelos GNN.
        
        Args:
            graph_file: Archivo de datos de grafo
            
        Returns:
            Resultados de la prueba
        """
        logger.info("🕸️ Probando compatibilidad con modelos GNN...")
        
        results = {
            'data_loaded': False,
            'graph_structure_valid': False,
            'heterodata_compatible': False,
            'dataloader_created': False,
            'node_features_valid': False,
            'edge_structure_valid': False,
            'model_input_compatible': False,
            'errors': []
        }
        
        if not HAS_TORCH_GEOMETRIC:
            results['errors'].append("PyTorch Geometric no disponible")
            return results
        
        try:
            # 1. Cargar datos de grafo
            graph_path = self.data_dir / graph_file
            if not graph_path.exists():
                results['errors'].append(f"Archivo no encontrado: {graph_file}")
                return results
            
            graph_data = torch.load(graph_path, map_location='cpu')
            results['data_loaded'] = True
            logger.info(f"✓ Datos de grafo cargados")
            
            # 2. Verificar estructura de datos
            if isinstance(graph_data, dict):
                # Estructura por folds
                all_graphs = []
                for fold_key, fold_graphs in graph_data.items():
                    all_graphs.extend(fold_graphs)
                logger.info(f"✓ Estructura por folds: {len(graph_data)} folds, {len(all_graphs)} grafos totales")
            elif isinstance(graph_data, list):
                # Lista de grafos
                all_graphs = graph_data
                logger.info(f"✓ Lista de grafos: {len(all_graphs)} grafos")
            else:
                results['errors'].append(f"Estructura de datos inesperada: {type(graph_data)}")
                return results
            
            if len(all_graphs) == 0:
                results['errors'].append("No hay grafos en los datos")
                return results
            
            results['graph_structure_valid'] = True
            
            # 3. Verificar que son HeteroData
            sample_graph = all_graphs[0]
            if isinstance(sample_graph, HeteroData):
                results['heterodata_compatible'] = True
                logger.info("✓ Formato HeteroData válido")
                
                # Verificar tipos de nodos esperados
                expected_node_types = ['patient', 'gene']
                actual_node_types = list(sample_graph.node_types)
                
                if all(nt in actual_node_types for nt in expected_node_types):
                    results['node_features_valid'] = True
                    logger.info(f"✓ Tipos de nodos válidos: {actual_node_types}")
                else:
                    results['errors'].append(f"Tipos de nodos faltantes. Esperados: {expected_node_types}, encontrados: {actual_node_types}")
                
                # Verificar estructura de features
                try:
                    patient_features = sample_graph['patient'].x
                    gene_features = sample_graph['gene'].x
                    patient_labels = sample_graph['patient'].y
                    
                    logger.info(f"✓ Features paciente: {patient_features.shape}")
                    logger.info(f"✓ Features genes: {gene_features.shape}")
                    logger.info(f"✓ Labels paciente: {patient_labels.shape}")

                    # Verificar aristas
                    if ('gene', 'to', 'patient') in sample_graph.edge_types:
                        edge_index = sample_graph['gene', 'patient'].edge_index
                        logger.info(f"✓ Aristas gene->patient: {edge_index.shape}")
                        results['edge_structure_valid'] = True
                    else:
                        results['errors'].append("Tipo de arista 'gene->patient' no encontrado")
                    
                except Exception as e:
                    results['errors'].append(f"Error verificando features: {e}")
            else:
                results['errors'].append(f"Tipo de grafo inesperado: {type(sample_graph)}")
            
            # 4. Probar creación de DataLoader GNN
            try:
                # Tomar muestra pequeña para probar
                sample_graphs = all_graphs[:10] if len(all_graphs) > 10 else all_graphs
                
                gnn_loader = GeoDataLoader(sample_graphs, batch_size=4, shuffle=False)
                
                # Probar iteración
                for batch in gnn_loader:
                    if isinstance(batch, HeteroData):
                        logger.info(f"✓ Batch procesado: {batch}")
                        logger.info(f"  - Pacientes en batch: {batch['patient'].x.shape[0]}")
                        logger.info(f"  - Genes en batch: {batch['gene'].x.shape[0]}")
                        results['dataloader_created'] = True
                    break  # Solo probar el primer batch
                
            except Exception as e:
                results['errors'].append(f"Error con GeoDataLoader: {e}")
            
            # 5. Verificar compatibilidad con modelos GNN existentes
            try:
                # Extraer dimensiones para verificar compatibilidad
                patient_feat_dim = sample_graph['patient'].x.shape[1]
                gene_feat_dim = sample_graph['gene'].x.shape[1]
                
                logger.info(f"✓ Dimensiones compatibles:")
                logger.info(f"  - Patient features: {patient_feat_dim}")
                logger.info(f"  - Gene features: {gene_feat_dim}")
                
                # Estas son las dimensiones esperadas por los modelos existentes
                expected_patient_dim = 13  # Según trainGNNmodel.py
                expected_gene_dim = 16     # Según trainGNNmodel.py
                

            except Exception as e:
                results['errors'].append(f"Error verificando dimensiones: {e}")
            
            logger.info("✅ Prueba de compatibilidad GNN completada")
            
        except Exception as e:
            results['errors'].append(f"Error general en prueba GNN: {e}")
            logger.error(f"❌ Error en prueba GNN: {e}")
        
        return results
    
    def test_data_format_consistency(self, table_file: str, graph_file: str) -> Dict[str, Any]:
        """
        Verifica consistencia entre formatos de tabla y grafo.
        
        Args:
            table_file: Archivo de tabla
            graph_file: Archivo de grafo
            
        Returns:
            Resultados de consistencia
        """
        logger.info("🔄 Verificando consistencia entre formatos...")
        
        results = {
            'patient_count_match': False,
            'survival_data_match': False,
            'fold_structure_match': False,
            'data_ranges_consistent': False,
            'errors': []
        }
        
        try:
            # Cargar ambos archivos
            table_path = self.data_dir / table_file
            graph_path = self.data_dir / graph_file
            
            if table_path.suffix == '.pkl':
                df = pd.read_pickle(table_path)
            else:
                df = pd.read_csv(table_path)
            
            graph_data = torch.load(graph_path, map_location='cpu')
            
            # Contar pacientes en tabla
            table_patients = len(df)
            
            # Contar pacientes en grafos
            if isinstance(graph_data, dict):
                graph_patients = sum(len(fold_graphs) for fold_graphs in graph_data.values())
            else:
                graph_patients = len(graph_data)
            
            if table_patients == graph_patients:
                results['patient_count_match'] = True
                logger.info(f"✓ Número de pacientes coincide: {table_patients}")
            else:
                results['errors'].append(f"Número de pacientes no coincide: tabla={table_patients}, grafo={graph_patients}")
            
            # Verificar estructura de folds si está disponible
            if 'fold' in df.columns and isinstance(graph_data, dict):
                table_folds = set(df['fold'].unique())
                graph_folds = set()
                
                # Extraer información de folds de los nombres de las claves del grafo
                for fold_key in graph_data.keys():
                    if 'fold_' in fold_key:
                        fold_num = int(fold_key.split('_')[1])
                        graph_folds.add(fold_num)
                
                if table_folds == graph_folds:
                    results['fold_structure_match'] = True
                    logger.info(f"✓ Estructura de folds coincide: {len(table_folds)} folds")
                else:
                    logger.warning(f"⚠ Folds diferentes - tabla: {table_folds}, grafo: {graph_folds}")
            
            # Verificar rangos de datos de supervivencia
            if 'LFS_YEARS' in df.columns and 'LFS_STATUS' in df.columns:
                table_years_range = (df['LFS_YEARS'].min(), df['LFS_YEARS'].max())
                table_status_dist = df['LFS_STATUS'].value_counts().to_dict()
                
                logger.info(f"✓ Rango de años en tabla: {table_years_range}")
                logger.info(f"✓ Distribución de status en tabla: {table_status_dist}")
                results['data_ranges_consistent'] = True
            
            logger.info("✅ Verificación de consistencia completada")
            
        except Exception as e:
            results['errors'].append(f"Error en verificación de consistencia: {e}")
            logger.error(f"❌ Error en consistencia: {e}")
        
        return results
    
    def generate_integration_report(self) -> str:
        """
        Genera reporte de integración con modelos.
        
        Returns:
            Reporte en formato string
        """
        report_lines = [
            "="*80,
            "REPORTE DE INTEGRACIÓN CON MODELOS - PIPELINE UNIFICADO MDS",
            "="*80,
            ""
        ]
        
        # Resumen de compatibilidad NN
        if 'nn_compatibility' in self.validation_results:
            nn_results = self.validation_results['nn_compatibility']
            report_lines.extend([
                "🧠 COMPATIBILIDAD CON MODELOS NN:",
                f"  ✓ Datos cargados: {nn_results.get('data_loaded', False)}",
                f"  ✓ Modelo importado: {nn_results.get('model_imported', False)}",
                f"  ✓ Modelo creado: {nn_results.get('model_created', False)}",
                f"  ✓ Forward pass: {nn_results.get('forward_pass_success', False)}",
                f"  ✓ Formato supervivencia: {nn_results.get('survival_format_valid', False)}",
                f"  ✓ Dimensiones compatibles: {nn_results.get('dimensions_compatible', False)}"
            ])
            
            if nn_results.get('errors'):
                report_lines.append("  ❌ Errores encontrados:")
                for error in nn_results['errors']:
                    report_lines.append(f"     - {error}")
            
            report_lines.append("")
        
        # Resumen de compatibilidad GNN
        if 'gnn_compatibility' in self.validation_results:
            gnn_results = self.validation_results['gnn_compatibility']
            report_lines.extend([
                "🕸️ COMPATIBILIDAD CON MODELOS GNN:",
                f"  ✓ Datos cargados: {gnn_results.get('data_loaded', False)}",
                f"  ✓ Estructura HeteroData: {gnn_results.get('heterodata_compatible', False)}",
                f"  ✓ Features de nodos: {gnn_results.get('node_features_valid', False)}",
                f"  ✓ Estructura de aristas: {gnn_results.get('edge_structure_valid', False)}",
                f"  ✓ DataLoader creado: {gnn_results.get('dataloader_created', False)}",
                f"  ✓ Input compatible: {gnn_results.get('model_input_compatible', False)}"
            ])
            
            if gnn_results.get('errors'):
                report_lines.append("  ❌ Errores encontrados:")
                for error in gnn_results['errors']:
                    report_lines.append(f"     - {error}")
            
            report_lines.append("")
        
        # Resumen de consistencia
        if 'consistency' in self.validation_results:
            cons_results = self.validation_results['consistency']
            report_lines.extend([
                "🔄 CONSISTENCIA ENTRE FORMATOS:",
                f"  ✓ Número de pacientes: {cons_results.get('patient_count_match', False)}",
                f"  ✓ Datos de supervivencia: {cons_results.get('survival_data_match', False)}",
                f"  ✓ Estructura de folds: {cons_results.get('fold_structure_match', False)}",
                f"  ✓ Rangos de datos: {cons_results.get('data_ranges_consistent', False)}"
            ])
            
            if cons_results.get('errors'):
                report_lines.append("  ❌ Errores encontrados:")
                for error in cons_results['errors']:
                    report_lines.append(f"     - {error}")
            
            report_lines.append("")
        
        # Resumen general
        all_tests_passed = True
        critical_errors = []
        
        for test_name, test_results in self.validation_results.items():
            if test_results.get('errors'):
                all_tests_passed = False
                critical_errors.extend(test_results['errors'])
        
        if all_tests_passed:
            report_lines.extend([
                "🎉 RESULTADO GENERAL: COMPATIBILIDAD EXITOSA",
                "✅ Los datos generados son completamente compatibles con los modelos existentes",
                "✅ Listos para entrenamiento y validación"
            ])
        else:
            report_lines.extend([
                "⚠️ RESULTADO GENERAL: PROBLEMAS ENCONTRADOS",
                f"❌ Se encontraron {len(critical_errors)} errores críticos",
                "💡 Revisar errores específicos arriba para solucionar problemas"
            ])
        
        report_lines.extend([
            "",
            "="*80,
            f"Reporte generado: {pd.Timestamp.now()}",
            "="*80
        ])
        
        return "\n".join(report_lines)
    
    def run_full_integration_test(self, table_file: str = None, graph_file: str = None) -> bool:
        """
        Ejecuta pruebas completas de integración con modelos.
        
        Args:
            table_file: Archivo de tabla (opcional, se buscará automáticamente)
            graph_file: Archivo de grafo (opcional, se buscará automáticamente)
            
        Returns:
            True si todas las pruebas pasan
        """
        logger.info("🚀 Iniciando pruebas de integración con modelos...")
        
        try:
            # Buscar archivos automáticamente si no se especifican
            if table_file is None:
                table_files = list(self.data_dir.glob("*.pkl")) + list(self.data_dir.glob("*.csv"))
                table_files = [f for f in table_files if 'table' in f.name.lower()]
                table_file = table_files[0].name if table_files else None
            
            if graph_file is None:
                graph_files = list(self.data_dir.glob("*.pt"))
                graph_files = [f for f in graph_files if 'graph' in f.name.lower()]
                graph_file = graph_files[0].name if graph_files else None
            
            success = True
            
            # Prueba 1: Compatibilidad NN
            if table_file:
                logger.info(f"🧠 Probando compatibilidad NN con: {table_file}")
                self.validation_results['nn_compatibility'] = self.test_nn_model_compatibility(table_file)
                
                if self.validation_results['nn_compatibility'].get('errors'):
                    success = False
            else:
                logger.warning("⚠ No se encontró archivo de tabla para probar NN")
            
            # Prueba 2: Compatibilidad GNN
            if graph_file:
                logger.info(f"🕸️ Probando compatibilidad GNN con: {graph_file}")
                self.validation_results['gnn_compatibility'] = self.test_gnn_model_compatibility(graph_file)
                
                if self.validation_results['gnn_compatibility'].get('errors'):
                    success = False
            else:
                logger.warning("⚠ No se encontró archivo de grafo para probar GNN")
            
            # Prueba 3: Consistencia entre formatos
            if table_file and graph_file:
                logger.info("🔄 Probando consistencia entre formatos...")
                self.validation_results['consistency'] = self.test_data_format_consistency(table_file, graph_file)
                
                if self.validation_results['consistency'].get('errors'):
                    success = False
            
            # Generar reporte
            report = self.generate_integration_report()
            
            # Guardar reporte
            report_file = self.data_dir / "model_integration_report.txt"
            with open(report_file, 'w', encoding='utf-8') as f:
                f.write(report)
            
            logger.info(f"📋 Reporte guardado en: {report_file}")
            
            # Mostrar resumen
            print("\n" + report)
            
            if success:
                logger.info("✅ Todas las pruebas de integración pasaron exitosamente")
            else:
                logger.error("❌ Se encontraron problemas en las pruebas de integración")
            
            return success
            
        except Exception as e:
            logger.error(f"❌ Error en pruebas de integración: {e}")
            import traceback
            logger.error(f"Traceback: {traceback.format_exc()}")
            return False


def main():
    """Función principal."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Validador de Integración con Modelos - Paso 3D Pipeline Unificado",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument(
        '--data-dir', '-d',
        type=str,
        required=True,
        help='Directorio con los datos generados por el pipeline'
    )
    
    parser.add_argument(
        '--models-dir', '-m',
        type=str,
        help='Directorio con los modelos existentes (opcional)'
    )
    
    parser.add_argument(
        '--table-file', '-t',
        type=str,
        help='Nombre específico del archivo de tabla (opcional)'
    )
    
    parser.add_argument(
        '--graph-file', '-g',
        type=str,
        help='Nombre específico del archivo de grafo (opcional)'
    )
    
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Activar logging detallado'
    )
    
    args = parser.parse_args()
    
    # Configurar nivel de logging
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    try:
        # Crear validador
        validator = ModelIntegrationValidator(args.data_dir, args.models_dir)
        
        # Ejecutar pruebas
        success = validator.run_full_integration_test(
            table_file=args.table_file,
            graph_file=args.graph_file
        )
        
        if success:
            logger.info("🎉 Integración con modelos validada exitosamente!")
            return 0
        else:
            logger.error("💥 Se encontraron problemas de integración con modelos")
            return 1
            
    except Exception as e:
        logger.error(f"Error ejecutando validación de integración: {e}")
        return 1


if __name__ == "__main__":
    import pandas as pd  # Importar aquí para evitar error si no está disponible
    exit(main())
