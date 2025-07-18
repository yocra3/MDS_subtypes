#!/usr/bin/env python3
"""
Pipeline Unificado de Preparación de Datos para Modelos de Supervivencia (MDS)

Este script unifica la preparación de datos para modelos GNN y CoxPH, permitiendo:
- Selección flexible de variables clínicas
- Generación de datos para modelos de grafos (GNN) y tablas (NN)
- Uso de folds de cross-validation pre-generados
- Configuración completamente flexible vía archivos YAML/JSON

Autor: Sistema automatizado
Fecha: 2025-07-04

# Ejemplo de uso desde la línea de comandos:
#
# Ejecutar el pipeline con configuración YAML y salida en un directorio específico:
#   python prepare_data_unified.py --config configs/default.yaml --output results/experiment1
#
# Ejecutar solo generación de datos tipo grafo:
#   python prepare_data_unified.py --config configs/gnn_only.yaml --output results/gnn_test --formats graph
#
# Mostrar la configuración sin ejecutar el pipeline (modo dry-run):
#   python prepare_data_unified.py --config configs/default.yaml --output results/test --dry-run
"""

import argparse
import logging
import sys
import os
from pathlib import Path
import yaml
import json
import pandas as pd
import pickle
from typing import Dict, Any, Optional, List
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline

# Añadir el directorio utils al path para importar módulos
sys.path.append(str(Path(__file__).parent / "utils"))

from data_loader import DataLoader
from clinical_processor import ClinicalProcessor
from data_generators import DataGenerators

# Configurar logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('prepare_data_unified.log'),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)


class UnifiedDataPipeline:
    """Pipeline unificado para preparación de datos de supervivencia MDS."""
    
    def __init__(self, config_path: str, output_dir: str):
        """
        Inicializar pipeline con configuración.
        
        Args:
            config_path: Ruta al archivo de configuración YAML/JSON
            output_dir: Directorio donde guardar los resultados
        """
        self.config_path = Path(config_path)
        self.output_dir = Path(output_dir)
        self.config = self._load_config()
        
        # Crear directorio de salida si no existe
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Inicializar componentes
        self.data_loader = DataLoader(self.config)
        self.clinical_processor = ClinicalProcessor(self.config)
        self.data_generators = DataGenerators(self.config)
        
        logger.info(f"Pipeline inicializado con configuración: {config_path}")
        logger.info(f"Directorio de salida: {output_dir}")
    
    def _load_config(self) -> Dict[str, Any]:
        """Cargar configuración desde archivo YAML o JSON."""
        try:
            with open(self.config_path, 'r', encoding='utf-8') as f:
                if self.config_path.suffix.lower() == '.yaml' or self.config_path.suffix.lower() == '.yml':
                    config = yaml.safe_load(f)
                elif self.config_path.suffix.lower() == '.json':
                    config = json.load(f)
                else:
                    raise ValueError(f"Formato de configuración no soportado: {self.config_path.suffix}")
            
            logger.info(f"Configuración cargada exitosamente desde {self.config_path}")
            return config
            
        except Exception as e:
            logger.error(f"Error cargando configuración desde {self.config_path}: {e}")
            raise
    
    def run_pipeline(self, output_formats: Optional[List[str]] = None) -> Dict[str, Any]:
        """Ejecutar pipeline completo de preparación de datos."""
        logger.info("=== Iniciando pipeline de preparación de datos ===")
        
        try:
            # 1. Cargar datos
            logger.info("Paso 1: Cargando datos...")
            df_clinical = self.data_loader.load_clinical_data()
            df_mutations = self.data_loader.load_mutation_vaf_data()
            gene_embeddings = self.data_loader.load_gene_embeddings()
            
            # 2. Seleccionar genes
            logger.info("Paso 2: Seleccionando genes...")
            selected_genes = self.config['variable_processing']['gene_selection']

            # 3. Procesar variables clínicas (incluye filtro obligatorio de casos completos)
            logger.info("Paso 3: Procesando variables clínicas...")
            df_processed, processing_metadata, fitted_transformer = self.clinical_processor.prepare_clinical_features(df_clinical, selected_genes)
            self.data_generators.fitted_transformer = fitted_transformer

            # Verificar que el filtro de casos completos fue aplicado
            if 'complete_cases_filtering' not in processing_metadata:
                raise ValueError("El filtro de casos completos no fue aplicado correctamente")
            
            # Log resumen del filtro de casos completos
            filtering_info = processing_metadata['complete_cases_filtering']
            logger.info(f"Filtro de casos completos aplicado: {filtering_info['initial_samples']} → {filtering_info['final_samples']} muestras")
            
            # 4. Generar datos en formatos requeridos
            output_formats = output_formats or self.config['data_generation']['output_formats']
            logger.info(f"Paso 4: Generando datos en formatos: {output_formats}")
            
            results = {}
            for format_type in output_formats:
                logger.info(f"Generando formato: {format_type}")
                
                if format_type == 'table':
                    format_results = self.data_generators.generate_table_data(
                        df_processed, processing_metadata, self.output_dir
                    )
                elif format_type == 'graph':
                                  
                    if len(processing_metadata['selected_genes']) == 0:
                        processing_metadata['selected_genes'] = df_mutations['GENE'].unique().tolist()
                    
                    format_results = self.data_generators.generate_graph_data(
                        df_processed, processing_metadata, df_mutations, gene_embeddings, self.output_dir
                    )
                else:
                    logger.warning(f"Formato desconocido: {format_type}")
                    continue
                
                results[format_type] = format_results
            
            # 5. Guardar metadatos del pipeline completo
            logger.info("Paso 5: Guardando metadatos del pipeline...")
            self._save_pipeline_metadata(results, processing_metadata)
            
            # 6. Guardar fitted transformers para aplicar a nuevos datos
            logger.info("Paso 6: Guardando fitted transformers...")
            self._save_fitted_transformers()
            
            logger.info("=== Pipeline completado exitosamente ===")
            return results
            
        except Exception as e:
            logger.error(f"Error en el pipeline: {e}")
            raise
    
    def _save_pipeline_metadata(self, results: Dict[str, Any], processing_metadata: Dict[str, Any]) -> None:
        """Guardar metadatos completos del pipeline."""
        metadata = {
            'pipeline_info': {
                'version': '1.0.0',
                'execution_timestamp': str(pd.Timestamp.now()),
                'config_file': str(self.config_path),
                'output_directory': str(self.output_dir)
            },
            'config_used': self.config,
            'processing_metadata': processing_metadata,
            'generated_files': results,
            'complete_cases_filtering': processing_metadata.get('complete_cases_filtering', {}),
            'validation_info': {
                'config_driven_filtering': True,
                'mandatory_filtering_applied': 'complete_cases_filtering' in processing_metadata,
                'variables_from_config': processing_metadata.get('complete_cases_filtering', {}).get('config_source', {})
            }
        }
        
        # Guardar metadatos principales
        metadata_path = self.output_dir / 'pipeline_metadata.json'
        with open(metadata_path, 'w', encoding='utf-8') as f:
            json.dump(metadata, f, indent=2, ensure_ascii=False, default=str)
        
        logger.info(f"Metadatos del pipeline guardados en: {metadata_path}")
        
        # Guardar resumen del filtro de casos completos
        filtering_summary_path = self.output_dir / 'complete_cases_filtering_summary.json'
        filtering_summary = {
            'summary': processing_metadata.get('complete_cases_filtering', {}),
            'timestamp': str(pd.Timestamp.now()),
            'config_source': str(self.config_path)
        }
        
        with open(filtering_summary_path, 'w', encoding='utf-8') as f:
            json.dump(filtering_summary, f, indent=2, ensure_ascii=False, default=str)
        
        logger.info(f"Resumen del filtro de casos completos guardado en: {filtering_summary_path}")
        
        # Guardar copia de la configuración utilizada
        config_copy_path = self.output_dir / 'config_used.yaml'
        with open(config_copy_path, 'w', encoding='utf-8') as f:
            yaml.dump(self.config, f, default_flow_style=False, allow_unicode=True)
        
        logger.info(f"Copia de la configuración guardada en: {config_copy_path}")
        
        # Log resumen para el usuario
        filtering_info = processing_metadata.get('complete_cases_filtering', {})
        logger.info("=== RESUMEN DEL FILTRO DE CASOS COMPLETOS ===")
        logger.info(f"Muestras iniciales: {filtering_info.get('initial_samples', 'N/A')}")
        logger.info(f"Muestras finales: {filtering_info.get('final_samples', 'N/A')}")
        logger.info(f"Muestras eliminadas: {filtering_info.get('removed_samples', 'N/A')}")
        logger.info(f"Tasa de eliminación: {filtering_info.get('removal_rate', 0):.1%}")
        logger.info(f"Variables utilizadas: {len(filtering_info.get('model_variables_used', []))}")
        logger.info("===============================================")

    def _save_fitted_transformers(self) -> None:
        """Guardar los fitted transformers del procesador clínico para aplicar a nuevos datos."""
        try:
            fitted_transformers = self.data_generators.fitted_transformer
            
            if not fitted_transformers:
                logger.info("No hay fitted transformers para guardar")
                return
            
            transformers_dir = self.output_dir / 'fitted_transformers'
            transformers_dir.mkdir(parents=True, exist_ok=True)
            
            pipeline_path = transformers_dir / 'preprocessing_pipeline.pkl'
            with open(pipeline_path, 'wb') as f:
                pickle.dump(fitted_transformers, f)
            logger.info(f"Pipeline unificada guardada en: {pipeline_path}")
            
                            
        except Exception as e:
            logger.error(f"Error guardando fitted transformers: {e}")
            raise



def main():
    """Función principal con interfaz de línea de comandos."""
    parser = argparse.ArgumentParser(
        description="Pipeline Unificado de Preparación de Datos para Modelos MDS",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Ejemplos de uso:
  python prepare_data_unified.py --config configs/default.yaml --output results/experiment1
  python prepare_data_unified.py --config configs/gnn_only.yaml --output results/gnn_test --formats graph
  python prepare_data_unified.py --config configs/full_vars.yaml --output results/full_analysis
        """
    )
    
    parser.add_argument(
        '--config', '-c',
        type=str,
        required=True,
        help='Ruta al archivo de configuración (YAML o JSON)'
    )
    
    parser.add_argument(
        '--output', '-o',
        type=str,
        required=True,
        help='Directorio donde guardar los resultados'
    )
    
    parser.add_argument(
        '--formats', '-f',
        nargs='+',
        choices=['table', 'graph'],
        help='Formatos de salida a generar (por defecto usa la configuración)'
    )
    
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Activar logging detallado'
    )
    
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Mostrar configuración sin ejecutar el pipeline'
    )
    
    args = parser.parse_args()
    
    # Configurar nivel de logging
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    try:
        # Crear pipeline
        pipeline = UnifiedDataPipeline(args.config, args.output)
        
        if args.dry_run:
            logger.info("=== MODO DRY-RUN: Mostrando configuración ===")
            print(json.dumps(pipeline.config, indent=2, ensure_ascii=False))
            return
        
        # Ejecutar pipeline
        results = pipeline.run_pipeline(output_formats=args.formats)
        
        logger.info("Pipeline ejecutado exitosamente!")
        logger.info(f"Resultados disponibles en: {args.output}")
        
        return 0
        
    except Exception as e:
        logger.error(f"Error ejecutando pipeline: {e}")
        return 1


if __name__ == "__main__":
    exit(main())
