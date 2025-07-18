#!/usr/bin/env python3
"""
Script para ejecutar cross-validation con modelos GNN.
=============
Este script ejecuta un análisis de cross-validation (CV) utilizando modelos de Graph Neural Networks (GNN) sobre datos de grafos preprocesados y particionados en folds. Es el equivalente en Python al script recomputeIPSSM_cv.R, pero adaptado para modelos GNN y el pipeline MDS_subtypes.
Características principales:
----------------------------
- Valida la existencia y formato de los archivos de entrada requeridos (datos de grafos y configuración).
- Permite seleccionar modelos específicos a ejecutar mediante argumentos de línea de comandos.
- Ejecuta el pipeline de cross-validation usando funciones definidas en `gnn_cv_functions.py`.
- Muestra y guarda los resultados de CV, incluyendo métricas como el C-index y dimensiones de los datos.
- Soporta ejecución tanto en entornos Docker como localmente.
- Incluye opciones para solo validar entradas sin ejecutar el análisis.
----
Argumentos principales:
-----------------------
--config_path         Ruta al archivo YAML de configuración de modelos GNN (opcional, por defecto: scripts/3.GeneticScores/configs/gnn_models_config.yaml).
--graph_data_path     Ruta al archivo .pt con los datos de grafos y folds (obligatorio).
--results_dir         Directorio donde se guardarán los resultados (obligatorio).
--lightning_logs_dir  Directorio para guardar logs de Lightning (opcional).
--models              Lista de modelos a ejecutar (deben estar definidos en el YAML, opcional).
--validate_only       Solo valida las entradas, sin ejecutar el análisis (flag).

Ejemplo de uso desde la línea de comando:
-----------------------------------------

python scripts/3.GeneticScores/run_gnn_cv.py \
    --graph_data_path results/gnn/preprocess/graphs_genesEncScGPT_folds.pt \
    --config_path scripts/3.GeneticScores/configs/gnn_models_config.yaml \
    --results_dir results/gnn/cv \
    --lightning_logs_dir results/gnn/lightning_logs \
    --models v1_boolean v4_boolean_optim

Autores:
--------
Adaptado del pipeline MDS_subtypes.
"""

import os
import sys
import logging
import argparse
from pathlib import Path
from typing import Optional, List

# Añadir paths necesarios
sys.path.append('scripts/3.GeneticScores/utils')
from gnn_cv_functions import run_gnn_cv_pipeline, GNNCrossValidator

# Configurar logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def setup_environment():
    """Configura el entorno de ejecución."""
    # Verificar que estamos en el directorio correcto
    if not Path("scripts/3.GeneticScores").exists():
        raise RuntimeError("Ejecutar desde el directorio raíz del proyecto")
    
    # Verificar CUDA si está disponible
    try:
        import torch
        if torch.cuda.is_available():
            logger.info(f"CUDA disponible: {torch.cuda.get_device_name(0)}")
        else:
            logger.info("Ejecutando en CPU")
    except ImportError:
        logger.warning("PyTorch no disponible")


def validate_inputs(graph_data_path: str) -> bool:
    """
    Valida que los archivos de entrada existan.
    
    Args:
        graph_data_path: Ruta al archivo de datos de grafos
        
    Returns:
        True si las validaciones pasan
    """
    if not os.path.exists(graph_data_path):
        logger.error(f"Archivo de grafos no encontrado: {graph_data_path}")
        return False
    
    try:
        import torch
        data = torch.load(graph_data_path)
        if not isinstance(data, dict):
            logger.error("Formato de datos de grafos inválido")
            return False
        
        n_folds = len(data)
        logger.info(f"Archivo de grafos válido: {n_folds} folds")
        
        # Verificar que cada fold tenga datos
        for fold_name, fold_data in data.items():
            if not isinstance(fold_data, list) or len(fold_data) == 0:
                logger.error(f"Fold {fold_name} vacío o inválido")
                return False
        
        return True
        
    except Exception as e:
        logger.error(f"Error validando archivo de grafos: {str(e)}")
        return False


def run_cv_analysis(graph_data_path: str,
                   config_path: str,
                   results_dir: str,
                   lightning_logs_dir: Optional[str] = None,
                   models_to_run: Optional[List[str]] = None) -> bool:
    """
    Ejecuta el análisis de cross-validation.
    
    Args:
        graph_data_path: Ruta a los datos de grafos
        config_path: Ruta al archivo de configuración YAML
        results_dir: Directorio para resultados
        models_to_run: Lista de modelos a ejecutar
        
    Returns:
        True si el análisis se completa exitosamente
    """
    logger.info("="*60)
    logger.info("INICIANDO CROSS-VALIDATION DE MODELOS GNN")
    logger.info("="*60)
    
    try:

        logger.info(f"Modelos a ejecutar: {models_to_run}")
        logger.info(f"Configuración desde: {config_path}")
        
        # Ejecutar CV
        results = run_gnn_cv_pipeline(
            graph_data_path=graph_data_path,
            config_path=config_path,
            results_dir=results_dir,
            lightning_logs_dir=lightning_logs_dir,
            models_to_run=models_to_run
        )
        
        # Mostrar resultados
        logger.info("\n" + "="*60)
        logger.info("RESULTADOS DE CROSS-VALIDATION")
        logger.info("="*60)
        
        for model_name, result in results.items():
            c_idx = result.get('mean_val_c_index', 'N/A')
            c_std = result.get('std_val_c_index', 'N/A')
            data_dims = result.get('data_dimensions', {})
            
            if c_idx != 'N/A' and c_std != 'N/A':
                logger.info(f"{model_name:20s} | C-index: {c_idx:.4f} ± {c_std:.4f} | Dims: {data_dims.get('patient_feat_dim', '?')}x{data_dims.get('gene_feat_dim', '?')}")
            else:
                logger.info(f"{model_name:20s} | C-index: Error en cálculo")
        
        logger.info("="*60)
        logger.info(f"Resultados guardados en: {results_dir}")
        
        return True
        
    except Exception as e:
        logger.error(f"Error durante el análisis de CV: {str(e)}")
        return False


def main():
    """Función principal."""
    parser = argparse.ArgumentParser(
        description="Ejecutar cross-validation con modelos GNN",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        "--config_path",
        type=str,
        default="scripts/3.GeneticScores/configs/gnn_models_config.yaml",
        help="Ruta al archivo de configuración YAML de modelos"
    )
    
    parser.add_argument(
        "--graph_data_path",
        type=str,
        help="Ruta al archivo de datos de grafos con folds"
    )
    
    parser.add_argument(
        "--results_dir",
        type=str,
        help="Directorio para guardar resultados"
    )
    
    parser.add_argument(
        "--lightning_logs_dir",
        type=str,
        help="Directorio para guardar logs de Lightning"
    )
    

    parser.add_argument(
        "--models",
        nargs="+",
        default=None,
        help="Modelos específicos a ejecutar (deben existir en configuración YAML)"
    )

    parser.add_argument(
        "--validate_only",
        action="store_true",
        help="Solo validar entradas sin ejecutar CV"
    )
    
    args = parser.parse_args()
    
    try:
        # Setup inicial
        setup_environment()
        
        # Validar entradas
        if not validate_inputs(args.graph_data_path):
            logger.error("Validación de entradas fallida")
            return 1
        
        # Validar archivo de configuración
        if not os.path.exists(args.config_path):
            logger.warning(f"Archivo de configuración no encontrado: {args.config_path}")
            logger.info("Se usará configuración por defecto")
        
        if args.validate_only:
            logger.info("Validación exitosa. Terminando sin ejecutar CV.")
            return 0
        
        # Crear directorio de resultados
        Path(args.results_dir).mkdir(parents=True, exist_ok=True)
        
        # Ejecutar análisis
        success = run_cv_analysis(
            graph_data_path=args.graph_data_path,
            config_path=args.config_path,
            results_dir=args.results_dir,
            lightning_logs_dir=args.lightning_logs_dir,
            models_to_run=args.models
        )
        
        if success:
            logger.info("Cross-validation completado exitosamente!")
            return 0
        else:
            logger.error("Cross-validation falló")
            return 1
            
    except KeyboardInterrupt:
        logger.info("Ejecución interrumpida por el usuario")
        return 1
    except Exception as e:
        logger.error(f"Error inesperado: {str(e)}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
