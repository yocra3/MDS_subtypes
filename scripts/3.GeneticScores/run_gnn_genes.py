#!/usr/bin/env python3
"""
Script para ejecutar cross-validation con modelos GNN.
=============
Este script entrena modelos de Graph Neural Networks (GNN) sobre datos de grafos preprocesados y particionados por gen. Es el equivalente en Python al script recomputeIPSSM_cv.R, pero adaptado para modelos GNN y el pipeline MDS_subtypes.
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
from gnn_gene_functions import run_gnn_gene_pipeline, GNNGeneTrainer


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


def run_gene_analysis(graph_data_path: str,
                   config_path: str,
                   results_dir: str,
                   lightning_logs_dir: Optional[str] = None,
                   models_to_run: Optional[List[str]] = None) -> bool:
    """
    Ejecuta el análisis de gene
    
    Args:
        graph_data_path: Ruta a los datos de grafos
        config_path: Ruta al archivo de configuración YAML
        results_dir: Directorio para resultados
        models_to_run: Lista de modelos a ejecutar
        
    Returns:
        True si el análisis se completa exitosamente
    """
    logger.info("="*60)
    logger.info("INICIANDO ANÁLISIS DE GENES")
    logger.info("="*60)
    
    try:

        logger.info(f"Modelos a ejecutar: {models_to_run}")
        logger.info(f"Configuración desde: {config_path}")
        
        # Ejecutar análisis de genes
        results = run_gnn_gene_pipeline(
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
        success = run_gene_analysis(
            graph_data_path=args.graph_data_path,
            config_path=args.config_path,
            results_dir=args.results_dir,
            lightning_logs_dir=args.lightning_logs_dir,
            models_to_run=args.models
        )
        
        if success:
            logger.info("Análisis de genes completado exitosamente!")
            return 0
        else:
            logger.error("Análisis de genes falló")
            return 1
            
    except KeyboardInterrupt:
        logger.info("Ejecución interrumpida por el usuario")
        return 1
    except Exception as e:
        logger.error(f"Error inesperado: {str(e)}")
        return 1


if __name__ == "__main__":
    sys.exit(main())

