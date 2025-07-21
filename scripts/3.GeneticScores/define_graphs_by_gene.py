from pathlib import Path
import yaml
import pickle
import logging
from typing import Dict, Any, Optional, List
import pandas as pd
import numpy as np
from itertools import combinations
import argparse
import torch
# Añadir el directorio utils al path para importar módulos
import sys
sys.path.append(str(Path(__file__).parent / "utils"))

from data_loader import DataLoader
from data_generators import DataGenerators
from clinical_processor import ClinicalProcessor
from gnn_inference_utils import load_config

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


class GraphGenePipeline:
    
    def __init__(self, config_path: str, output_dir: str):
        """
        Inicializar pipeline con configuración.
        
        Args:
            config_path: Ruta al archivo de configuración YAML/JSON
            output_dir: Directorio donde guardar los resultados
        """
        self.config_path = Path(config_path)
        self.output_dir = Path(output_dir)
        self.config = load_config(str(self.config_path))
        
        # Crear directorio de salida si no existe
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Inicializar componentes
        self.data_loader = DataLoader(self.config)
        self.clinical_processor = ClinicalProcessor(self.config)
        self.data_generators = DataGenerators(self.config)
        self.data_generators.fitted_transformer =  pickle.load(open(self.config["data_paths"]["transformer_file"], "rb"))

        logger.info(f"Pipeline inicializado con configuración: {config_path}")
        logger.info(f"Directorio de salida: {output_dir}")
    
   
    def run_pipeline(self) -> Dict[str, Any]:
        """Ejecutar pipeline completo de preparación de datos."""
        logger.info("=== Iniciando pipeline de preparación de datos ===")
        
        try:
            # 1. Cargar datos
            logger.info("Paso 1: Cargando datos...")
            df_clinical = self.data_loader.load_clinical_data()
            df_mutations = self.data_loader.load_mutation_vaf_data()
            gene_embeddings = self.data_loader.load_gene_embeddings()
            
            ## Seleccionar casos completos:
            df_processed, metadata = self.clinical_processor._filter_complete_cases_from_config(df_clinical)
        

            # 2. Seleccionar genes
            logger.info("Paso 2: Seleccionando genes...")
            selected_genes = self.config['variable_processing']['gene_selection']

            logger.info("Paso 3: Generando grafos...")
            graph_list, gene_dict = self.generate_graph_data_individual(df_processed, self.config['variable_processing'],
                                                df_mutations, gene_embeddings)


            # Guardar archivo
            output_file = self.output_dir / f"{self.config['data_generation']['output_filename']}.pt"
            torch.save((graph_list, gene_dict), output_file)
        except Exception as e:
            logger.error(f"Error en la ejecución de la pipeline: {e}")
            raise

    def generate_graph_data_individual(self, df_processed: pd.DataFrame, metadata: Dict[str, Any],
            df_mutations: Optional[pd.DataFrame], gene_embeddings: Optional[Dict[str, np.ndarray]]) -> Dict[str, Any]:
        """Genera datos en formato grafo usando HeteroData."""
        logger.info("Generando datos en formato grafo...")
        
        # Preparar datos de genes
        selected_genes = metadata['gene_selection']

        if len(selected_genes) == 0:
            selected_genes = df_mutations['GENE'].unique().tolist()
        gene_features = self.data_generators._prepare_gene_features(selected_genes, gene_embeddings)

        
        graph_list = []
        gene_dict = {}
        idx = 0
        for patient_idx in range(len(df_processed)):
            graph = self.data_generators._create_hetero_graph_for_patient(
                df_processed, patient_idx, gene_features, selected_genes, df_mutations
            )
            graph_list.append(graph)

            ## Save genes present in the graph            
            patient_id = df_processed['ID'].iloc[patient_idx]
            patient_mutations = df_mutations[df_mutations['ID'] == patient_id]
            if len(patient_mutations) > 0:
                genes = patient_mutations['GENE'].unique()
                for gene in genes:
                    if gene not in gene_dict:
                        gene_dict[gene] = [idx]
                    else:
                        gene_dict[gene].append(idx)
            idx += 1

        return graph_list, gene_dict


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
    
    args = parser.parse_args()
    

    try:
        # Crear pipeline
        pipeline = GraphGenePipeline(args.config, args.output)
        

        # Ejecutar pipeline
        pipeline.run_pipeline()

        
        logger.info("Pipeline ejecutado exitosamente!")
        logger.info(f"Resultados disponibles en: {args.output}")
        
        return 0
        
    except Exception as e:
        logger.error(f"Error ejecutando pipeline: {e}")
        return 1


if __name__ == "__main__":
    exit(main())
