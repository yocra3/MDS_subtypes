## Evaluate effect of genes in GNN models
import sys
import torch
import yaml

# Agregar ruta para importar clases de evaluación GNN
sys.path.append('./scripts/3.GeneticScores/utils/')
from evalGNN_classes import EvalGNN
from gnn_inference_utils import find_checkpoint, instantiate_gnn_model
from typing import Dict, Any, Tuple
from torch_geometric.data.batch import Batch
import numpy as np
import pandas as pd
import argparse
import os

class GNN_Inference_individual:
    """
    Clase para ejecutar inferencia de modelos GNN con validación cruzada.
    
    Esta clase permite cargar datos de grafos, configuraciones de modelos y checkpoints
    para realizar inferencia en múltiples folds de validación cruzada.
    """
        
    def __init__(self, 
        graph_data_path: str,
        config_path: str = "scripts/3.GeneticScores/configs/gnn_models_config.yaml",
        check_point_path: str = "results/gnn/full_input_model/model_checkpoints",
        results_dir: str = "results/gnn/cv_results"):
        """
        Inicializa la clase GNN_Inference_individual.
        
        Args:
            graph_data_path (str): Ruta al archivo con los datos de grafos (.pt)
            config_path (str): Ruta al archivo de configuración YAML
            check_point_path (str): Ruta al directorio con los checkpoints de modelos
            results_dir (str): Directorio donde guardar los resultados
        """

        self.graph_data_path = graph_data_path
        self.config_path = config_path
        self.check_point_path = check_point_path
        self.results_dir = results_dir
        self.graph_data = torch.load(graph_data_path)
        self.config = yaml.safe_load(open(config_path))

        self.patient_feat_dim, self.gene_feat_dim = self.get_input_dimensions(self.graph_data)

    def get_input_dimensions(self, graph_data: Dict[str, Any]) -> Tuple[int, int]:
        """
        Extrae las dimensiones de características de pacientes y genes del primer grafo.

        Args:
            graph_data (Dict[str, Any]): Datos de grafos cargados

        Returns:
            Tuple[int, int]: Dimensiones de características (paciente, gene)
        """
        first_patient_key = list(graph_data.keys())[0]
        first_graph = graph_data[first_patient_key]['base']

        patient_feat_dim = first_graph['patient'].x.shape[1]
        gene_feat_dim = first_graph['gene'].x.shape[1]

        return patient_feat_dim, gene_feat_dim


    def run_inference(self, model, reference_score: float):

        predictions_list = []
        # Cargar los datos de prueba para el fold actual
        for patient in self.graph_data.keys():
            patient_graph = self.graph_data[patient]

            test_graph = list(patient_graph.values())
            test_graph_keys = list(patient_graph.keys())
            # Asegurarse de que el grafo tenga datos
            if not test_graph:
                continue
            # Convertir a Batch para inferencia
            if isinstance(test_graph, list):
                test_graph_list = Batch.from_data_list(test_graph)
            else:
                continue
            # Realizar inferencia
            predictions = model.predict(test_graph_list)

            ## Normalize
            predictions = (predictions - reference_score)/ np.log(2)
            if len(test_graph_keys) == 1:
                predictions = predictions.squeeze()
                # Si solo hay un grafo, usar el ID del paciente directamente
                pred_df = pd.DataFrame({
                    'ID': [test_graph_list['patient'].ID[0]], 
                    'Score': [predictions], 
                    'Gene Combination': test_graph_keys
                })
            else:
                # Si hay múltiples grafos, incluir la combinación de genes
                pred_df = pd.DataFrame(list(zip(test_graph_list['patient'].ID, predictions.squeeze(), test_graph_keys)), 
                    columns=['ID', 'Score', 'Gene Combination'])
            predictions_list.append(pred_df)
        return predictions_list

def parse_arguments():
    """
    Parsea los argumentos de línea de comandos para el script.
    
    Returns:
        argparse.Namespace: Argumentos parseados
    """
    parser = argparse.ArgumentParser(description='Run GNN inference with cross-validation')
    
    parser.add_argument('--graph_data_path', type=str, required=True,
                        help='Path to the graph data file (.pt)')
    parser.add_argument('--config_path', type=str, 
                        default="scripts/3.GeneticScores/configs/gnn_models_config.yaml",
                        help='Path to the configuration file')
    parser.add_argument('--model_name', type=str, required=True,
                        help='Name of the GNN model to use for inference')
    parser.add_argument('--reference_score', type=float, required=True,
                        help='Reference score for normalization')
    parser.add_argument('--check_point_path', type=str,
                        default="results/gnn/full_input_model/model_checkpoints",
                        help='Path to the checkpoint directory')
    parser.add_argument('--results_dir', type=str,
                        default="results/gnn/cv_results",
                        help='Directory to save results')
    
    return parser.parse_args()


def main():

    args = parse_arguments()

    gnn_inference = GNN_Inference_individual(
        graph_data_path=args.graph_data_path,
        config_path=args.config_path,
        check_point_path=args.check_point_path,
        results_dir=args.results_dir
    )
    
    model = instantiate_gnn_model(args.model_name, gnn_inference.config, 
                                gnn_inference.patient_feat_dim, gnn_inference.gene_feat_dim, 
                                args.check_point_path)

    predictions_list = gnn_inference.run_inference(model, reference_score=args.reference_score)
    # Combinar resultados de todos los folds
    out_df = pd.concat(predictions_list, ignore_index=True)
    # Guardar resultados en el directorio de resultados
    os.makedirs(args.results_dir, exist_ok=True)
    out_df.to_csv(f"{args.results_dir}/{args.model_name}_individual_gene_predictions.csv", index=False, sep = '\t')

if __name__ == "__main__":
    sys.exit(main())
