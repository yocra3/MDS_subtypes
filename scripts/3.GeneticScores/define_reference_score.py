# Define a baseline reference score for GNN model of fold 0. To do so,
# compute the score for all patients and take the median as the reference.

import torch
import yaml
import sys
import argparse
from torch_geometric.data.batch import Batch
import pandas as pd
import glob
import os

# Agregar ruta para importar clases de evaluación GNN
sys.path.append('./scripts/3.GeneticScores/utils/')
from evalGNN_classes import EvalGNN
from gnn_inference_utils import get_graph_dimensions, find_checkpoint, instantiate_gnn_model


class GNN_Inference:
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
        Inicializa la clase GNN_Inference.
        
        Args:
            graph_data_path (str): Ruta al archivo con los datos de grafos (.pt)
            config_path (str): Ruta al archivo de configuración YAML
            check_point_path (str): Ruta al directorio con los checkpoints de modelos
            results_dir (str): Directorio donde guardar los resultados
        """

        self.graph_data_path = graph_data_path
        self.config_path = config_path
        self.check_point_path = check_point_path
        self.graph_data = torch.load(graph_data_path)
        self.config = yaml.safe_load(open(config_path))

        self.patient_feat_dim, self.gene_feat_dim = get_graph_dimensions(self.graph_data)


    def instantiate_models(self, model_name: str, fold: int = 0):
        """
        Instancia un modelo GNN específico para un fold de validación cruzada.
        
        Args:
            model_name (str): Nombre del modelo a instanciar
            fold (int): Número del fold de validación cruzada
            
        Returns:
            model: Instancia del modelo GNN cargado con el checkpoint correspondiente
            
        Raises:
            FileNotFoundError: Si no se encuentran checkpoints para el modelo especificado
            ValueError: Si el tipo de modelo no es reconocido
        """
        print(f"Instantiating model {model_name} for fold {fold}")
        
        checkpoint_path = find_checkpoint(self.check_point_path, model_name, fold)
        model = instantiate_gnn_model(model_name, self.config, 
                                    self.patient_feat_dim, self.gene_feat_dim, 
                                    checkpoint_path)
        return model

    def run_GNN_inference(self, model_name: str, fold: int = 0):
        """
        Ejecuta inferencia para un fold específico de validación cruzada.
        
        Args:
            model_name (str): Nombre del modelo a usar
            fold (int): Número del fold a procesar
            
        Returns:
            pandas.DataFrame: DataFrame con IDs de pacientes y predicciones
        """

        model = self.instantiate_models(model_name, fold)
        
        predictions_list = []
        # Cargar los datos de prueba para el fold actual
        for key in self.graph_data.keys():
            test_graph = self.graph_data[key]
            test_graph_list = Batch.from_data_list(test_graph)
            # Realizar inferencia
            predictions = model.predict(test_graph_list)
            pred_df = pd.DataFrame(list(zip(test_graph_list['patient'].ID, predictions.squeeze())), 
                  columns=['ID', model_name])
            predictions_list.append(pred_df)
        # Combinar resultados de todos los folds
        pred_df = pd.concat(predictions_list, ignore_index=True)
        # Guardar resultados en el directorio de resultados
        return pred_df



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
    parser.add_argument('--fold', required=True,
                        help='List of fold to define the model (e.g., 0 1 2)')
    parser.add_argument('--check_point_path', type=str,
                        default="results/gnn/full_input_model/model_checkpoints",
                        help='Path to the checkpoint directory')
    parser.add_argument('--results_dir', type=str,
                        default="results/gnn/cv_results",
                        help='Directory to save results')
    
    return parser.parse_args()

def main():
    """
    Función principal del script.
    
    Ejecuta la inferencia de modelos GNN para validación cruzada y guarda los resultados.
    """
    args = parse_arguments()
    
    # Inicializar la clase de inferencia
    gnn_inference = GNN_Inference(
        graph_data_path=args.graph_data_path,
        config_path=args.config_path,
        check_point_path=args.check_point_path,
        results_dir=args.results_dir
    )
    
    results_df = gnn_inference.run_GNN_inference(args.model_name, fold = args.fold)
    print(results_df[args.model_name].describe())
    print(f"Median: {results_df[args.model_name].median()}")
    
    print("Inference completed for all models.")

if __name__ == "__main__":
    sys.exit(main())

# count    2336.000000
# mean        1.032704
# std         0.947058
# min        -0.430854
# 25%         0.266818
# 50%         0.838476
# 75%         1.694185
# max         4.319068
# Name: v1_boolean, dtype: float64
# Median: 0.8384760022163391