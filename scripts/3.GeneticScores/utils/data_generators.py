"""
Módulo para generar datos en diferentes formatos.
"""
import pandas as pd
import numpy as np
import logging
import pickle
import json
import torch
from pathlib import Path
from typing import Dict, List, Any, Optional, Tuple
from torch_geometric.data import HeteroData
from scipy.sparse import csr_matrix
from sklearn.compose import ColumnTransformer

logger = logging.getLogger(__name__)

class DataGenerators:
    """Clase para generar datos en formatos tabla y grafo."""
    
    def __init__(self, config: Dict[str, Any]):
        """Inicializa el generador de datos."""
        self.config = config
        self.data_generation = config['data_generation']
        self.fitted_transformer = None
        
        logger.info("DataGenerators inicializado")
    
    def generate_table_data(self, df_processed: pd.DataFrame, metadata: Dict[str, Any], 
                           output_path: Path) -> Dict[str, Any]:
        """Genera datos en formato tabla."""
        logger.info("Generando datos en formato tabla...")
        
        table_config = self.data_generation.get('table_format', {})
        output_filename = table_config.get('output_filename', 'table_data')
        
        # Extraer features
        feature_columns = metadata['feature_columns']
        X = df_processed[feature_columns].values
        
        # Extraer variables de supervivencia
        survival_vars = self.config['variable_processing']['survival_variables']
        time_var = survival_vars['time_variable']
        status_var = survival_vars['status_variable']
        
        y_time = df_processed[time_var].values
        y_event = df_processed[status_var].values
        
        # Extraer folds de CV
        cv_config = self.config.get('cross_validation', {})
        fold_column = cv_config.get('fold_column', 'fold')
        folds = df_processed[fold_column].values
        
        # Crear DataFrame final
        df_final = pd.DataFrame(X, columns=feature_columns)
        df_final[time_var] = y_time
        df_final[status_var] = y_event
        df_final[fold_column] = folds
        
        # Guardar archivo como pickle de pandas
        output_file = output_path / f"{output_filename}.pkl"
        
        output_config = self.config.get('output', {})
        
        # Usar pandas.to_pickle para guardar el DataFrame completo
        df_final.to_pickle(output_file)
        
        logger.info(f"Datos tabla guardados en formato pickle: {output_file}")
        logger.info(f"Shape final: {df_final.shape}")
        logger.info(f"Columnas: {list(df_final.columns)}")
        
        return {
            'file_path': str(output_file),
            'shape': df_final.shape,
            'feature_columns': feature_columns,
            'n_features': len(feature_columns),
            'dataframe_columns': list(df_final.columns)
        }
    
    def generate_graph_data(self, df_processed: pd.DataFrame, metadata: Dict[str, Any], 
                           df_mutations: Optional[pd.DataFrame], gene_embeddings: Optional[Dict[str, np.ndarray]],
                           output_path: Path) -> Dict[str, Any]:
        """Genera datos en formato grafo usando HeteroData."""
        logger.info("Generando datos en formato grafo...")
        
        graph_config = self.data_generation.get('graph_format', {})
        output_filename = graph_config.get('output_filename', 'graph_data')

        # Preparar datos de genes
        selected_genes = metadata['selected_genes']
        gene_features = self._prepare_gene_features(selected_genes, gene_embeddings)

        # Crear lista de grafos por paciente (como en prepareGraph.py)
        graph_list = []
        for patient_idx in range(len(df_processed)):
            graph = self._create_hetero_graph_for_patient(
                df_processed, patient_idx, gene_features, selected_genes, df_mutations
            )
            graph_list.append(graph)
        
        # Separar por folds (simulando train/test)
        cv_config = self.config.get('cross_validation', {})
        fold_column = cv_config.get('fold_column', 'fold')
        folds = df_processed[fold_column].values
        
        # Crear diccionario por fold
        graph_data = {}
        for fold in np.unique(folds):
            fold_indices = np.where(folds == fold)[0]
            graph_data[f'fold_{fold}'] = [graph_list[i] for i in fold_indices]
        
        # Guardar archivo
        output_file = output_path / f"{output_filename}.pt"
        torch.save(graph_data, output_file)
        
        logger.info(f"Datos grafo guardados en: {output_file}")
        logger.info(f"Total grafos: {len(graph_list)}")
        logger.info(f"Folds: {len(graph_data)}")
        
        return {
            'file_path': str(output_file),
            'n_graphs': len(graph_list),
            'n_folds': len(graph_data),
            'selected_genes': selected_genes
        }
    
    def _prepare_gene_features(self, selected_genes: List[str], 
                              gene_embeddings: Optional[Dict[str, np.ndarray]]) -> torch.Tensor:
        """Prepara features de genes como tensor."""
        if gene_embeddings is not None:
            # Usar embeddings
            embedding_dim = list(gene_embeddings.values())[0].shape[0]
            gene_features_matrix = np.zeros((len(selected_genes) + 1, embedding_dim))  # +1 para dummy gene
            
            for i, gene in enumerate(selected_genes):
                if gene in gene_embeddings:
                    gene_features_matrix[i] = gene_embeddings[gene]
                else:
                    # Embedding aleatorio para genes sin embedding
                    gene_features_matrix[i] = np.random.normal(0, 0.1, embedding_dim)
            
            # Gene dummy (como en prepareGraph.py)
            gene_features_matrix[-1] = np.zeros(embedding_dim)
            
            logger.info(f"Usando embeddings para genes: {embedding_dim}D")
        else:
            # Features one-hot
            n_genes = len(selected_genes) + 1
            gene_features_matrix = np.eye(n_genes)
            logger.info(f"Usando features one-hot para genes: {n_genes}D")
        
        return torch.from_numpy(gene_features_matrix).to(torch.float)
    
    def _create_hetero_graph_for_patient(self, df_processed: pd.DataFrame, patient_idx: int,
                                        gene_features: torch.Tensor, selected_genes: List[str],
                                        df_mutations: Optional[pd.DataFrame]) -> HeteroData:
        """Crea un HeteroData graph para un paciente específico."""
        # Extraer features del paciente (excluyendo genes)
        patient_row = df_processed.iloc[[patient_idx]]
        patient_x = self.fitted_transformer.transform(patient_row)
        patient_x = torch.tensor(patient_x, dtype=torch.float)

        # Extraer variables de supervivencia
        survival_vars = self.config['variable_processing']['survival_variables']
        time_var = survival_vars['time_variable']
        status_var = survival_vars['status_variable']
        
        patient_y = torch.tensor([[patient_row[time_var].squeeze(), patient_row[status_var].squeeze()]], dtype=torch.float)
        # Obtener genes mutados del paciente
        patient_genes_data = []
 
        if df_mutations is not None:
            # Usar datos VAF si están disponibles
            patient_id = patient_row.get('ID', f'Patient_{patient_idx:03d}').squeeze()
            patient_mutations = df_mutations[df_mutations['ID'] == patient_id]
            
            for gene in selected_genes:
                gene_mutation = patient_mutations[patient_mutations['GENE'] == gene]
                if not gene_mutation.empty:
                    vaf = gene_mutation.iloc[0]['VAF']
                    patient_genes_data.append({'GENE': gene, 'VAF': vaf})

        # Añadir gene dummy (como en prepareGraph.py)
        patient_genes_data.append({'GENE': 'Gene_0', 'VAF': 0.0})

        # Crear índices y pesos
        gene_indices = []
        vafs = []
        
        for gene_data in patient_genes_data:
            gene_name = gene_data['GENE']
            if gene_name == 'Gene_0':
                gene_idx = len(selected_genes)  # Último índice para dummy
            else:
                gene_idx = selected_genes.index(gene_name)
            
            gene_indices.append(gene_idx)
            vafs.append(gene_data['VAF'])
        # Seleccionar features de genes
        gene_sel_x = gene_features[gene_indices]
        
        # Crear edges (paciente -> genes)
        n_genes_patient = len(gene_indices)
        src = [0] * n_genes_patient  # Todos los edges vienen del paciente (nodo 0)
        dst = list(range(n_genes_patient))  # Hacia cada gen
        edge_index = torch.tensor([dst, src], dtype=torch.long)
        
        # Pesos de edges
        weights = torch.tensor(vafs, dtype=torch.float)
        
        # Crear HeteroData object
        graph = HeteroData()
        graph['patient'].x = patient_x
        graph['patient'].y = patient_y
        graph['patient'].ID = patient_id
        graph['gene'].x = gene_sel_x
        graph['gene'].weights = weights
        graph['gene', 'patient'].edge_index = edge_index
       
        return graph
    
    