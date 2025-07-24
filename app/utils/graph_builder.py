"""
Constructor de grafos heterogéneos para modelos GNN de predicción de riesgo MDS.

Este módulo se encarga de crear estructuras de grafos heterogéneos que representan
la relación entre un paciente y sus genes mutados, para ser usados por modelos GNN.
"""

import torch
import numpy as np
import pandas as pd
import logging
from typing import Dict, List, Any, Optional, Tuple
from torch_geometric.data import HeteroData
from sklearn.compose import ColumnTransformer

from models.base_model import PatientData

logger = logging.getLogger(__name__)


class GraphBuilder:
    """
    Constructor de grafos heterogéneos paciente-genes.
    
    Crea estructuras de grafo donde:
    - Nodos paciente: contienen features clínicas
    - Nodos genes: contienen embeddings de genes
    - Edges: conectan paciente con genes mutados
    """
    
    def __init__(self, config: Optional[Dict[str, Any]] = None, 
                 transformer: Optional[ColumnTransformer] = None):
        """
        Inicializa el constructor de grafos.
        
        Args:
            config: Configuración opcional para el constructor
            transformer: ColumnTransformer pre-entrenado para features clínicas
        """
        self.config = config or {}
        self.transformer = transformer
        logger.info("GraphBuilder inicializado")
    
    def build_patient_graph(self, clinical_features: torch.Tensor, 
                           mutated_genes: List[str],
                           gene_embeddings: Dict[str, np.ndarray]) -> HeteroData:
        """
        Construye un grafo heterogéneo para un paciente.
        
        Args:
            clinical_features: Tensor con features clínicas del paciente [1, n_features]
            mutated_genes: List[str] - Lista de genes mutados
            gene_embeddings: Dict[str, np.ndarray] - Diccionario con embeddings por gen
        Returns:
            HeteroData: Grafo heterogéneo paciente-genes
        """
        try:
            # 1. Preparar datos del paciente
            patient_x = clinical_features.clone()
            mutated_genes.append("Gene_0")
            # 2. Obtener embeddings de genes mutados
            gene_sel_x = self._get_gene_embeddings(mutated_genes, gene_embeddings)
            
            # 3. Crear estructura de edges
            edge_index = self._create_edges(len(mutated_genes))
            
            # 4. Construir HeteroData object
            graph = self._create_hetero_data(
                patient_x, gene_sel_x, edge_index, mutated_genes
            )
            
            logger.debug(f"Grafo creado: {len(mutated_genes)} genes mutados")
            return graph
            
        except Exception as e:
            logger.error(f"Error construyendo grafo: {e}")
            raise
    
    def build_patient_graph_from_data(self, patient_data: PatientData,
                                     gene_embeddings: Dict[str, np.ndarray]) -> HeteroData:
        """
        Construye grafo directamente desde PatientData.
        
        Args:
            patient_data: Datos estructurados del paciente
            gene_embeddings: Dict[str, np.ndarray] - Diccionario con embeddings por gen
            
        Returns:
            HeteroData: Grafo heterogéneo
        """
        # Convertir PatientData a tensores
        clinical_features = self._patient_data_to_tensor(patient_data)
        mutated_genes = patient_data.mutations
        # Usar mutaciones del paciente directamente
        return self.build_patient_graph(clinical_features, mutated_genes, 
                                      gene_embeddings)
    
    def _get_gene_embeddings(self, mutated_genes: List[str], 
                            gene_embeddings: Dict[str, np.ndarray]) -> torch.Tensor:
        """
        Obtiene embeddings de genes mutados desde diccionario.
        
        Args:
            mutated_genes: Lista de genes mutados
            gene_embeddings: Diccionario con embeddings por gen
            
        Returns:
            torch.Tensor: Embeddings concatenados de genes mutados
        """
        embedding_dim = 16 ## Hard-coded

        gene_embedding_list = []
        for gene in mutated_genes:
            if gene in gene_embeddings:
                embedding = torch.from_numpy(gene_embeddings[gene]).float()
                gene_embedding_list.append(embedding)
            else:
                logger.warning(f"Gen {gene} no encontrado en embeddings")

        # Concatenar embeddings
        gene_sel_x = torch.stack(gene_embedding_list)
        return gene_sel_x
    
    def _create_edges(self, num_genes: int) -> torch.Tensor:
        """
        Crea conexiones unidireccionales entre paciente y genes mutados.
        
        Args:
            num_genes: Número de genes mutados
            
        Returns:
            edge_index tensor [2, num_edges]
        """
        if num_genes == 0:
            return torch.empty((2, 0), dtype=torch.long)

        # Conexiones unidireccionales: paciente -> genes
        target_nodes = [0] * num_genes
        source_nodes = list(range(num_genes))
        
        edge_index = torch.tensor([source_nodes, target_nodes], dtype=torch.long)
        
        return edge_index
    
    def _create_hetero_data(self, patient_x: torch.Tensor, gene_x: torch.Tensor,
                           edge_index: torch.Tensor,
                           mutated_genes: List[str]) -> HeteroData:
        """
        Crea el objeto HeteroData final.
        
        Args:
            patient_x: Features del paciente
            gene_x: Features de genes
            edge_index: Indices de edges
            mutated_genes: Lista de genes mutados
            
        Returns:
            HeteroData: Grafo heterogéneo completo
        """
        graph = HeteroData()
        
        # Nodos paciente
        graph['patient'].x = patient_x
        
        # Nodos genes
        graph['gene'].x = gene_x
        graph['gene'].gene_names = mutated_genes
        
        # Edges paciente -> genes
        graph['gene', 'patient'].edge_index = edge_index
        
        return graph
    
    def _patient_data_to_tensor(self, patient_data: PatientData) -> torch.Tensor:
        """
        Convierte PatientData a tensor de features clínicas usando transformer pre-entrenado.
        
        Args:
            patient_data: Datos del paciente
            
        Returns:
            torch.Tensor: Features clínicas transformadas [1, n_features]
        """
        # Crear DataFrame temporal para el transformer
        patient_df = pd.DataFrame([{
            'AGE': patient_data.AGE,
            'BM_BLAST': patient_data.BM_BLAST,
            'WBC': patient_data.WBC,
            'ANC': patient_data.ANC,
            'HB': patient_data.HB,
            'PLT': patient_data.PLT,
            'MONOCYTES': patient_data.MONOCYTES,
            'SEX': patient_data.SEX,
            'plus8': patient_data.plus8,
            'del7': patient_data.del7,
            'del20q': patient_data.del20q,
            'del7q': patient_data.del7q,
            'delY': patient_data.delY,
            'del5q': patient_data.del5q,
            'complex': patient_data.complex,
            'CYTO_IPSSR': patient_data.CYTO_IPSSR
        }])
        
        # Aplicar transformer pre-entrenado
        features_transformed = self.transformer.transform(patient_df)
        
        return torch.tensor(features_transformed, dtype=torch.float32)
    
    def validate_graph(self, graph: HeteroData) -> bool:
        """
        Valida que el grafo esté bien formado.
        
        Args:
            graph: Grafo a validar
            
        Returns:
            bool: True si el grafo es válido
        """
        try:
            # Verificar que existan los nodos requeridos
            assert 'patient' in graph.node_types
            assert 'gene' in graph.node_types
            
            # Verificar dimensiones
            assert graph['patient'].x.shape[0] == 1  # Un paciente
            assert graph['gene'].x.shape[0] > 0  # Al menos un gen
           
            # Verificar edges
            assert ('gene', 'to', 'patient') in graph.edge_types
            edge_index = graph['gene', 'patient'].edge_index
            
            assert edge_index.shape[0] == 2  # [src, dst]
            assert edge_index.shape[1] == graph['gene'].x.shape[0]  # Un edge por gen
            
            logger.debug("Grafo validado correctamente")
            return True
            
        except Exception as e:
            logger.error(f"Grafo inválido: {e}")
            return False


# Funciones de utilidad

def create_dummy_gene_embeddings(gene_list: List[str], embedding_dim: int = 16) -> Dict[str, np.ndarray]:
    """
    Crea embeddings dummy para testing.
    
    Args:
        gene_list: Lista de genes
        embedding_dim: Dimensión de embeddings
        
    Returns:
        Dict[str, np.ndarray]: Diccionario con embeddings dummy por gen
    """
    embeddings_dict = {}
    for gene in gene_list:
        embeddings_dict[gene] = np.random.randn(embedding_dim)
    
    # Agregar gene dummy
    embeddings_dict['Gene_0'] = np.zeros(embedding_dim)
    
    return embeddings_dict


def load_gene_embeddings_from_file(file_path: str) -> Dict[str, np.ndarray]:
    """
    Carga embeddings de genes desde archivo.
    
    Args:
        file_path: Ruta al archivo de embeddings
        
    Returns:
        Dict[str, np.ndarray]: Diccionario con embeddings cargados
    """
    try:

        embeddings_df = pd.read_csv(file_path, sep="\t", index_col=0)
            
        gene_embeddings = {}
        for gene_name in embeddings_df.index:
            gene_embeddings[gene_name] = embeddings_df.loc[gene_name].values.astype(np.float32)
            
        # Añadir embedding dummy para "Gene_0" (gene ficticio)
        embedding_dim = embeddings_df.shape[1]
        gene_embeddings["Gene_0"] = np.zeros(embedding_dim, dtype=np.float32)        
               
        logger.info(f"Embeddings cargados desde {file_path}: {len(gene_embeddings)} genes")
        return gene_embeddings
    except Exception as e:
        logger.error(f"Error cargando embeddings: {e}")
        raise
