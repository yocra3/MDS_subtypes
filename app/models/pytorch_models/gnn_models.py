"""
Modelos GNN para predicción de riesgo MDS usando PyTorch Lightning.

Este módulo contiene las implementaciones de los modelos GNN (Graph Neural Network)
basados en GraphSAGE para predicción de supervivencia en pacientes MDS.
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import SAGEConv
import lightning as L
from pycox.models.loss import cox_ph_loss
from lifelines.utils import concordance_index

class PatientGNNSAGE(L.LightningModule):
    """
    Modelo GNN para predicción de supervivencia en pacientes MDS.
    
    Utiliza GraphSAGE para procesar grafos heterogéneos paciente-genes
    y predice riesgo de supervivencia usando Cox Proportional Hazards.
    """
    
    def __init__(self, patient_feat_dim, gene_feat_dim, hidden_gene_dim, 
                 hidden_dim,  learning_rate=0.001):
        """
        Inicializa el modelo GNN.
        
        Args:
            patient_feat_dim: Dimensión de features de paciente
            gene_feat_dim: Dimensión de features de genes
            hidden_gene_dim: Dimensión oculta para procesamiento de genes
            hidden_dim: Dimensión oculta general
            use_vaf: Si usar VAF (Variant Allele Frequency)
            hidden_vaf: Dimensión oculta para procesamiento VAF
            learning_rate: Tasa de aprendizaje
        """
        super().__init__()
        
        # MLP for processing patient features
        self.patient_mlp = nn.Sequential(
            nn.Linear(patient_feat_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim),
        )
        
        # MLP for processing gene features
        self.gene_mlp = nn.Sequential(
            nn.Linear(gene_feat_dim, hidden_gene_dim),
            nn.ReLU(),
            nn.Linear(hidden_gene_dim, hidden_dim),
        )

        # GraphSAGE layers with sum aggregation
        self.sage1 = SAGEConv(in_channels=(hidden_dim, hidden_dim), 
                              out_channels=hidden_dim, 
                              aggr='sum')  # Sum aggregation
        self.sage2 = SAGEConv(in_channels=(hidden_dim, hidden_dim), 
                              out_channels=hidden_dim, 
                              aggr='sum')  # Sum aggregation

        # Fully connected layer for prediction
        self.fc = nn.Linear(hidden_dim, 1)
        self.learning_rate = learning_rate

    def forward(self, batch):
        """
        Forward pass del modelo.
        
        Args:
            batch: Batch de datos con estructura HeteroData
            
        Returns:
            torch.Tensor: Predicciones de riesgo
        """
        ## Get data from the batch
        patient_x = batch['patient'].x
        gene_x = batch['gene'].x
 
        # Apply MLP to patient features
        patient_x = self.patient_mlp(patient_x)  # Shape: [1, hidden_dim]
               
        gene_x = self.gene_mlp(gene_x)  # Shape: [num_genes, hidden_dim]
        
        # Apply GraphSAGE layers
        patient_x = F.relu(self.sage1((gene_x, patient_x), batch['gene', 'patient'].edge_index))
        patient_x = F.relu(self.sage2((gene_x, patient_x), batch['gene', 'patient'].edge_index))
        
        # Final prediction using updated patient embeddings
        out = self.fc(patient_x)
        return out

    def configure_optimizers(self):
        """Configura el optimizador."""
        optimizer = torch.optim.Adam(self.parameters(), lr=self.learning_rate)
        return optimizer


class BaseEvalGNN(L.LightningModule):
    """
    Clase base para evaluación de modelos GNN pre-entrenados.
    
    Carga un modelo pre-entrenado desde checkpoint y proporciona
    métodos para realizar predicciones.
    """
    
    def __init__(self, ckpt_path, model_class, **model_kwargs):
        """
        Inicializa el evaluador base.
        
        Args:
            ckpt_path: Ruta al checkpoint del modelo
            model_class: Clase del modelo a cargar
            **model_kwargs: Argumentos para inicializar el modelo
        """
        super().__init__()
        # init the pretrained LightningModule
        self.model = model_class(**model_kwargs)
        checkpoint = torch.load(ckpt_path, map_location='cpu')
        self.model.load_state_dict(checkpoint['state_dict'])
        self.eval()
    
    def forward(self, x):
        """Forward pass usando el modelo cargado."""
        return self.model(x)
    
    def predict(self, data):
        """
        Método unificado para hacer predicciones en grafos.
        
        Args:
            data: Datos de entrada (HeteroData)
            
        Returns:
            numpy.ndarray: Predicciones del modelo
        """
        with torch.no_grad():
            return self(data).detach().numpy()


class EvalGNN(BaseEvalGNN):
    """
    Evaluador específico para modelo PatientGNNSAGE.
    
    Wrapper para cargar y usar modelos PatientGNNSAGE pre-entrenados
    para realizar predicciones de riesgo.
    """
    
    def __init__(self, ck_path, patient_feat_dim, gene_feat_dim, hidden_gene_dim, 
                 hidden_dim, learning_rate=0.001):
        """
        Inicializa el evaluador GNN.
        
        Args:
            ck_path: Ruta al checkpoint del modelo
            patient_feat_dim: Dimensión de features de paciente
            gene_feat_dim: Dimensión de features de genes
            hidden_gene_dim: Dimensión oculta para genes
            hidden_dim: Dimensión oculta general
            out_dim: Dimensión de salida
            use_vaf: Si usar VAF weights
            hidden_vaf: Dimensión oculta VAF
            learning_rate: Tasa de aprendizaje
        """
        super().__init__(ck_path, PatientGNNSAGE,
                        patient_feat_dim=patient_feat_dim, gene_feat_dim=gene_feat_dim,
                        hidden_gene_dim=hidden_gene_dim, hidden_dim=hidden_dim,
                        learning_rate=learning_rate)
