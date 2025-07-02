import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import SAGEConv
import lightning as L
from pycox.models.loss import cox_ph_loss
from lifelines.utils import concordance_index

def safe_concordance_index(times, predictions, events, default_value=0.5):
    try:
        return concordance_index(times, predictions, events)
    except Exception as e:
        return default_value

def focal_cox_loss(predictions, times, events, alpha=0.25, gamma=2.0):
    """
    Cox loss con focal weighting para enfocarse en casos difíciles
    
    Args:
        predictions: Tensor de predicciones del modelo
        times: Tensor de tiempos de supervivencia
        events: Tensor de eventos (1=evento, 0=censurado)
        alpha: Factor de balance (típicamente 0.25)
        gamma: Factor de enfoque (típicamente 2.0)
        
    Returns:
        Focal Cox Loss que da más peso a casos difíciles
    """
    # Calcular Cox loss base
    base_loss = cox_ph_loss(predictions, times, events)
    
    # Calcular pesos focales
    # pt representa qué tan "fácil" es la predicción (mayor pt = más fácil)
    # Si base_loss es alta → pt es bajo → caso difícil
    # Si base_loss es baja → pt es alto → caso fácil
    pt = torch.exp(-base_loss)
    
    # Factor focal: (1 - pt)^gamma
    # Casos fáciles (pt alto): (1-pt) es pequeño → peso pequeño
    # Casos difíciles (pt bajo): (1-pt) es grande → peso grande
    focal_weight = alpha * (1 - pt) ** gamma
    
    # Loss final con peso focal
    focal_loss = focal_weight * base_loss
    
    return focal_loss

class PatientGNNSAGEFocal(L.LightningModule):
    """
    Patient GNN con Focal Cox Loss para enfocarse en casos difíciles de predecir.
    
    La Focal Loss ayuda cuando:
    - Hay casos atípicos difíciles de predecir
    - El modelo se estanca en casos "fáciles"
    - Queremos mejor performance en casos complejos
    """
    
    def __init__(self, patient_feat_dim, gene_feat_dim, hidden_gene_dim, hidden_dim, out_dim, 
                 use_vaf=False, hidden_vaf=0, learning_rate=0.001, 
                 focal_alpha=0.25, focal_gamma=2.0):
        super().__init__()
        
        # Guardar hiperparámetros automáticamente para logging y checkpointing
        self.save_hyperparameters()
        
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

        # VAF processing si se usa
        if use_vaf:
            if hidden_vaf > 0:
                # MLP for processing VAF features
                self.vaf_mal = nn.Sequential(
                    nn.Linear(1, hidden_vaf),
                    nn.ReLU(),
                    nn.Linear(hidden_vaf, gene_feat_dim)
                )
            else:
                self.vaf_mal = nn.Linear(1, gene_feat_dim)
        else:
            self.vaf_mal = None

        # GraphSAGE layers with sum aggregation
        self.sage1 = SAGEConv(in_channels=(hidden_dim, hidden_dim), 
                              out_channels=hidden_dim, 
                              aggr='sum')
        self.sage2 = SAGEConv(in_channels=(hidden_dim, hidden_dim), 
                              out_channels=hidden_dim, 
                              aggr='sum')

        # Fully connected layer for prediction
        self.fc = nn.Linear(hidden_dim, out_dim)
        
        # Hyperparameters - accesibles como self.hparams.learning_rate, etc.
        self.learning_rate = learning_rate
        self.use_vaf = use_vaf
        self.focal_alpha = focal_alpha
        self.focal_gamma = focal_gamma

    def forward(self, batch):
        """
        Forward pass del modelo
        
        Args:
            batch: Batch de datos heterogéneos con patient y gene features
            
        Returns:
            out: Predicciones de supervivencia [batch_size, 1]
        """
        # Get data from the batch
        patient_x = batch['patient'].x
        gene_x = batch['gene'].x
        
        # Apply MLP to patient features
        patient_x = self.patient_mlp(patient_x)  # Shape: [batch_size, hidden_dim]
        
        # Process gene features con VAF si está habilitado
        if self.use_vaf:
            weights = batch['gene'].weights
            vaf_transformed = self.vaf_mal(weights.unsqueeze(1))
            gene_x = gene_x + vaf_transformed
        
        gene_x = self.gene_mlp(gene_x)  # Shape: [num_genes, hidden_dim]
        
        # GraphSAGE layers: Update patient embeddings using gene information
        patient_x = F.relu(self.sage1((gene_x, patient_x), batch['gene', 'patient'].edge_index))
        patient_x = F.relu(self.sage2((gene_x, patient_x), batch['gene', 'patient'].edge_index))
        
        # Final prediction
        out = self.fc(patient_x)
        return out

    def training_step(self, batch, batch_idx):
        """
        Training step con Focal Cox Loss
        """
        # Get survival data
        patient_y = batch['patient'].y
        times = patient_y[:, 0]
        events = patient_y[:, 1]
        
        # Forward pass
        out = self.forward(batch)
        
        # Calcular Focal Cox Loss (pérdida principal)
        focal_loss = focal_cox_loss(out, times, events, self.focal_alpha, self.focal_gamma)
        
        # Calcular Cox Loss estándar para comparación
        cox_loss = cox_ph_loss(out, times, events)
        
        # Calcular métrica de rendimiento
        c_index = safe_concordance_index(times.cpu().numpy(), -out.cpu().detach().numpy(), events.cpu().numpy())
        
        # Logging
        batch_size = out.size(0)
        self.log("train_loss", focal_loss, prog_bar=True, on_epoch=True, batch_size=batch_size)
        self.log("train_focal_loss", focal_loss, prog_bar=False, on_epoch=True, batch_size=batch_size)
        self.log("train_cox_loss", cox_loss, prog_bar=False, on_epoch=True, batch_size=batch_size)
        self.log("train_c_index", c_index, prog_bar=True, on_epoch=True, batch_size=batch_size)
        
        # Logging adicional: diferencia entre focal y cox
        focal_cox_diff = focal_loss - cox_loss
        self.log("train_focal_cox_diff", focal_cox_diff, prog_bar=False, on_epoch=True, batch_size=batch_size)
        
        return focal_loss

    def validation_step(self, batch, batch_idx):
        """
        Validation step - usa Cox loss estándar para comparabilidad
        """
        # Get survival data
        patient_y = batch['patient'].y
        times = patient_y[:, 0]
        events = patient_y[:, 1]
        
        # Forward pass
        out = self.forward(batch)
        
        # En validación usar Cox loss estándar para comparar con otros modelos
        val_loss = cox_ph_loss(out, times, events)
        c_index = safe_concordance_index(times.cpu().numpy(), -out.cpu().detach().numpy(), events.cpu().numpy())
        
        # Logging
        batch_size = out.size(0)
        self.log("val_loss", val_loss, prog_bar=True, batch_size=batch_size)
        self.log("val_c_index", c_index, prog_bar=True, batch_size=batch_size)
        
        # También calcular focal loss en validación para monitoring
        val_focal_loss = focal_cox_loss(out, times, events, self.focal_alpha, self.focal_gamma)
        self.log("val_focal_loss", val_focal_loss, prog_bar=False, batch_size=batch_size)
        
        return val_loss

    def configure_optimizers(self):
        """
        Configurar optimizador - usa self.hparams.learning_rate si está disponible
        """
        # Usar self.hparams si está disponible, sino usar self.learning_rate
        lr = getattr(self.hparams, 'learning_rate', self.learning_rate)
        optimizer = torch.optim.Adam(self.parameters(), lr=lr)
        return optimizer