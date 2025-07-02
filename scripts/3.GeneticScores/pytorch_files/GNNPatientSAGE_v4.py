import torch
import torch.nn as nn
import torch.nn.functional as F
import lightning as L
from torch_geometric.nn import SAGEConv
from pycox.models.loss import cox_ph_loss
from lifelines.utils import concordance_index

def safe_concordance_index(times, predictions, events, default_value=0.5):
    try:
        return concordance_index(
            times.cpu().numpy(), 
            -predictions.cpu().detach().numpy(), 
            events.cpu().numpy()
        )
    except Exception as e:
        print(f"Error al calcular concordance index: {str(e)}")
        return default_value

class PatientGNNSAGE_v4(L.LightningModule):
    def __init__(self, patient_feat_dim, gene_feat_dim, hidden_gene_dim, hidden_dim, out_dim, 
                 dropout=0.3, l2_reg=1e-4, use_vaf=False, hidden_vaf=0, learning_rate=0.001):
        super().__init__()
        self.save_hyperparameters()
        
        self.use_vaf = use_vaf
        self.learning_rate = learning_rate
        self.l2_reg = l2_reg
        self.dropout = dropout
        
        # MLP para patient features con Dropout y BatchNorm
        self.patient_mlp = nn.Sequential(
            nn.Linear(patient_feat_dim, hidden_dim),
            nn.BatchNorm1d(hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, hidden_dim),
            nn.BatchNorm1d(hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout)
        )
        
        # MLP para gene features con Dropout y BatchNorm
        self.gene_mlp = nn.Sequential(
            nn.Linear(gene_feat_dim, hidden_gene_dim),
            nn.BatchNorm1d(hidden_gene_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_gene_dim, hidden_dim),
            nn.BatchNorm1d(hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout)
        )
        
        # VAF processing si se usa
        if self.use_vaf:
            self.vaf_mlp = nn.Sequential(
                nn.Linear(1, hidden_vaf),
                nn.ReLU(),
                nn.Dropout(dropout),
                nn.Linear(hidden_vaf, hidden_dim)
            )
        
        # GraphSAGE layers - MANTENER agregación 'sum' como en el original
        self.sage1 = SAGEConv((hidden_dim, hidden_dim), hidden_dim, aggr='sum')
        self.sage2 = SAGEConv((hidden_dim, hidden_dim), hidden_dim, aggr='sum')
        
        # Final prediction layer con Dropout
        self.fc = nn.Sequential(
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, out_dim)
        )

    def forward(self, batch):
        # Get data from the batch
        patient_x = batch['patient'].x
        gene_x = batch['gene'].x
        
        # Apply MLP to patient features
        patient_x = self.patient_mlp(patient_x)
        
        # Process gene features
        if self.use_vaf:
            weights = batch['gene'].weights.unsqueeze(1)
            vaf_transformed = self.vaf_mlp(weights)
            gene_x = self.gene_mlp(gene_x)
            gene_x = gene_x + vaf_transformed
        else:
            gene_x = self.gene_mlp(gene_x)
        
        # GraphSAGE layers con Dropout aplicado después
        patient_x = F.relu(self.sage1((gene_x, patient_x), batch['gene', 'patient'].edge_index))
        patient_x = F.dropout(patient_x, p=self.dropout, training=self.training)
        
        patient_x = F.relu(self.sage2((gene_x, patient_x), batch['gene', 'patient'].edge_index))
        patient_x = F.dropout(patient_x, p=self.dropout, training=self.training)
        
        # Final prediction
        out = self.fc(patient_x)
        
        return out

    def training_step(self, batch, batch_idx):
        # Get survival data
        patient_y = batch['patient'].y
        times = patient_y[:, 0]
        events = patient_y[:, 1]
        
        # Forward pass
        out = self.forward(batch)
        
        # Cox loss principal
        cox_loss = cox_ph_loss(out, times, events)
        
        # L2 Regularization - calcular la norma L2 de todos los parámetros
        l2_loss = 0
        for param in self.parameters():
            l2_loss += torch.norm(param, 2)
        
        # Loss total = Cox loss + L2 regularization
        total_loss = cox_loss + self.l2_reg * l2_loss
        
        # Calculate metrics
        c_index = safe_concordance_index(times, out, events)
        
        # Log metrics
        batch_size = patient_y.size(0)
        self.log("train_loss", total_loss, prog_bar=True, on_epoch=True, batch_size=batch_size)
        self.log("train_cox_loss", cox_loss, prog_bar=False, on_epoch=True, batch_size=batch_size)
        self.log("train_l2_loss", l2_loss, prog_bar=False, on_epoch=True, batch_size=batch_size)
        self.log("train_c_index", c_index, prog_bar=True, on_epoch=True, batch_size=batch_size)
        
        return total_loss

    def validation_step(self, batch, batch_idx):
        # Get survival data
        patient_y = batch['patient'].y
        times = patient_y[:, 0]
        events = patient_y[:, 1]
        
        # Forward pass
        out = self.forward(batch)
        
        # Calculate loss (solo Cox loss en validación)
        val_loss = cox_ph_loss(out, times, events)
        c_index = safe_concordance_index(times, out, events)
        
        # Log metrics
        batch_size = patient_y.size(0)
        self.log("val_loss", val_loss, prog_bar=True, batch_size=batch_size)
        self.log("val_c_index", c_index, prog_bar=True, batch_size=batch_size)
        
        return val_loss

    def configure_optimizers(self):
        # Usar Adam estándar (mantenemos el optimizador original)
        optimizer = torch.optim.Adam(self.parameters(), lr=self.learning_rate)
        return optimizer