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

class PatientGNNSAGE_v2(L.LightningModule):
    def __init__(self, patient_feat_dim, gene_feat_dim, hidden_gene_dim, hidden_dim, out_dim, use_vaf = False, hidden_vaf = 0, learning_rate=0.001):
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
                              aggr='sum')  # Sum aggregation
        self.sage2 = SAGEConv(in_channels=(hidden_dim, hidden_dim), 
                              out_channels=hidden_dim, 
                              aggr='sum')  # Sum aggregation
        self.sage_genes = SAGEConv(in_channels=(hidden_dim, hidden_dim), 
                              out_channels=hidden_dim, 
                              aggr='sum')
        # Fully connected layer for prediction
        self.fc = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, out_dim)
        )
        self.learning_rate = learning_rate

        self.use_vaf = use_vaf 
    def training_step(self, batch, batch_idx):
        ## Get data from the batch
        patient_x = batch['patient'].x
        patient_y = batch['patient'].y
        times = patient_y[:, 0]
        events = patient_y[:, 1]
        gene_x = batch['gene'].x          
        # Apply MLP to patient features
        patient_x = self.patient_mlp(patient_x)  # Shape: [1, hidden_dim]
        # Apply MLP to gene features
        if (self.use_vaf):
            weights = batch['gene'].weights
            vaf_transformed = self.vaf_mal(weights.unsqueeze(1))  # Transform VAF to hidden_gene_dim
            gene_x = gene_x + vaf_transformed
        gene_x = self.gene_mlp(gene_x)  # Shape: [num_genes, hidden_dim]
        # First GraphSAGE layer: Update patient embeddings using sum aggregation
        patient_x1 = F.relu(self.sage1((gene_x, patient_x), batch['gene', 'patient'].edge_index))
        # Update gene embeddings using sum aggregation
        gene_x = F.relu(self.sage_genes((patient_x1, gene_x), batch['gene', 'patient'].edge_index[[1, 0], :]))
        # Second GraphSAGE layer: Update patient embeddings using sum aggregation
        patient_x = F.relu(self.sage2((gene_x, patient_x1), batch['gene', 'patient'].edge_index))
        # Final prediction using updated patient embeddings
        out = self.fc(patient_x)
        loss = cox_ph_loss(out, times, events)
        c_index = safe_concordance_index(times.cpu().numpy(), -out.cpu().detach().numpy(), events.cpu().numpy())
        batch_size = patient_x.size(0)
        self.log("train_loss", loss, prog_bar=True, on_epoch=True, batch_size=batch_size)
        self.log("train_c_index", c_index, prog_bar=True, on_epoch=True, batch_size=batch_size)
        return loss

    def validation_step(self, batch, batch_idx):
        ## Get data from the batch
        patient_x = batch['patient'].x
        patient_y = batch['patient'].y
        times = patient_y[:, 0]
        events = patient_y[:, 1]
        gene_x = batch['gene'].x
        
        # Apply MLP to patient features
        patient_x = self.patient_mlp(patient_x)  # Shape: [1, hidden_dim]

        if (self.use_vaf):
            weights = batch['gene'].weights
            vaf_transformed = self.vaf_mal(weights.unsqueeze(1))  # Transform VAF to hidden_gene_dim
            gene_x = gene_x + vaf_transformed
        gene_x = self.gene_mlp(gene_x)  # Shape: [num_genes, hidden_dim]
       # First GraphSAGE layer: Update patient embeddings using sum aggregation
        patient_x1 = F.relu(self.sage1((gene_x, patient_x), batch['gene', 'patient'].edge_index))
        # Update gene embeddings using sum aggregation
        gene_x = F.relu(self.sage_genes((patient_x1, gene_x), batch['gene', 'patient'].edge_index[[1, 0], :]))
        # Second GraphSAGE layer: Update patient embeddings using sum aggregation
        patient_x = F.relu(self.sage2((gene_x, patient_x1), batch['gene', 'patient'].edge_index))
     # Final prediction using updated patient embeddings
        out = self.fc(patient_x)
        val_loss = cox_ph_loss(out, times, events)
        c_index = safe_concordance_index(times.cpu().numpy(), -out.cpu().detach().numpy(), events.cpu().numpy())
        
        # Especificar explícitamente batch_size
        batch_size = patient_x.size(0)
        self.log("val_loss", val_loss, prog_bar=True, batch_size=batch_size)
        self.log("val_c_index", c_index, prog_bar=True, batch_size=batch_size)
        return val_loss

    def configure_optimizers(self):
        optimizer = torch.optim.Adam(self.parameters(), lr=self.learning_rate)
        return optimizer

    def forward(self, batch):
        ## Get data from the batch
        patient_x = batch['patient'].x
        gene_x = batch['gene'].x
        # Apply MLP to patient features
        patient_x = self.patient_mlp(patient_x)  # Shape: [1, hidden_dim]
        if (self.use_vaf):
            weights = batch['gene'].weights
            vaf_transformed = self.vaf_mal(weights.unsqueeze(1))  # Transform VAF to hidden_gene_dim
            gene_x = gene_x + vaf_transformed
        gene_x = self.gene_mlp(gene_x)  # Shape: [num_genes, hidden_dim]
        patient_x1 = F.relu(self.sage1((gene_x, patient_x), batch['gene', 'patient'].edge_index))
        # Update gene embeddings using sum aggregation
        gene_x = F.relu(self.sage_genes((patient_x1, gene_x), batch['gene', 'patient'].edge_index[[1, 0], :]))
        # Second GraphSAGE layer: Update patient embeddings using sum aggregation
        patient_x = F.relu(self.sage2((gene_x, patient_x1), batch['gene', 'patient'].edge_index))
        # Final prediction using updated patient embeddings
        out = self.fc(patient_x)
        return out