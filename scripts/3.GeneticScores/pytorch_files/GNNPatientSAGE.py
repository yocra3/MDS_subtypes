import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import SAGEConv

class PatientGNNSAGE(nn.Module):
    def __init__(self, patient_feat_dim, gene_feat_dim, hidden_dim, out_dim):
        super(PatientGNNSAGE, self).__init__()
        
        # MLP for processing patient features
        self.patient_mlp = nn.Sequential(
            nn.Linear(patient_feat_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim),
        )
        
        # MLP for processing gene features
        self.gene_mlp = nn.Sequential(
            nn.Linear(gene_feat_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim),
        )
        
        # GraphSAGE layers with sum aggregation
        self.sage1 = SAGEConv(in_channels=(hidden_dim, hidden_dim), 
                              out_channels=hidden_dim, 
                              aggr='sum')  # Sum aggregation
        self.sage2 = SAGEConv(in_channels=(hidden_dim, hidden_dim), 
                              out_channels=hidden_dim, 
                              aggr='sum')  # Sum aggregation

        # Fully connected layer for prediction
        self.fc = nn.Linear(hidden_dim, out_dim)

    def forward(self, data):
        # Apply MLP to patient features
        patient_x = self.patient_mlp(data.x_patient)  # Shape: [1, hidden_dim]
        print(patient_x.shape)
        # Apply MLP to gene features
        gene_x = self.gene_mlp(data.x_gene)  # Shape: [num_genes, hidden_dim]

        # First GraphSAGE layer: Update patient embeddings using sum aggregation
        patient_x = F.relu(self.sage1((gene_x, patient_x), data.edge_index))
        print(patient_x.shape)
        # Second GraphSAGE layer: Further update patient embeddings using sum aggregation
        patient_x = F.relu(self.sage2((gene_x, patient_x), data.edge_index))
        print(patient_x.shape)
       # patient_x = patient_x.view(patient_x.shape[0], -1, hidden_dim).sum(dim=1)  # Shape: [num_patients, hidden_dim]

        # Final prediction using updated patient embeddings
        out = self.fc(patient_x)
        return out
