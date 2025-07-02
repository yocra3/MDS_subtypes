import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import SAGEConv

class PatientGNNSAGE_noGene(nn.Module):
    def __init__(self, patient_feat_dim, hidden_dim, out_dim):
        super(PatientGNNSAGE_noGene, self).__init__()
        
        # MLP for processing patient features
        self.patient_mlp = nn.Sequential(
            nn.Linear(patient_feat_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim),
        )
        
        
        # Fully connected layer for prediction
        self.fc = nn.Linear(hidden_dim, out_dim)

    def forward(self, data):
        # Apply MLP to patient features
        patient_x = self.patient_mlp(data['patient'].x)  # Shape: [1, hidden_dim]
        # Apply MLP to gene features

        # Final prediction using updated patient embeddings
        out = self.fc(patient_x)
        return out
