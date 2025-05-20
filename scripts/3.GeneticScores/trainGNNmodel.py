"""
Train model
docker run --gpus all -v $PWD:$PWD -w $PWD -it mds_subtypes_python:1.2 python
"""

import torch
import torch.nn.functional as F
from torch_geometric.nn import GraphConv
from torch_geometric.data import Data
from torch_geometric.loader import DataLoader
from torch.utils.data import random_split

from torchinfo import summary
import sys
from pycox.models import CoxPH
import random
import torchtuples as tt

sys.path.append('./scripts/3.GeneticScores/pytorch_files/')
from GNNPatientSAGE import PatientGNNSAGE


# Create a single graph for a patient
def create_patient_graph(patient_feat, gene_feats, edge_index, label):
    """
    Create a graph representing a single patient with connections to genes.
    Args:
        patient_feat: Tensor of shape [1, patient_feat_dim]
        gene_feats: Tensor of shape [num_genes, gene_feat_dim]
        edge_index: Tensor of shape [2, num_edges], edges between patient and genes
        label: Scalar label for the patient (e.g., survival in months)
    Returns:
        Data object for the patient graph.
    """
    data = Data()
    data.x_patient = patient_feat  # Feature for the patient node
    data.x_gene = gene_feats       # Features for gene nodes
    data.edge_index = edge_index   # Edges connecting the patient to genes
    data.y = torch.tensor([label], dtype=torch.float)  # Patient label
    return data

# Create a list of graphs for all patients
def create_patient_graphs():
    num_patients = 5
    num_genes = 6
    patient_feat_dim = 10
    gene_feat_dim = 8
    graphs = []
    for i in range(num_patients):
        # Generate random features and labels for demonstration
        patient_feat = torch.rand((1, patient_feat_dim))  # Single patient node
        num_genes = random.randint(3, 6)
        gene_feats = torch.rand((num_genes, gene_feat_dim))  # Gene nodes
        edge_index = torch.tensor([
            list(range(num_genes)),  
            [0] * num_genes, 
        ], dtype=torch.long)
        label = [torch.randint(10, 50, (1,)), torch.randint(0, 1, (1, ))]  # Random survival label
        # Create a graph for the patient
        graphs.append(create_patient_graph(patient_feat, gene_feats, edge_index, label))   
    return graphs

# Create a list of patient graphs
patient_graphs = create_patient_graphs()

# Define the split ratio
train_ratio = 0.8  # 80% for training
test_ratio = 0.2   # 20% for testing

# Compute the sizes for training and test sets
train_size = int(len(patient_graphs) * train_ratio)
test_size = len(patient_graphs) - train_size

# Split the dataset into training and testing subsets
train_dataset, test_dataset = random_split(patient_graphs, [train_size, test_size])

# Create DataLoader instances for both training and test sets
train_loader = DataLoader(train_dataset, batch_size=32, shuffle=True)
test_loader = DataLoader(test_dataset, batch_size=32, shuffle=False)

# Hyperparameters
patient_feat_dim = 10
gene_feat_dim = 8
hidden_dim = 16
out_dim = 1  # Single output for survival regression
learning_rate = 0.01
num_epochs = 100
batch_size = 2

model = PatientGNNSAGE(patient_feat_dim, gene_feat_dim, hidden_dim, out_dim)
summary(model = model)

model(patient_graphs[0])    
for batch in train_loader:
    batch.edge_index[1] = batch.batch
    batch

model(batch)


def coxph_loss(preds, times, events):
    """
    Computes the Cox Proportional Hazards loss.
    Args:
        preds (torch.Tensor): Predicted risk scores (higher means higher risk).
        times (torch.Tensor): Observed survival times.
        events (torch.Tensor): Event indicators (1 if event occurred, 0 if censored).

    Returns:
        torch.Tensor: Negative partial log-likelihood loss.
    """
    # Sort by descending time
    idx = torch.argsort(times, descending=True)
    preds, events = preds[idx], events[idx]  
    # Compute the log partial likelihood
    log_cumulative_hazard = torch.log(torch.cumsum(torch.exp(preds), dim=0))
    log_likelihood = preds - log_cumulative_hazard
    # Mask with event indicator
    loss = -torch.sum(log_likelihood * events)
    return loss

# Initialize model, optimizer, and loss
optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)

# Training loop
for epoch in range(num_epochs):
    model.train()
    total_loss = 0
    for batch in train_loader:
        optimizer.zero_grad()
        batch.edge_index[1] = batch.batch
        # Forward pass
        risk_scores = model(batch)
        # Compute loss
        times = batch.y[:,0]  # Observed survival times
        events = batch.y[:,1]  # Event indicators
        print(risk_scores, times, events)
        # Compute CoxPH loss
        loss = coxph_loss(risk_scores, times, events)
        loss.backward()
        optimizer.step()
        total_loss += loss.item()
    # Print average loss per epoch
    if epoch % 10 == 0:
        print(f"Epoch {epoch}, Loss: {total_loss / len(patient_graphs):.4f}")

    