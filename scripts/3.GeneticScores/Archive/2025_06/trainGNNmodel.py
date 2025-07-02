"""
Train model
docker run --gpus all -v $PWD:$PWD -w $PWD -it mds_subtypes_python:1.4 python
"""
import lightning as L
import sys
import torch
from torch.utils.data import random_split
from torch_geometric.loader import DataLoader
from lightning.pytorch.loggers import TensorBoardLogger

sys.path.append('./scripts/3.GeneticScores/pytorch_files/')
from GNNPatientSAGE import PatientGNNSAGE
import json

patient_graphs = torch.load("results/gnn/preprocess/boolean_scGPT.pt")
patient_graphs_red = torch.load("results/gnn/preprocess/boolean_scGPT_reduced.pt")


# Define the split ratio
train_ratio = 0.8  # 80% for training
test_ratio = 0.2   # 20% for testing

# Compute the sizes for training and test sets
train_size = int(len(patient_graphs) * train_ratio)
test_size = len(patient_graphs) - train_size

# Split the dataset into training and testing subsets
train_dataset, test_dataset = random_split(patient_graphs, [train_size, test_size])
train_dataset_red, test_dataset_red = random_split(patient_graphs_red, [train_size, test_size])

# Hyperparameters
patient_feat_dim = 15
gene_feat_dim = 512
hidden_gene_dim = 32
hidden_dim = 16
out_dim = 1  # Single output for survival regression
learning_rate = 0.001
num_epochs = 20

gnn_model = PatientGNNSAGE(patient_feat_dim, gene_feat_dim, hidden_gene_dim, hidden_dim, out_dim, learning_rate)

# Create DataLoader instances for both training and test sets
train_loader = DataLoader(train_dataset, batch_size=32, shuffle=True)
test_loader = DataLoader(test_dataset, batch_size=128, shuffle=False)

# train model
logger_gnn = TensorBoardLogger("lightning_logs", name="gnn_all_features")
trainer_gnn = L.Trainer(max_epochs = num_epochs, logger = logger_gnn)
trainer_gnn.fit(model=gnn_model, train_dataloaders=train_loader, val_dataloaders=test_loader)

## Try to freeze patient MLP layer at the beggning of training
gnn_model2 = PatientGNNSAGE(patient_feat_dim, gene_feat_dim, hidden_gene_dim, hidden_dim, out_dim, learning_rate)
for param in gnn_model2.patient_mlp.parameters():
    param.requires_grad = False

logger_gnn2 = TensorBoardLogger("lightning_logs", name="gnn_freeze_20epochs")
trainer_gnn2 = L.Trainer(max_epochs = 50, logger = logger_gnn2)
trainer_gnn2.fit(model=gnn_model2, train_dataloaders=train_loader, val_dataloaders=test_loader)

for param in gnn_model2.patient_mlp.parameters():
    param.requires_grad = True

logger_gnn2b = TensorBoardLogger("lightning_logs", name="gnn_freeze_20epochs_unfreeze")
trainer_gnn2b = L.Trainer(max_epochs = num_epochs, logger = logger_gnn2b)
trainer_gnn2b.fit(model=gnn_model2, train_dataloaders=train_loader, val_dataloaders=test_loader)

## Evalute reduced model
train_red_loader = DataLoader(train_dataset_red, batch_size=32, shuffle=True)
test_red_loader = DataLoader(test_dataset_red, batch_size=128, shuffle=False)

gene_red_feat_dim = 16  # Reduced feature dimension for genes
gnn_model_red = PatientGNNSAGE(patient_feat_dim, gene_red_feat_dim, hidden_gene_dim, hidden_dim, out_dim, learning_rate)

# train model
logger_gnn_red = TensorBoardLogger("lightning_logs", name="gnn_all_features_reduced")
trainer_gnn_red = L.Trainer(max_epochs = num_epochs, logger = logger_gnn_red)
trainer_gnn_red.fit(model=gnn_model_red, train_dataloaders=train_red_loader, val_dataloaders=test_red_loader)

## This is the best model 

## Evalute reduced model with optimized hyperparameters
hidden_gene_dim = 12
hidden_dim = 54
learning_rate = 1e-4

gnn_model_red_opt = PatientGNNSAGE(patient_feat_dim, gene_red_feat_dim, hidden_gene_dim, hidden_dim, out_dim, learning_rate)
# train model
logger_gnn_red_opt = TensorBoardLogger("lightning_logs", name="gnn_all_features_reduced_opt")
trainer_gnn_red_opt = L.Trainer(max_epochs = num_epochs, logger = logger_gnn_red_opt)
trainer_gnn_red_opt.fit(model=gnn_model_red_opt, train_dataloaders=train_red_loader, val_dataloaders=test_red_loader)   
