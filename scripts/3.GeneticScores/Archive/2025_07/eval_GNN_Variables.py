"""
Module Summary:
We will evaluate the effect of the GNN variables on the model performance.
docker run --gpus all -v $PWD:$PWD -w $PWD -it mds_subtypes_python:1.5 python

"""

import os
import torch
import lightning as L
import sys
import pandas as pd
from scipy.sparse import csr_matrix
from torch_geometric.data import HeteroData

import torch.nn.functional as F
from torch import nn
from torch.utils.data import TensorDataset, DataLoader
import random
import torchtuples as tt
import numpy as np
from typing import Optional, Union, Tuple
from torch_geometric.data.batch import Batch

sys.path.append('./scripts/3.GeneticScores/pytorch_files/')
from GNNPatientSAGE import PatientGNNSAGE

## Load data
patients = pd.read_csv("results/mutations/patients_gnn_test.tsv", 
    sep = "\t")

gene_encoding_red = pd.read_csv("results/preprocess/gene_embedding_reduced.tsv", 
    sep = "\t", index_col=0)

scale_factors_df = pd.read_csv("results/mutations/scale_factors_gnn.tsv", 
    sep = "\t", index_col=0) 


## Select genes in common between embedding and mutation map
common_genes = set(patients.columns).intersection(gene_encoding_red.index)

gene_embed = gene_encoding_red.loc[list(common_genes), :].copy()
gene_embed.loc["Gene_0"] = 0  # Add dummy gene row with zeros
gene_x = torch.from_numpy(csr_matrix(gene_embed).todense()).to(torch.float)
gene_mapping = {index: i for i, index in enumerate(gene_embed.index)}


def define_graph(patient_vals, gene_embedding, gene_mapping):
    # Define patient features
    patient_feats = patient_vals.drop(['OS_YEARS', 'OS_STATUS', 'ID'] + list(common_genes), axis=1)
    patient_x = torch.from_numpy(csr_matrix(patient_feats).todense()).to(torch.float)
    ## Get genes in patient
    patient_genes_vars = patient_vals[list(common_genes)]
    patient_genes = patient_genes_vars.columns[patient_genes_vars.loc[0] == 1]
    patient_genes = pd.concat([pd.Series('Gene_0'), pd.Series(patient_genes)])
    indeces = list()
    for gene in patient_genes:
        indeces.append(gene_mapping[gene])
    gene_sel_x = gene_x[indeces, :]
    # Index
    src = [0] * gene_sel_x.shape[0]
    dst = list(range(gene_sel_x.shape[0]))
    edge_index = torch.tensor([dst, src])
   # Generate Object
    graph = HeteroData()
    graph['patient'].x = patient_x
    graph['gene'].x = gene_sel_x 
    graph['gene', 'patient'].edge_index = edge_index
    return graph

pat1 = define_graph(patients.loc[[0]], gene_embed, gene_mapping)
combined_graph = Batch.from_data_list(patient_graphs)

class EvalGNN(L.LightningModule):
    def __init__(self, ck_path, patient_feat_dim, gene_feat_dim, hidden_gene_dim, hidden_dim, out_dim, learning_rate=0.001):
        super().__init__()
        # init the pretrained LightningModule
        self.model =  PatientGNNSAGE(patient_feat_dim, gene_feat_dim, hidden_gene_dim, hidden_dim, out_dim, learning_rate)
        checkpoint = torch.load(ck_path, map_location='cpu')
        self.model.load_state_dict(checkpoint['state_dict'])       
    def forward(self, x):
        out = self.model(x)
        return out

patient_feat_dim = 15
gene_feat_dim = 16
hidden_gene_dim = 32
hidden_dim = 16
out_dim = 1  # Single output for survival regression
learning_rate = 0.001
num_epochs = 20

# Predict with GNN model
ck_gnn = "lightning_logs/gnn_all_features_reduced/version_0/checkpoints/epoch=19-step=1140.ckpt"

gnn_model = EvalGNN(ck_gnn, patient_feat_dim, gene_feat_dim, hidden_gene_dim, hidden_dim, out_dim, learning_rate)
gnn_model.eval()
gnn_model(pat1).detach().numpy()

def create_patient_graph_list(patient, gene_embed, gene_mapping):
    """
    Create a list of patient graphs for all variables different than 0
    """
    risk_vars = patient.columns[patient.loc[0] == 1] 
    risk_vars = pd.concat([pd.Series(["BM_BLAST",  "SEX"]), pd.Series(risk_vars)]).unique()
    return risk_vars

create_patient_graph_list(patients.loc[[0]], gene_embed, gene_mapping)
