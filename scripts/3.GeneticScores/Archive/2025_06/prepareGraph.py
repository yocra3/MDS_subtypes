"""
This is the data preprocessing that should return 
the .pt with the data to run the model 

docker run --gpus all -v $PWD:$PWD -w $PWD -it mds_subtypes_python:1.4 python

"""

import os
import pandas as pd

import torch
from torch_geometric.data import HeteroData
from scipy.sparse import csr_matrix

## Read data
patient_features = pd.read_csv("results/mutations/patients_features_gnn.tsv", 
    sep = "\t")

mutation_map = pd.read_csv("results/mutations/mutation_map_gnn.tsv", 
    sep = "\t")

mutation_maf = pd.read_csv("results/mutations/mutation_maf_gnn.tsv", 
    sep = "\t")

gene_encoding = pd.read_csv("results/preprocess/gene_embedding_scgpt.tsv", 
    sep = "\t", index_col=0)

gene_encoding_red = pd.read_csv("results/preprocess/gene_embedding_reduced.tsv", 
    sep = "\t", index_col=0)


## Remove rows with NaaN
mutation_maf = mutation_maf.dropna(subset=["VAF"])

### Generate heterodata object
# Features
patient_features_sel = patient_features.drop(['OS_YEARS', 'OS_STATUS', 'ID'], axis=1)
patients_names = patient_features.ID

## Select genes in common between embedding and mutation map
common_genes = set(mutation_map.Gene.unique()).intersection(gene_encoding.index)

mutation_map_filt = mutation_map[mutation_map.Gene.isin(common_genes)]

gene_embed = gene_encoding.loc[list(common_genes), :].copy()
gene_embed.loc["Gene_0"] = 0  # Add dummy gene row with zeros
gene_x = torch.from_numpy(csr_matrix(gene_embed).todense()).to(torch.float)
gene_mapping = {index: i for i, index in enumerate(gene_embed.index)}

## Define graph list
graph_list = [None] * len(patients_names)

for i, patient in enumerate(patients_names):
    # Define patient features
    patient_feat = patient_features.loc[patient_features.ID == patient]
    patient_feat_sel = patient_feat.drop(['OS_YEARS', 'OS_STATUS', 'ID'], axis=1)
    patient_x = torch.from_numpy(csr_matrix(patient_feat_sel).todense()).to(torch.float)
    patient_y = torch.from_numpy(csr_matrix(patient_feat[['OS_YEARS', 'OS_STATUS']]).todense()).to(torch.float)
    # Define gene features   
    patient_genes = pd.concat([pd.Series('Gene_0'), mutation_map_filt.loc[mutation_map_filt.ID == patient].Gene])
    indeces = list()
    for gene in patient_genes:
        indeces.append(gene_mapping[gene])
    gene_sel_x = gene_x[indeces, :]
    # Index
    src = [0] * gene_sel_x.shape[0]
    dst = list(range(gene_sel_x.shape[0]))
    edge_index = torch.tensor([dst, src])
   # Generate Object
    data = HeteroData()
    data['patient'].x = patient_x
    data['patient'].y = patient_y
    data['gene'].x = gene_sel_x 
    data['gene', 'patient'].edge_index = edge_index
    graph_list[i] = data
# Save graph list
torch.save(graph_list, "results/gnn/preprocess/boolean_scGPT.pt" )



## Define graph list with reduced gene embeddings
gene_encoding_red = gene_encoding_red.loc[list(common_genes), :].copy()
gene_encoding_red.loc["Gene_0"] = 0  # Add dummy gene row with zeros
gene_red_x = torch.from_numpy(csr_matrix(gene_encoding_red).todense()).to(torch.float)
gene_red_mapping = {index: i for i, index in enumerate(gene_encoding_red.index)}
# Define graph list
graph_list_red = [None] * len(patients_names)

for i, patient in enumerate(patients_names):
    # Define patient features
    patient_feat = patient_features.loc[patient_features.ID == patient]
    patient_feat_sel = patient_feat.drop(['OS_YEARS', 'OS_STATUS', 'ID'], axis=1)
    patient_x = torch.from_numpy(csr_matrix(patient_feat_sel).todense()).to(torch.float)
    patient_y = torch.from_numpy(csr_matrix(patient_feat[['OS_YEARS', 'OS_STATUS']]).todense()).to(torch.float)
    # Define gene features   
    patient_genes = pd.concat([pd.Series('Gene_0'), mutation_map_filt.loc[mutation_map_filt.ID == patient].Gene])
    indeces = list()
    for gene in patient_genes:
        indeces.append(gene_red_mapping[gene])
    gene_sel_x = gene_red_x[indeces, :]
    # Index
    src = [0] * gene_sel_x.shape[0]
    dst = list(range(gene_sel_x.shape[0]))
    edge_index = torch.tensor([dst, src])
   # Generate Object
    data = HeteroData()
    data['patient'].x = patient_x
    data['patient'].y = patient_y
    data['gene'].x = gene_sel_x 
    data['gene', 'patient'].edge_index = edge_index
    graph_list_red[i] = data
# Save graph list
torch.save(graph_list_red, "results/gnn/preprocess/boolean_scGPT_reduced.pt" )


