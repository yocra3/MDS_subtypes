"""
This is the data preprocessing that should return 
the .pt with the data to run the model 

docker run --rm --gpus all -v $PWD:$PWD -w $PWD -it mds_subtypes_python:1.5 python

"""

import os
import pandas as pd

import torch
from torch_geometric.data import HeteroData
from scipy.sparse import csr_matrix

## Read data
patient_vals = pd.read_csv("results/gnn/preprocess/patient_variables.tsv", sep = "\t")

mutation_map = pd.read_csv("results/gnn/preprocess/mutation_maf_gnn.tsv", 
    sep = "\t")

gene_encoding_red = pd.read_csv("results/preprocess/gene_embedding_reduced.tsv", 
    sep = "\t", index_col=0)


# Process data
## Select columns for defining patient features
patient_clin_vars = patient_vals[["SEX", "BM_BLAST", "WBC", "ANC", 
    "HB", "PLT", "plus8", "del7", "del20q", "del7q", "delY",
    "del5q", "complex", "ID", "LFS_YEARS", "LFS_STATUS", "train"]].copy()

## Remove rows with NaaN
mutation_maf = mutation_map.dropna(subset=["VAF"])

## Process genetic data
common_genes = set(mutation_maf.GENE.unique()).intersection(gene_encoding_red.index)

mutation_map_filt = mutation_map[mutation_map.GENE.isin(common_genes)]

gene_embed = gene_encoding_red.loc[list(common_genes), :].copy()
gene_embed.loc["Gene_0"] = 0  # Add dummy gene row with zeros
gene_x = torch.from_numpy(csr_matrix(gene_embed).todense()).to(torch.float)
gene_mapping = {index: i for i, index in enumerate(gene_embed.index)}


def define_graph(patient_df, patient_id, mutation_map, gene_embedding, gene_mapping):
    # Define patient features
    patient_feats = patient_df.loc[patient_df.ID == patient_id]
    patient_feats_sel = patient_feats.drop(['LFS_YEARS', 'LFS_STATUS', 'ID', 'train'], axis=1)
    patient_x = torch.from_numpy(csr_matrix(patient_feats_sel).todense()).to(torch.float)
    patient_y = torch.from_numpy(csr_matrix(patient_feats[['LFS_YEARS', 'LFS_STATUS']]).todense()).to(torch.float)
    ## Get genes in patient
    patient_genes = mutation_map.loc[mutation_map.ID == patient_id]
    patient_genes = patient_genes.drop("ID", axis=1)
    gene_0 = pd.DataFrame([{'GENE': 'Gene_0', 'VAF': 0}])
    patient_genes = pd.concat([patient_genes, gene_0], ignore_index=True)
    indeces = [gene_mapping[gene] for gene in patient_genes.GENE]
    gene_sel_x = gene_x[indeces, :]
    # Index
    src = [0] * gene_sel_x.shape[0]
    dst = list(range(gene_sel_x.shape[0]))
    edge_index = torch.tensor([dst, src])
    # Weights
    weights = torch.tensor(patient_genes.VAF, dtype = torch.float)
    # Generate Object
    graph = HeteroData()
    graph['patient'].x = patient_x
    graph['patient'].y = patient_y
    graph['patient'].ID = patient_id
    graph['gene'].x = gene_sel_x 
    graph['gene'].weights = weights
    graph['gene', 'patient'].edge_index = edge_index
    return graph

### Generate heterodata object
# Features
patients_names = patient_clin_vars.ID
graph_list = [define_graph(patient_clin_vars, patient, mutation_map_filt, gene_embed, gene_mapping) for patient in patients_names]

train_list = [x for x, m in zip(graph_list, patient_clin_vars['train']) if m]
test_list = [x for x, m in zip(graph_list, patient_clin_vars['train']) if not m]

# Save graph list
torch.save([train_list, test_list], "results/gnn/preprocess/graphs_genesEncScGPT.pt" )

## Same but using logPLT
patient_vals2 = pd.read_csv("results/gnn/preprocess/patient_variables_logPLT.tsv", sep = "\t")
patient_clin_vars2 = patient_vals2[["SEX", "BM_BLAST", "WBC", "ANC", 
    "HB", "logPLT", "plus8", "del7", "del20q", "del7q", "delY",
    "del5q", "complex", "ID", "LFS_YEARS", "LFS_STATUS", "train"]].copy()

graph_list2 = [define_graph(patient_clin_vars2, patient, mutation_map_filt, gene_embed, gene_mapping) for patient in patients_names]

