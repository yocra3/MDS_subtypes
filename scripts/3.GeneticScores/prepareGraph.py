"""
This is the data preprocessing that should return 
the .pt with the data to run the model 

docker run --gpus all -v $PWD:$PWD -w $PWD -it pyg/gennius python

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


### Generate heterodata object
# Features
patient_features_sel = patient_features.drop(['OS_YEARS', 'OS_STATUS', 'ID'], axis=1)
patients_names = patient_features.ID


genes = pd.DataFrame({'Gene' : mutation_map['Gene'].unique()})
gene_encoding = pd.get_dummies(genes, columns=['Gene']) 
gene_x = torch.from_numpy(csr_matrix(gene_encoding).todense()).to(torch.float)
gene_mapping = {index: i for i, index in enumerate(genes.Gene.sort_values())}

## Define graph list
graph_list = [None] * len(patients_names)

for i, patient in enumerate(patients_names):
    # Define patient features
    patient_feat = patient_features.loc[patient_features.ID == patient]
    patient_feat_sel = patient_feat.drop(['OS_YEARS', 'OS_STATUS', 'ID'], axis=1)
    patient_x = torch.from_numpy(csr_matrix(patient_feat_sel).todense()).to(torch.float)
    patient_y = torch.from_numpy(csr_matrix(patient_feat[['OS_YEARS', 'OS_STATUS']]).todense()).to(torch.float)
    # Define gene features   
    patient_genes = mutation_map.loc[mutation_map.ID == patient].Gene
    indeces = list()
    for gene in patient_genes:
        indeces.append(gene_mapping[gene])
    genes_sel = gene_encoding.iloc[:, indeces]
    genes_sel['Gene_0'] = False
    gene_x = torch.from_numpy(csr_matrix(genes_sel).todense()).to(torch.float)
    # Index
    src = [0] * gene_x.shape[1]
    dst = list(range(gene_x.shape[1]))
    edge_index = torch.tensor([src, dst])
   # Generate Object
    data = HeteroData()
    data['patient'].x = patient_x
    data['patient'].y = patient_y
    data['gene'].x = gene_x 
    data['patient', 'interaction', 'gene'].edge_index = edge_index
    graph_list[i] = data


torch.save(graph_list, "results/gnn/preprocess/boolean_dummy.pt" )







# patient_x = torch.from_numpy(csr_matrix(patient_features_sel).todense()).to(torch.float)
# patient_mapping = {index: i for i, index in enumerate(patient_features.ID)}


# # Index
# src = [patient_mapping[index] for index in mutation_map['ID']]
# dst = [gene_mapping[index] for index in mutation_map['Gene']]
# edge_index = torch.tensor([src, dst])

# # Generate Object
# data = HeteroData()
# data['patient'].x = patient_x
# data['gene'].x = gene_x 
# data['patient', 'interaction', 'gene'].edge_index = edge_index
    
# torch.save(data, os.path.join('results/mutations/gnn_graph.pt'))

   