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

patient_x = torch.from_numpy(csr_matrix(patient_features_sel).todense()).to(torch.float)
patient_mapping = {index: i for i, index in enumerate(patient_features.ID)}

genes = pd.DataFrame({'Gene' : mutation_map['Gene'].unique()})
gene_encoding = pd.get_dummies(genes, columns=['Gene']) 
gene_x = torch.from_numpy(csr_matrix(gene_encoding).todense()).to(torch.float)
gene_mapping = {index: i for i, index in enumerate(genes.Gene)}

# Index
src = [patient_mapping[index] for index in mutation_map['ID']]
dst = [gene_mapping[index] for index in mutation_map['Gene']]
edge_index = torch.tensor([src, dst])

# Generate Object
data = HeteroData()
data['patient'].x = patient_x
data['gene'].x = gene_x 
data['patient', 'interaction', 'gene'].edge_index = edge_index
    
torch.save(data, os.path.join('results/mutations/gnn_graph.pt'))

   