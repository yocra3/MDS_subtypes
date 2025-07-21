"""
Compare the performance of different models 
docker run --rm --gpus all -v $PWD:$PWD -w $PWD -it mds_subtypes_python:1.5 python
"""
import os
import torch
import lightning as L
from lightning.pytorch.loggers import TensorBoardLogger
import torch.nn.functional as F
from torch import nn
from torch.utils.data import TensorDataset, DataLoader,random_split
from torchinfo import summary
import sys
from pycox.models.loss import cox_ph_loss
import random
import torchtuples as tt
import numpy as np
from lifelines.utils import concordance_index
import pandas as pd
from typing import Optional, Union, Tuple
from torch_geometric.data.batch import Batch
from torch_geometric.data import HeteroData

sys.path.append('./scripts/3.GeneticScores/pytorch_files/')
from GNNPatientSAGE import PatientGNNSAGE
from GNNPatientSAGE_v2 import PatientGNNSAGE_v2
from GNNPatientSAGE_v3 import PatientGNNSAGE_v3
from GNNPatientSAGE_v4 import PatientGNNSAGE_v4

from BasicNN import BasicNN

## NN models
patient_vals = pd.read_csv("results/gnn/preprocess/patient_variables.tsv", sep = "\t")
new_cindex = pd.read_csv("results/gnn/patient_recomputedIPSSM.tsv", sep = "\t")

patient_vals = patient_vals.merge(new_cindex[['ID', 'IPSSM_SCORE_train']], on="ID", how="left")

test_patients = patient_vals[patient_vals.train == False]

clin_kar_vars = test_patients[["SEX", "BM_BLAST", "WBC", "ANC", 
    "HB", "PLT", "plus8", "del7", "del20q", "del7q", "delY",
    "del5q", "complex", "LFS_YEARS", 
    "LFS_STATUS", "train"]].copy()
all_vars = test_patients.drop(['OS_YEARS', 'OS_STATUS', 'AMLt_YEARS', 'AMLt_STATUS', 'ID', "IPSSM_SCORE"], axis=1).copy()

clin_kar_dim = clin_kar_vars.shape[1] - 3 # Exclude 'LFS_YEARS', 'LFS_STATUS' and train columns
all_dim = all_vars.shape[1] - 3 
hidden_dim = 16
learning_rate = 0.001

class EvalNN(L.LightningModule):
    def __init__(self, ck_path, inp_dim, hidden_dim, learning_rate=0.001):
        super().__init__()
        # init the pretrained LightningModule
        self.model = BasicNN(inp_dim, hidden_dim, learning_rate)
        checkpoint = torch.load(ck_path, map_location='cpu')
        self.model.load_state_dict(checkpoint['state_dict'])       
    def forward(self, x):
        out = self.model(x)
        return out

# Predict with clinical and karyotype NN
ck_clin_kar = "results/gnn/model_checkpoints/NN_clin_kar-v1.ckpt"
clin_NN = EvalNN(ck_clin_kar, clin_kar_dim, hidden_dim, learning_rate)
clin_NN.eval()


clin_kar_vars_X = clin_kar_vars[clin_kar_vars.train == False].drop(columns=['LFS_YEARS', 'LFS_STATUS', 'train']).values
clin_kar_vars_X = torch.tensor(clin_kar_vars_X, dtype=torch.float32)
test_patients['clin_kar_pred'] = clin_NN(clin_kar_vars_X).detach().numpy()

# Predict with all variables NN
ck_all = "results/gnn/model_checkpoints/NN_all_features-v1.ckpt"
all_NN = EvalNN(ck_all, all_dim, hidden_dim, learning_rate)
all_NN.eval()

all_vars_X = all_vars.drop(columns=['LFS_YEARS', 'LFS_STATUS', 'train']).values
all_vars_X = torch.tensor(all_vars_X, dtype=torch.float32)
test_patients['all_pred'] = all_NN(all_vars_X).detach().numpy()


# GNN models
train_graphs, val_graphs = torch.load("results/gnn/preprocess/graphs_genesEncScGPT.pt")
val_graph_list = Batch.from_data_list(val_graphs)

class EvalGNN(L.LightningModule):
    def __init__(self, ck_path, patient_feat_dim, gene_feat_dim, hidden_gene_dim, hidden_dim, out_dim, use_vaf = False, hidden_vaf = 0, learning_rate=0.001):
        super().__init__()
        # init the pretrained LightningModule
        self.model =  PatientGNNSAGE(patient_feat_dim, gene_feat_dim, hidden_gene_dim, hidden_dim, out_dim, use_vaf, hidden_vaf, learning_rate)
        checkpoint = torch.load(ck_path, map_location='cpu')
        self.model.load_state_dict(checkpoint['state_dict'])       
    def forward(self, x):
        out = self.model(x)
        return out

class EvalGNN2(L.LightningModule):
    def __init__(self, ck_path, patient_feat_dim, gene_feat_dim, hidden_gene_dim, hidden_dim, out_dim, use_vaf = False, hidden_vaf = 0, learning_rate=0.001):
        super().__init__()
        # init the pretrained LightningModule
        self.model =  PatientGNNSAGE_v2(patient_feat_dim, gene_feat_dim, hidden_gene_dim, hidden_dim, out_dim, use_vaf, hidden_vaf, learning_rate)
        checkpoint = torch.load(ck_path, map_location='cpu')
        self.model.load_state_dict(checkpoint['state_dict'])       
    def forward(self, x):
        out = self.model(x)
        return out


class EvalGNN3(L.LightningModule):
    def __init__(self, ck_path, patient_feat_dim, gene_feat_dim, hidden_gene_dim, hidden_dim, out_dim, use_vaf = False, hidden_vaf = 0, learning_rate=0.001):
        super().__init__()
        # init the pretrained LightningModule
        self.model =  PatientGNNSAGE_v3(patient_feat_dim, gene_feat_dim, hidden_gene_dim, hidden_dim, out_dim, use_vaf, hidden_vaf, learning_rate)
        checkpoint = torch.load(ck_path, map_location='cpu')
        self.model.load_state_dict(checkpoint['state_dict'])       
    def forward(self, x):
        out = self.model(x)
        return out

class EvalGNN4(L.LightningModule):
    def __init__(self, ck_path,patient_feat_dim, gene_feat_dim, hidden_gene_dim, hidden_dim, out_dim, 
                 dropout=0.3, l2_reg=1e-4, use_vaf=False, hidden_vaf=0, learning_rate=0.001):
        super().__init__()
        # init the pretrained LightningModule
        self.model =  PatientGNNSAGE_v4(patient_feat_dim, gene_feat_dim, hidden_gene_dim, hidden_dim, out_dim, 
                 dropout=dropout, l2_reg=l2_reg, use_vaf=use_vaf, hidden_vaf=hidden_vaf, learning_rate=learning_rate)
        checkpoint = torch.load(ck_path, map_location='cpu')
        self.model.load_state_dict(checkpoint['state_dict'])       
    def forward(self, x):
        out = self.model(x)
        return out

patient_feat_dim = 13
gene_feat_dim = 16
hidden_gene_dim = 32
hidden_dim = 16
out_dim = 1  # Single output for survival regression
learning_rate = 0.001
hidden_vaf = 8

# Predict with GNN model
## Boolean
ck_gnn_bool = "results/gnn/model_checkpoints/GNN_boolean-v2.ckpt"

gnn_model_bool = EvalGNN(ck_gnn_bool, patient_feat_dim, gene_feat_dim, hidden_gene_dim, hidden_dim, out_dim, learning_rate)
gnn_model_bool.eval()
test_patients['gnn_bool_pred'] = gnn_model_bool(val_graph_list).detach().numpy()

## Optimal parameters
ck_gnn_bool_optim = "results/gnn/model_checkpoints/GNN_boolean_optim.ckpt"

hidden_gene_dim = 10
hidden_dim = 60
learning_rate = 0.0001
gnn_model_bool_optim = EvalGNN(ck_gnn_bool_optim, patient_feat_dim, gene_feat_dim, hidden_gene_dim, hidden_dim, out_dim, learning_rate)
gnn_model_bool_optim.eval()
test_patients['gnn_bool_optim_pred'] = gnn_model_bool_optim(val_graph_list).detach().numpy()

## Boolean v2
ck_gnn_bool_v2 = "results/gnn/model_checkpoints/GNN_boolean_v2-v1.ckpt"

gnn_model_bool_v2 = EvalGNN2(ck_gnn_bool_v2, patient_feat_dim, gene_feat_dim, hidden_gene_dim, hidden_dim, out_dim, learning_rate)
gnn_model_bool_v2.eval()
test_patients['gnn_bool_v2_pred'] = gnn_model_bool_v2(val_graph_list).detach().numpy()



## Boolean v3
ck_gnn_bool_v3 = "results/gnn/model_checkpoints/GNN_boolean_v3.ckpt"

gnn_model_bool_v3 = EvalGNN3(ck_gnn_bool_v3, patient_feat_dim, gene_feat_dim, hidden_gene_dim, hidden_dim, out_dim, learning_rate)
gnn_model_bool_v3.eval()
test_patients['gnn_bool_v3_pred'] = gnn_model_bool_v3(val_graph_list).detach().numpy()


## Boolean v4 optim
ck_gnn_bool_v4_optim = "results/gnn/model_checkpoints/GNN_boolean_v4_optim-v1.ckpt"

hidden_gene_dim_optim = 20
hidden_dim_optim = 80
learning_rate_optim = 0.0001
dropout_optim = 0.4     # Más dropout para modelo más grande
l2_reg_optim = 5e-5     # L2 un poco menor
patient_feat_dim = 13
gene_feat_dim = 16
out_dim = 1

gnn_model_bool_v4_optim = EvalGNN4(ck_gnn_bool_v4_optim, patient_feat_dim, gene_feat_dim, hidden_gene_dim_optim, hidden_dim_optim, out_dim, dropout_optim, 
    l2_reg_optim, False, 0, learning_rate_optim)
gnn_model_bool_v4_optim.eval() 
test_patients['gnn_bool_v4_optim_pred'] = gnn_model_bool_v4_optim(val_graph_list).detach().numpy()


# VAF
ck_gnn_vaf = "results/gnn/model_checkpoints/GNN_vaf-v2.ckpt"

gnn_model_vaf = EvalGNN(ck_gnn_vaf, patient_feat_dim, gene_feat_dim, hidden_gene_dim, hidden_dim, out_dim, use_vaf = True, hidden_vaf = hidden_vaf, learning_rate = learning_rate)
gnn_model_vaf.eval()
test_patients['gnn_vaf_pred'] = gnn_model_vaf(val_graph_list).detach().numpy()

## Optimal hyperparameters
ck_gnn_vaf_optim = "results/gnn/model_checkpoints/GNN_vaf_optim.ckpt"

hidden_gene_dim = 15
hidden_dim = 70
learning_rate = 0.0001
hidden_vaf = 4


gnn_model_vaf_optim = EvalGNN(ck_gnn_vaf_optim, patient_feat_dim, gene_feat_dim, hidden_gene_dim, hidden_dim, out_dim, use_vaf = True, hidden_vaf = hidden_vaf, learning_rate = learning_rate)
gnn_model_vaf_optim.eval()
test_patients['gnn_vaf_optim_pred'] = gnn_model_vaf_optim(val_graph_list).detach().numpy()

## Version 4
ck_gnn_vaf_v4 = "results/gnn/model_checkpoints/GNN_vaf_v4_optim-v2.ckpt"

gnn_model_vaf_v4 = EvalGNN4(ck_gnn_vaf_v4, patient_feat_dim, gene_feat_dim, hidden_gene_dim_optim, hidden_dim_optim, out_dim, use_vaf = True, hidden_vaf = 6, learning_rate = learning_rate)
gnn_model_vaf_v4.eval()
test_patients['gnn_vaf_v4_optim_pred'] = gnn_model_vaf_v4(val_graph_list).detach().numpy()

concordance_index(test_patients['LFS_YEARS'].values, -test_patients['clin_kar_pred'].values, test_patients['LFS_STATUS'].values)
# 0.7435147854859695
concordance_index(test_patients['LFS_YEARS'].values, -test_patients['all_pred'].values, test_patients['LFS_STATUS'].values)
# 0.7568787445588447
concordance_index(test_patients['LFS_YEARS'].values, -test_patients['gnn_bool_pred'].values, test_patients['LFS_STATUS'].values)
# 0.7675281494450423
concordance_index(test_patients['LFS_YEARS'].values, -test_patients['gnn_bool_optim_pred'].values, test_patients['LFS_STATUS'].values)
# 0.7659540292657854
concordance_index(test_patients['LFS_YEARS'].values, -test_patients['gnn_vaf_pred'].values, test_patients['LFS_STATUS'].values)
# 0.7561238093708338
concordance_index(test_patients['LFS_YEARS'].values, -test_patients['gnn_vaf_optim_pred'].values, test_patients['LFS_STATUS'].values)
# 0.7605570457940473

concordance_index(test_patients['LFS_YEARS'].values, -test_patients['gnn_bool_v2_pred'].values, test_patients['LFS_STATUS'].values)
# 0.7631912877266813
concordance_index(test_patients['LFS_YEARS'].values, -test_patients['gnn_bool_v3_pred'].values, test_patients['LFS_STATUS'].values)
# 0.7613762307852933
concordance_index(test_patients['LFS_YEARS'].values, -test_patients['gnn_bool_v4_optim_pred'].values, test_patients['LFS_STATUS'].values)
# 0.764219284578441
concordance_index(test_patients['LFS_YEARS'].values, -test_patients['gnn_vaf_v4_optim_pred'].values, test_patients['LFS_STATUS'].values)
# 0.7545657516423856


concordance_index(test_patients['LFS_YEARS'].values, -test_patients['IPSSM_SCORE'].values, test_patients['LFS_STATUS'].values)
# 0.7752220633824309
concordance_index(test_patients['LFS_YEARS'].values, -test_patients['IPSSM_SCORE_train'].values, test_patients['LFS_STATUS'].values)
# 0.7656970300528455

test_os = test_patients.dropna(subset = ['OS_YEARS', 'OS_STATUS'])
concordance_index(test_os['OS_YEARS'].values, -test_os['gnn_bool_pred'].values, test_os['OS_STATUS'].values)
# 0.7704047085785718
concordance_index(test_os['OS_YEARS'].values, -test_os['gnn_vaf_pred'].values, test_os['OS_STATUS'].values)
# 0.7593624667951456

concordance_index(test_os['OS_YEARS'].values, -test_os['gnn_bool_optim_pred'].values, test_os['OS_STATUS'].values)
# 0.7670712016250847
concordance_index(test_os['OS_YEARS'].values, -test_os['gnn_vaf_optim_pred'].values, test_os['OS_STATUS'].values)
# 0.7616195287949025

concordance_index(test_os['OS_YEARS'].values, -test_os['IPSSM_SCORE'].values, test_os['OS_STATUS'].values)
# 0.7781308054238936
concordance_index(test_os['OS_YEARS'].values, -test_os['IPSSM_SCORE_train'].values, test_os['OS_STATUS'].values)
# 0.7713248953938573

concordance_index(test_patients['AMLt_YEARS'].values, -test_patients['gnn_bool_pred'].values, test_patients['AMLt_STATUS'].values)
# 0.8202218304686353
concordance_index(test_patients['AMLt_YEARS'].values, -test_patients['gnn_vaf_pred'].values, test_patients['AMLt_STATUS'].values)
# 0.8058248861466137

concordance_index(test_patients['AMLt_YEARS'].values, -test_patients['gnn_bool_optim_pred'].values, test_patients['AMLt_STATUS'].values)
# 0.8124724548259145
concordance_index(test_patients['AMLt_YEARS'].values, -test_patients['gnn_vaf_optim_pred'].values, test_patients['AMLt_STATUS'].values)
# 0.8121051858381079

concordance_index(test_patients['AMLt_YEARS'].values, -test_patients['IPSSM_SCORE'].values, test_patients['AMLt_STATUS'].values)
# 0.8307624504186867
concordance_index(test_patients['AMLt_YEARS'].values, -test_patients['IPSSM_SCORE_train'].values, test_patients['AMLt_STATUS'].values)
# 0.8157044219186131

## Compute for all patients
all_graphs = train_graphs + val_graphs
all_graph_list = Batch.from_data_list(all_graphs)

# Predict with GNN model
all_pred_bool = gnn_model_bool(all_graph_list).detach().numpy()
all_pred_vaf = gnn_model_vaf(all_graph_list).detach().numpy()
all_pred_bool_optim = gnn_model_bool_optim(all_graph_list).detach().numpy()
all_pred_vaf_optim = gnn_model_vaf_optim(all_graph_list).detach().numpy()
all_pred_bool_v4_optim = gnn_model_bool_v4_optim(all_graph_list).detach().numpy()

all_pred_pd = pd.DataFrame({
    'ID': all_graph_list['patient'].ID,
    #'gnn_bool_pred': all_pred_bool.flatten(),
  #  'gnn_vaf_pred': all_pred_vaf.flatten(),
    'gnn_bool_optim_pred': all_pred_bool_optim.flatten(),
    'gnn_vaf_optim_pred': all_pred_vaf_optim.flatten(),
    'gnn_bool_v4_optim_pred': all_pred_bool_v4_optim.flatten()
})

## Merge predictions with patient values
all_patient_vals = patient_vals.merge(all_pred_pd, on='ID', how='left').dropna(subset=[ 'IPSSM_SCORE'])
all_patient_vals.to_csv('results/gnn/all_patients_predictions.tsv', sep='\t', index=False)


concordance_index(all_patient_vals['LFS_YEARS'].values, -all_patient_vals['gnn_bool_pred'].values, all_patient_vals['LFS_STATUS'].values)
# 0.7566120080609822
concordance_index(all_patient_vals['LFS_YEARS'].values, -all_patient_vals['gnn_vaf_pred'].values, all_patient_vals['LFS_STATUS'].values)        
# 0.7588433852159366
concordance_index(all_patient_vals['LFS_YEARS'].values, -all_patient_vals['gnn_bool_optim_pred'].values, all_patient_vals['LFS_STATUS'].values)        
# 0.7588651966632451
concordance_index(all_patient_vals['LFS_YEARS'].values, -all_patient_vals['gnn_vaf_optim_pred'].values, all_patient_vals['LFS_STATUS'].values)        
# 0.7502305073408739
concordance_index(all_patient_vals['LFS_YEARS'].values, -all_patient_vals['gnn_bool_v4_optim_pred'].values, all_patient_vals['LFS_STATUS'].values)        
# 0.7503554274481863
concordance_index(all_patient_vals['LFS_YEARS'].values, -all_patient_vals['IPSSM_SCORE'].values, all_patient_vals['LFS_STATUS'].values)
# 0.7526063027151286


all_os = all_patient_vals.dropna(subset = ['OS_YEARS', 'OS_STATUS'])


concordance_index(all_os['OS_YEARS'].values, -all_os['gnn_bool_pred'].values, all_os['OS_STATUS'].values)
# 0.7534992800009596
concordance_index(all_os['OS_YEARS'].values, -all_os['gnn_vaf_pred'].values, all_os['OS_STATUS'].values)
# 0.7546380143971587
concordance_index(all_os['OS_YEARS'].values, -all_os['gnn_bool_optim_pred'].values, all_os['OS_STATUS'].values)
# 0.7570664677997098
concordance_index(all_os['OS_YEARS'].values, -all_os['gnn_vaf_optim_pred'].values, all_os['OS_STATUS'].values)
# 0.7469582598106483

concordance_index(all_os['OS_YEARS'].values, -all_os['IPSSM_SCORE'].values, all_os['OS_STATUS'].values)
# 0.7482917220218109

all_amlt = all_patient_vals.dropna(subset = ['AMLt_YEARS', 'AMLt_STATUS'])

concordance_index(all_amlt['AMLt_YEARS'].values, -all_amlt['gnn_bool_pred'].values, all_amlt['AMLt_STATUS'].values)
# 0.8093647398406544
concordance_index(all_amlt['AMLt_YEARS'].values, -all_amlt['gnn_vaf_pred'].values, all_amlt['AMLt_STATUS'].values)
# 0.8126110796340624
concordance_index(all_amlt['AMLt_YEARS'].values, -all_amlt['gnn_bool_optim_pred'].values, all_amlt['AMLt_STATUS'].values)
# 0.8057289869168262
concordance_index(all_amlt['AMLt_YEARS'].values, -all_amlt['gnn_vaf_optim_pred'].values, all_amlt['AMLt_STATUS'].values)
# 0.8004783691156631

concordance_index(all_amlt['AMLt_YEARS'].values, -all_amlt['IPSSM_SCORE'].values, all_amlt['AMLt_STATUS'].values)
# 0.8217355275809693


## Define a PD with the results

models = ['gnn_bool', 'gnn_vaf', 'gnn_bool_optim', 'gnn_vaf_optim', 'IPSSM_SCORE']
measures = ['OS', 'LFS', 'AMLt']
datasets = ['validation', 'all']
results = []
for measure in measures:
    for model in models:
        for dataset in datasets:
            if dataset == 'validation':
                if measure == "OS":
                    df = test_os
                else:
                    df = test_patients
            else:
                if measure == "OS":
                    df = all_os
                elif measure == "LFS":
                    df = all_patient_vals
                elif measure == "AMLt":
                    df = all_amlt
            if model == "IPSSM_SCORE":
                suffix = ""
            else:
                suffix = "_pred"
            c_index = concordance_index(df[measure + '_YEARS'].values, -df[model + suffix].values, df[measure + '_STATUS'].values)
            results.append({
                "model": model,
                "measure": measure,
                "dataset": dataset,
                "c_index": c_index
             })
results_df = pd.DataFrame(results)
results_df.to_csv('results/gnn/models_cindex.tsv', sep='\t', index=False)

import matplotlib.pyplot as plt
import seaborn as sns


results_df['model'] = results_df['model'].astype(str)
results_df['measure'] = results_df['measure'].astype(str)
results_df['dataset'] = results_df['dataset'].astype(str)

# Crear el FacetGrid
g = sns.FacetGrid(
    results_df,
    col="measure",
    sharey=True,
    height=4,
    aspect=1.2
)

# Dibujar los puntos por dataset en cada panel
g.map_dataframe(
    sns.barplot,
    x="dataset",
    y="c_index",
    hue="model",
    dodge=True,
    errorbar = None,
    palette="Set1"
)

# Ajustar leyenda y títulos
g.add_legend(title="Dataset")
g.set_axis_labels("Model", "C-index")
for ax in g.axes.flat:
    ax.set_ylim(0.7, 0.85)


g.set_titles(col_template="{col_name}")

g.add_legend(title="Model")
plt.legend(
    title="Model",
    bbox_to_anchor=(0.5, -0.15),
    loc='upper center',
    ncol=len(results_df['model'].unique())
)

plt.tight_layout()
plt.show()
plt.savefig("figures/panel_cindex_gnn_models.png", dpi=300)