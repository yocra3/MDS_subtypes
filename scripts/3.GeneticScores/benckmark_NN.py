"""
Module Summary:
This code will evaluate training a model using just a NN with or without genes.
Use LFS for training
docker run --rm --gpus all -v $PWD:$PWD -w $PWD -it mds_subtypes_python:1.5 python

"""

import os
import torch
import lightning as L
from lightning.pytorch.loggers import TensorBoardLogger
import torch.nn.functional as F
from torch.utils.data import TensorDataset, DataLoader,random_split
from torchinfo import summary
import sys
import random
import torchtuples as tt
import numpy as np
import pandas as pd
from typing import Optional, Union, Tuple

sys.path.append('./scripts/3.GeneticScores/pytorch_files/')
from BasicNN import BasicNN


# Load data
patient_vals = pd.read_csv("results/gnn/preprocess/patient_variables.tsv", sep = "\t")

clin_kar_vars = patient_vals[["SEX", "BM_BLAST", "WBC", "ANC", 
     "HB", "PLT", "plus8", "del7", "del20q", "del7q", "delY",
    "del5q", "complex", "LFS_YEARS", 
    "LFS_STATUS", "train"]].copy()
all_vars = patient_vals.drop(['OS_YEARS', 'OS_STATUS', 'AMLt_YEARS', 'AMLt_STATUS', 'ID', 'IPSSM_SCORE'], axis=1).copy()


## Function to convert DataFrame to DataLoader
def df_to_dataloader(
    df: pd.DataFrame,
    batch_size: int = 32) -> Union[DataLoader, Tuple[DataLoader, DataLoader]]:
    """
    Converts a Pandas DataFrame to one or two PyTorch DataLoaders.
    """
    X = df.drop(columns=['LFS_YEARS', 'LFS_STATUS', 'train']).values
    y = df[['LFS_YEARS', 'LFS_STATUS']].values
    X_tensor = torch.tensor(X, dtype=torch.float32)
    y_tensor = torch.tensor(y, dtype=torch.float32)
    train_dataset = TensorDataset(X_tensor[df['train']], y_tensor[df['train']])
    train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
    val_dataset = TensorDataset(X_tensor[df['train'] == False], y_tensor[df['train'] == False])
    val_loader = DataLoader(val_dataset, batch_size=64, shuffle=False)
    return train_loader, val_loader

## Define dataset
train_loader, test_loader = df_to_dataloader(clin_kar_vars, batch_size=32)
train_all_loader, test_all_loader = df_to_dataloader(all_vars, batch_size=32)

## Define model
clin_kar_dim = clin_kar_vars.shape[1] - 3 # Exclude 'LFS_YEARS', 'LFS_STATUS' and train columns
all_dim = all_vars.shape[1] - 3 
hidden_dim = 16
learning_rate = 0.001

clin_NN = BasicNN(clin_kar_dim, hidden_dim, learning_rate)
all_NN = BasicNN(all_dim, hidden_dim, learning_rate)


# train model
logger_base = TensorBoardLogger("lightning_logs/current", name="clinical_kar")
logger_all = TensorBoardLogger("lightning_logs/current", name="all_features")


def define_checkpoint_callback(model_name: str):
    """
    Define a callback to save the best model based on validation loss.
    """
    from lightning.pytorch.callbacks import ModelCheckpoint
    return ModelCheckpoint(
        dirpath="results/gnn/model_checkpoints/",
        filename=model_name,
        save_top_k=1,
        monitor="val_loss", 
        mode="min"
    )


# Create trainer with the callback
trainer_base = L.Trainer(max_epochs = 30, logger = logger_base, callbacks=[define_checkpoint_callback("NN_clin_kar")])
trainer_all = L.Trainer(max_epochs = 30, logger = logger_all, callbacks=[define_checkpoint_callback("NN_all_features")])


trainer_base.fit(model=clin_NN, train_dataloaders=train_loader, val_dataloaders=test_loader)
trainer_all.fit(model=all_NN, train_dataloaders=train_all_loader, val_dataloaders=test_all_loader)


## Train model with logPLT
patient_vals_log = pd.read_csv("results/gnn/preprocess/patient_variables_logPLT.tsv", sep = "\t")

all_vars_log = patient_vals_log.drop(['OS_YEARS', 'OS_STATUS', 'AMLt_YEARS', 'AMLt_STATUS', 'ID', 'IPSSM_SCORE'], axis=1).copy()

train_log_loader, test_log_loader = df_to_dataloader(all_vars_log, batch_size=32)
all_NN2 = BasicNN(all_dim, hidden_dim, learning_rate)

logger_all_log = TensorBoardLogger("lightning_logs/current", name="all_features_logPLT")
trainer_all_log = L.Trainer(max_epochs = 30, logger = logger_all_log, callbacks=[define_checkpoint_callback("NN_all_features_logPLT")])
trainer_all_log.fit(model=all_NN2, train_dataloaders=train_log_loader, val_dataloaders=test_log_loader)