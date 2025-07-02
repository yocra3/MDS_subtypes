"""
Reduce dimensionality of gene embeddings from scGPT before passing them to the GNN model.
docker run --gpus all -v $PWD:$PWD -w $PWD -it mds_subtypes_python:1.4 python

"""

import os
import torch
import lightning as L
from lightning.pytorch.loggers import TensorBoardLogger
from torch import nn
from torch.utils.data import TensorDataset, DataLoader,random_split
import sys
import random
import numpy as np
import pandas as pd

### Load data
embeddings = pd.read_csv("results/preprocess/gene_embedding_scgpt.tsv", sep="\t", index_col=0)
coding_map = pd.read_csv("results/preprocess/coding_genes.txt", sep="\t", index_col=0)

## Select only coding genes
coding_genes = coding_map.hgnc_symbol.values
intersection_genes = set(coding_genes).intersection(set(embeddings.index))
coding_embeddings = embeddings.loc[list(intersection_genes), :]

# Define autoencoder to reduce dimensionality
## Define the encoder and decoder
# define the LightningModule
class LitAutoEncoder(L.LightningModule):
    def __init__(self, encoder, decoder):
        super().__init__()
        self.encoder = nn.Sequential(nn.Linear(512, 128), nn.ReLU(), nn.Linear(128, 16))
        self.decoder = nn.Sequential(nn.Linear(16, 128), nn.ReLU(), nn.Linear(128, 512))
    def training_step(self, batch, batch_idx):
        x, _ = batch
        x = x.view(x.size(0), -1)
        z = self.encoder(x)
        x_hat = self.decoder(z)
        loss = nn.functional.mse_loss(x_hat, x)
        self.log("train_loss", loss, prog_bar=True, on_epoch=True)
        return loss
    def validation_step(self, batch, batch_idx):
        x, _ = batch
        x = x.view(x.size(0), -1)
        z = self.encoder(x)
        x_hat = self.decoder(z)
        loss = nn.functional.mse_loss(x_hat, x)
        self.log("val_loss", loss, prog_bar=True)
        return loss
    def configure_optimizers(self):
        optimizer = torch.optim.Adam(self.parameters(), lr=1e-3)
        return optimizer
# Convert pandas to tensor dataset
X = torch.FloatTensor(coding_embeddings.values)
dataset = TensorDataset(X, X)  

# Split dataset
train_size = int(0.8 * len(dataset))
val_size = len(dataset) - train_size
train_dataset, val_dataset = random_split(dataset, [train_size, val_size])

# Create data loaders
batch_size = 32
train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
val_loader = DataLoader(val_dataset, batch_size=batch_size)

# train model
logger_ae = TensorBoardLogger("lightning_logs", name="gene_embedding_ae")
trainer_ae = L.Trainer(max_epochs = 100, logger = logger_ae)

autoencoder = LitAutoEncoder(encoder, decoder)

trainer_ae.fit(model=autoencoder, train_dataloaders=train_loader, val_dataloaders=val_loader)

## Final ae (50 epochs)
logger_ae = TensorBoardLogger("lightning_logs", name="gene_embedding_ae")
trainer_ae = L.Trainer(max_epochs = 50, logger = logger_ae)

autoencoder = LitAutoEncoder(encoder, decoder)

trainer_ae.fit(model=autoencoder, train_dataloaders=train_loader, val_dataloaders=val_loader)

## Predict with the trained autoencoder
encoder = nn.Sequential(nn.Linear(512, 128), nn.ReLU(), nn.Linear(128, 16))
decoder = nn.Sequential(nn.Linear(16, 128), nn.ReLU(), nn.Linear(128, 512))
# Load the trained autoencoder
checkpoint = "lightning_logs/gene_embedding_ae/version_6/checkpoints/epoch=49-step=23950.ckpt"
autoencoder_pred = LitAutoEncoder.load_from_checkpoint(checkpoint, encoder=encoder, decoder=decoder)

trained_encoder = autoencoder.encoder
trained_encoder.eval()

reduced_embeddings = trained_encoder(X)
reduced_embeddings_df = pd.DataFrame(reduced_embeddings.detach().numpy(), index=coding_embeddings.index)
reduced_embeddings_df.to_csv("results/preprocess/gene_embedding_reduced.tsv", sep="\t")