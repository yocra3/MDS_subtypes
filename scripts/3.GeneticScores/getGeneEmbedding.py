"""
Get gene embeddings using the pre-trained model from scGPT
docker run --gpus all -v $PWD:$PWD -w $PWD -it mds_subtypes_python:1.4 python

"""
import copy
import json
import os
from pathlib import Path
import sys
import warnings

import torch
from anndata import AnnData
import numpy as np
import pandas as pd

from torchtext.vocab import Vocab
from torchtext._torchtext import (
    Vocab as VocabPybind,
)

import scgpt as scg
from scgpt.tasks import GeneEmbedding
from scgpt.tokenizer.gene_tokenizer import GeneVocab
from scgpt.model import TransformerModel
from scgpt.preprocess import Preprocessor
from scgpt.utils import set_seed

# Specify model path; here we load the pre-trained scGPT blood model
model_dir = Path("data/scGPT/")
model_config_file = model_dir / "args.json"
model_file = model_dir / "best_model.pt"
vocab_file = model_dir / "vocab.json"

# Retrieve model parameters from config files
with open(model_config_file, "r") as f:
    model_configs = json.load(f)

embsize = model_configs["embsize"]
nhead = model_configs["nheads"]
d_hid = model_configs["d_hid"]
nlayers = model_configs["nlayers"]
n_layers_cls = model_configs["n_layers_cls"]
n_input_bins = 51
pad_token = "<pad>"
pad_value = model_configs["pad_value"]

special_tokens = [pad_token, "<cls>", "<eoc>"]
vocab = GeneVocab.from_file(vocab_file)
for s in special_tokens:
    if s not in vocab:
        vocab.append_token(s)
        
gene2idx = vocab.get_stoi()

## Load model
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")


ntokens = len(vocab)  # size of vocabulary
model = TransformerModel(
    ntokens,
    embsize,
    nhead,
    d_hid,
    nlayers,
    vocab=vocab,
    pad_value=pad_value,
    n_input_bins=n_input_bins,
)

try:
    model.load_state_dict(torch.load(model_file))
    print(f"Loading all model params from {model_file}")
except:
    # only load params that are in the model and match the size
    model_dict = model.state_dict()
    pretrained_dict = torch.load(model_file)
    pretrained_dict = {
        k: v
        for k, v in pretrained_dict.items()
        if k in model_dict and v.shape == model_dict[k].shape
    }
    for k, v in pretrained_dict.items():
        print(f"Loading params {k} with shape {v.shape}")
        model_dict.update(pretrained_dict)
        model.load_state_dict(model_dict)

model.to(device)


## Get embeddings
gene_ids = np.array([id for id in gene2idx.values()])
gene_embeddings = model.encoder(torch.tensor(gene_ids, dtype=torch.long).to(device))
gene_embeddings = gene_embeddings.detach().cpu().numpy()

gene_embeddings = {gene: gene_embeddings[i] for i, gene in enumerate(gene2idx.keys())}

embedding_df = pd.DataFrame(gene_embeddings)
gene_embedding_df = embedding_df.T
gene_embedding_df.to_csv('results/preprocess/gene_embedding_scgpt.tsv', sep='\t', index=True)
