"""
Train model
docker run --gpus all -v $PWD:$PWD -w $PWD -it mds_subtypes_python:1.5 python
"""
import lightning as L
import sys
import torch
from torch.utils.data import random_split
from torch_geometric.loader import DataLoader
from lightning.pytorch.loggers import TensorBoardLogger

sys.path.append('./scripts/3.GeneticScores/pytorch_files/')
from GNNPatientSAGE import PatientGNNSAGE
from GNNPatientSAGE_v2 import PatientGNNSAGE_v2
from GNNPatientSAGE_v3 import PatientGNNSAGE_v3
from GNNPatientSAGE_v4 import PatientGNNSAGE_v4
from GNNPatientSAGEFocal import PatientGNNSAGEFocal

import json

train_graphs, val_graphs = torch.load("results/gnn/preprocess/graphs_genesEncScGPT.pt")
# Create DataLoader instances for both training and test sets
train_loader = DataLoader(train_graphs, batch_size=32, shuffle=True)
test_loader = DataLoader(val_graphs, batch_size=64, shuffle=False)



# train model
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

# Hyperparameters
patient_feat_dim = 13
gene_feat_dim = 16
hidden_gene_dim = 32
hidden_dim = 16
out_dim = 1  # Single output for survival regression
learning_rate = 0.001
hidden_vaf = 8

gnn_model_bool = PatientGNNSAGE(patient_feat_dim, gene_feat_dim, hidden_gene_dim, hidden_dim, out_dim, learning_rate)
gnn_model_vaf = PatientGNNSAGE(patient_feat_dim, gene_feat_dim, hidden_gene_dim, hidden_dim, out_dim, use_vaf=True, hidden_vaf = hidden_vaf, learning_rate = learning_rate)


# train model
logger_gnn_bool = TensorBoardLogger("lightning_logs/current", name="gnn_boolean")
logger_gnn_vaf = TensorBoardLogger("lightning_logs/current", name="gnn_vaf")

trainer_bool = L.Trainer(max_epochs = 30, logger = logger_gnn_bool, callbacks=[define_checkpoint_callback("GNN_boolean")])
trainer_vaf = L.Trainer(max_epochs = 30, logger = logger_gnn_vaf, callbacks=[define_checkpoint_callback("GNN_vaf")])

trainer_bool.fit(model=gnn_model_bool, train_dataloaders=train_loader, val_dataloaders=test_loader)
trainer_vaf.fit(model=gnn_model_vaf, train_dataloaders=train_loader, val_dataloaders=test_loader)

## Train again with optimized hyperparameters
### Boolean
hidden_gene_dim = 10
hidden_dim = 60
out_dim = 1  # Single output for survival regression
learning_rate = 0.0001

gnn_model_bool_optim = PatientGNNSAGE(patient_feat_dim, gene_feat_dim, hidden_gene_dim, hidden_dim, out_dim, learning_rate)
# train model
logger_gnn_bool_optim = TensorBoardLogger("lightning_logs/current", name="gnn_boolean_optim")
trainer_bool_optim = L.Trainer(max_epochs = 30, logger = logger_gnn_bool_optim, callbacks=[define_checkpoint_callback("GNN_boolean_optim")])
trainer_bool_optim.fit(model=gnn_model_bool_optim, train_dataloaders=train_loader, val_dataloaders=test_loader)

##
hidden_gene_dim = 15
hidden_dim = 70
out_dim = 1  # Single output for survival regression
learning_rate = 0.0001
hidden_vaf = 4

gnn_model_vaf_optim = PatientGNNSAGE(patient_feat_dim, gene_feat_dim, hidden_gene_dim, hidden_dim, out_dim, use_vaf=True, hidden_vaf = hidden_vaf, learning_rate = learning_rate)
logger_gnn_vaf_optim = TensorBoardLogger("lightning_logs/current", name="gnn_vaf_optim")
trainer_vaf_optim = L.Trainer(max_epochs = 30, logger = logger_gnn_vaf_optim, callbacks=[define_checkpoint_callback("GNN_vaf_optim")])
trainer_vaf_optim.fit(model=gnn_model_vaf_optim, train_dataloaders=train_loader, val_dataloaders=test_loader)

## Version 2: add gene update
gnn_model_v2_bool = PatientGNNSAGE_v2(patient_feat_dim, gene_feat_dim, hidden_gene_dim, hidden_dim, out_dim, learning_rate)
gnn_model_v2_vaf = PatientGNNSAGE_v2(patient_feat_dim, gene_feat_dim, hidden_gene_dim, hidden_dim, out_dim, use_vaf=True, hidden_vaf = hidden_vaf, learning_rate = learning_rate)


# train model
logger_gnn_bool = TensorBoardLogger("lightning_logs/current", name="gnn_boolean_v2")
logger_gnn_vaf = TensorBoardLogger("lightning_logs/current", name="gnn_vaf_v2")

trainer_bool_v2 = L.Trainer(max_epochs = 30, logger = logger_gnn_bool, callbacks=[define_checkpoint_callback("GNN_boolean_v2")])
trainer_vaf_v2 = L.Trainer(max_epochs = 30, logger = logger_gnn_vaf, callbacks=[define_checkpoint_callback("GNN_vaf_v2")])

trainer_bool_v2.fit(model=gnn_model_v2_bool, train_dataloaders=train_loader, val_dataloaders=test_loader)
trainer_vaf_v2.fit(model=gnn_model_v2_vaf, train_dataloaders=train_loader, val_dataloaders=test_loader)

## Version 3: add NN after second SAGE layer
gnn_model_v3_bool = PatientGNNSAGE_v3(patient_feat_dim, gene_feat_dim, hidden_gene_dim, hidden_dim, out_dim, learning_rate)
gnn_model_v3_vaf = PatientGNNSAGE_v3(patient_feat_dim, gene_feat_dim, hidden_gene_dim, hidden_dim, out_dim, use_vaf=True, hidden_vaf = hidden_vaf, learning_rate = learning_rate)

logger_gnn_bool_v3 = TensorBoardLogger("lightning_logs/current", name="gnn_boolean_v3")
logger_gnn_vaf_v3 = TensorBoardLogger("lightning_logs/current", name="gnn_vaf_v3")

trainer_bool_v3 = L.Trainer(max_epochs = 30, logger = logger_gnn_bool_v3, callbacks=[define_checkpoint_callback("GNN_boolean_v3")])
trainer_vaf_v3 = L.Trainer(max_epochs = 30, logger = logger_gnn_vaf_v3, callbacks=[define_checkpoint_callback("GNN_vaf_v3")])

trainer_bool_v3.fit(model=gnn_model_v3_bool, train_dataloaders=train_loader, val_dataloaders=test_loader)
trainer_vaf_v3.fit(model=gnn_model_v3_vaf, train_dataloaders=train_loader, val_dataloaders=test_loader)

## Version 4: add Dropout, BatchNorm, and L2 Regularization
# Hiperparámetros para v4 con regularización
dropout_v4 = 0.3        # 30% dropout
l2_reg_v4 = 1e-4        # L2 regularization
hidden_dim_v4 = 64      # Dimensión más grande para aprovechar la regularización

# Boolean model v4
gnn_model_v4_bool = PatientGNNSAGE_v4(
    patient_feat_dim=patient_feat_dim, 
    gene_feat_dim=gene_feat_dim, 
    hidden_gene_dim=hidden_gene_dim, 
    hidden_dim=hidden_dim_v4, 
    out_dim=out_dim, 
    dropout=dropout_v4,
    l2_reg=l2_reg_v4,
    learning_rate=learning_rate,
    use_vaf=False
)

# VAF model v4
gnn_model_v4_vaf = PatientGNNSAGE_v4(
    patient_feat_dim=patient_feat_dim, 
    gene_feat_dim=gene_feat_dim, 
    hidden_gene_dim=hidden_gene_dim, 
    hidden_dim=hidden_dim_v4, 
    out_dim=out_dim, 
    dropout=dropout_v4,
    l2_reg=l2_reg_v4,
    learning_rate=learning_rate,
    use_vaf=True,
    hidden_vaf=hidden_vaf
)

# Loggers para v4
logger_gnn_bool_v4 = TensorBoardLogger("lightning_logs/current", name="gnn_boolean_v4")
logger_gnn_vaf_v4 = TensorBoardLogger("lightning_logs/current", name="gnn_vaf_v4")

# Trainers para v4 (más épocas para aprovechar la regularización)
trainer_bool_v4 = L.Trainer(
    max_epochs=100, 
    logger=logger_gnn_bool_v4, 
    callbacks=[define_checkpoint_callback("GNN_boolean_v4")]
)
trainer_vaf_v4 = L.Trainer(
    max_epochs=100, 
    logger=logger_gnn_vaf_v4, 
    callbacks=[define_checkpoint_callback("GNN_vaf_v4")]
)

# Entrenar modelos v4
trainer_bool_v4.fit(model=gnn_model_v4_bool, train_dataloaders=train_loader, val_dataloaders=test_loader)
trainer_vaf_v4.fit(model=gnn_model_v4_vaf, train_dataloaders=train_loader, val_dataloaders=test_loader)




## Version 4 con hiperparámetros optimizados
# Hiperparámetros optimizados para v4
hidden_gene_dim_optim = 44
hidden_dim_optim = 55
learning_rate_optim = 7e-5
dropout_optim = 0.1     
l2_reg_optim = 4.4e-5     # L2 un poco menor

# Boolean model v4 optimizado
gnn_model_v4_bool_optim = PatientGNNSAGE_v4(
    patient_feat_dim=patient_feat_dim, 
    gene_feat_dim=gene_feat_dim, 
    hidden_gene_dim=hidden_gene_dim_optim, 
    hidden_dim=hidden_dim_optim, 
    out_dim=out_dim, 
    dropout=dropout_optim,
    l2_reg=l2_reg_optim,
    learning_rate=learning_rate_optim,
    use_vaf=False
)

hidden_gene_dim_optim_vaf = 26
hidden_dim_optim_vaf = 57
learning_rate_optim_vaf = 6e-04
dropout_optim_vaf = 0.36
l2_reg_optim_vaf = 8.5e-04
hidden_vaf_optim = 4


# VAF model v4 optimizado
gnn_model_v4_vaf_optim = PatientGNNSAGE_v4(
    patient_feat_dim=patient_feat_dim,
    gene_feat_dim=gene_feat_dim,
    hidden_gene_dim=hidden_gene_dim_optim_vaf,
    hidden_dim=hidden_dim_optim_vaf,
    out_dim=out_dim,
    dropout=dropout_optim_vaf,
    l2_reg=l2_reg_optim_vaf,
    learning_rate=learning_rate_optim_vaf,
    use_vaf=True,
    hidden_vaf=hidden_vaf_optim
)

# Loggers para v4 optimizados
logger_gnn_bool_v4_optim = TensorBoardLogger("lightning_logs/current", name="gnn_boolean_v4_optim")
logger_gnn_vaf_v4_optim = TensorBoardLogger("lightning_logs/current", name="gnn_vaf_v4_optim")

# Trainers para v4 optimizados
trainer_bool_v4_optim = L.Trainer(
    max_epochs=150, 
    logger=logger_gnn_bool_v4_optim, 
    callbacks=[define_checkpoint_callback("GNN_boolean_v4_optim")]
)
trainer_vaf_v4_optim = L.Trainer(
    max_epochs=150, 
    logger=logger_gnn_vaf_v4_optim, 
    callbacks=[define_checkpoint_callback("GNN_vaf_v4_optim")]
)

# Entrenar modelos v4 optimizados
trainer_bool_v4_optim.fit(model=gnn_model_v4_bool_optim, train_dataloaders=train_loader, val_dataloaders=test_loader)
trainer_vaf_v4_optim.fit(model=gnn_model_v4_vaf_optim, train_dataloaders=train_loader, val_dataloaders=test_loader)


## Version Focal: Focal Cox Loss
# Hiperparámetros para Focal Loss
focal_alpha = 0.25
focal_gamma = 2.0

# Boolean model con Focal Loss
gnn_model_focal_bool = PatientGNNSAGEFocal(
    patient_feat_dim=patient_feat_dim,
    gene_feat_dim=gene_feat_dim,
    hidden_gene_dim=hidden_gene_dim,
    hidden_dim=hidden_dim,
    out_dim=out_dim,
    learning_rate=learning_rate,
    use_vaf=False,
    focal_alpha=focal_alpha,
    focal_gamma=focal_gamma
)

# VAF model con Focal Loss
gnn_model_focal_vaf = PatientGNNSAGEFocal(
    patient_feat_dim=patient_feat_dim,
    gene_feat_dim=gene_feat_dim,
    hidden_gene_dim=hidden_gene_dim,
    hidden_dim=hidden_dim,
    out_dim=out_dim,
    learning_rate=learning_rate,
    use_vaf=True,
    hidden_vaf=hidden_vaf,
    focal_alpha=focal_alpha,
    focal_gamma=focal_gamma
)

# Loggers para Focal
logger_focal_bool = TensorBoardLogger("lightning_logs/current", name="gnn_focal_boolean")
logger_focal_vaf = TensorBoardLogger("lightning_logs/current", name="gnn_focal_vaf")

# Trainers para Focal
trainer_focal_bool = L.Trainer(max_epochs=30, logger=logger_focal_bool, callbacks=[define_checkpoint_callback("GNN_focal_boolean")])
trainer_focal_vaf = L.Trainer(max_epochs=30, logger=logger_focal_vaf, callbacks=[define_checkpoint_callback("GNN_focal_vaf")])

# Entrenar modelos con Focal Loss
trainer_focal_bool.fit(model=gnn_model_focal_bool, train_dataloaders=train_loader, val_dataloaders=test_loader)
trainer_focal_vaf.fit(model=gnn_model_focal_vaf, train_dataloaders=train_loader, val_dataloaders=test_loader)

