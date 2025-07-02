"""
Optimize GNN model with reduced gene embeddings using Optuna.
docker run --gpus all -v $PWD:$PWD -w $PWD -it mds_subtypes_python:1.4 python
"""

import optuna
import sys
from optuna.integration import PyTorchLightningPruningCallback
import lightning as L
from lightning.pytorch.loggers import TensorBoardLogger
from torch_geometric.loader import DataLoader
import torch
from torch.utils.data import random_split
import json
from sklearn.model_selection import KFold
import numpy as np
sys.path.append('./scripts/3.GeneticScores/pytorch_files/')
from GNNPatientSAGE import PatientGNNSAGE

def objective(trial):
    # Definir el espacio de búsqueda
    params = {
        'learning_rate': trial.suggest_float('learning_rate', 1e-4, 1e-1, log=True),
        'hidden_dim': trial.suggest_int('hidden_dim', 8, 64),
        'hidden_gene_dim': trial.suggest_int('hidden_gene_dim', 8, 64),
        'batch_size': trial.suggest_categorical('batch_size', [16, 32, 64])
    }  
    # Cargar datos
    patient_graphs = torch.load("results/gnn/preprocess/boolean_scGPT_reduced.pt")  
    # Split datos
    # Configurar k-fold cross validation
    k_folds = 5
    kfold = KFold(n_splits=k_folds, shuffle=True, random_state=42)
    # Lista para almacenar los resultados de cada fold
    fold_scores = []
    # Realizar k-fold cross validation
    for fold, (train_idx, val_idx) in enumerate(kfold.split(range(len(patient_graphs)))):
        # Crear datasets para este fold
        train_dataset = torch.utils.data.Subset(patient_graphs, train_idx)
        val_dataset = torch.utils.data.Subset(patient_graphs, val_idx)      
        # Crear dataloaders
        train_loader = DataLoader(train_dataset, batch_size=params['batch_size'], shuffle=True)
        val_loader = DataLoader(val_dataset, batch_size=len(val_dataset), shuffle=False)
        # Inicializar modelo
        model = PatientGNNSAGE(
            patient_feat_dim=15,
            gene_feat_dim=16,
            hidden_gene_dim=params['hidden_gene_dim'],
            hidden_dim=params['hidden_dim'],
            out_dim=1,
            learning_rate=params['learning_rate']
        )    
        # Configurar trainer
        logger = TensorBoardLogger(
            "lightning_logs/gnn_reduced_2/",
            name=f"optuna_trial_fold_{fold}_gd_{params['hidden_gene_dim']}_hd_{params['hidden_dim']}_bs_{params['batch_size']}"
        )
        trainer = L.Trainer(
            logger=logger,
            max_epochs=15,
            enable_progress_bar=True
        )   
        # Entrenar modelo
        trainer.fit(model, train_loader, val_loader)
        fold_scores.append(trainer.callback_metrics["val_c_index"].item())
    # Retornar la media de los scores
    mean_score = np.mean(fold_scores)
    
    # Registrar scores individuales de cada fold
    for fold, score in enumerate(fold_scores):
        trial.set_user_attr(f'fold_{fold}_score', score)
    return mean_score

study = optuna.create_study(
    direction="maximize", 
    study_name="gnn_optimization_2"
)
    
study.optimize(objective, n_trials=50)  # Ajustar número de trials según necesites
    
print("Número de trials completados:", len(study.trials))
print("Mejores hiperparámetros:", study.best_params)
print("Mejor valor:", study.best_value)
    
# Guardar resultados
results = {
    "best_params": study.best_params,
    "best_value": study.best_value,
    "all_trials": [
        {
        "params": trial.params,
        "value": trial.value
        }
    for trial in study.trials
]
}
    
with open("results/gnn/gnn_reduc_optimization/optuna_results_gnn_reduc2.json", "w") as f:
        json.dump(results, f, indent=4)



# Graficar importancia de parámetros
fig = optuna.visualization.plot_param_importances(study)
fig.write_html("results/gnn/gnn_reduc_optimization/param_importance2.html")
    
# Graficar historia de optimización
fig = optuna.visualization.plot_optimization_history(study)
fig.write_html("results/gnn/gnn_reduc_optimization/optimization_history2.html")

# Graficar relaciones entre parámetros
fig = optuna.visualization.plot_parallel_coordinate(study)
fig.write_html("results/gnn/gnn_reduc_optimization/parallel_coordinate2.html")

