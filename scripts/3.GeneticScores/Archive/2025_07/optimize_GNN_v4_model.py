"""
Optimize GNN v4 model (with Dropout, BatchNorm, L2) using Optuna.
docker run --gpus all -v $PWD:$PWD -w $PWD -it mds_subtypes_python:1.5 python
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
from GNNPatientSAGE_v4 import PatientGNNSAGE_v4

class GNNOptimizer:
    def __init__(self, batch_size=32, k_folds=5, max_epochs=20):
        self.batch_size = batch_size
        self.k_folds = k_folds
        self.max_epochs = max_epochs
        self.train_graphs, _ = torch.load("results/gnn/preprocess/graphs_genesEncScGPT.pt")
        
    def run_cross_validation(self, trial, model_params, use_vaf=False, log_dir="lightning_logs"):
        """
        Ejecuta validación cruzada para un conjunto de hiperparámetros
        """
        kfold = KFold(n_splits=self.k_folds, shuffle=True, random_state=42)
        fold_scores = []
        
        for fold, (train_idx, val_idx) in enumerate(kfold.split(range(len(self.train_graphs)))):
            # Crear datasets para este fold
            train_dataset = torch.utils.data.Subset(self.train_graphs, train_idx)
            val_dataset = torch.utils.data.Subset(self.train_graphs, val_idx)      
            
            # Crear dataloaders
            train_loader = DataLoader(train_dataset, batch_size=self.batch_size, shuffle=True)
            val_loader = DataLoader(val_dataset, batch_size=len(val_dataset), shuffle=False)
            
            # Inicializar modelo v4
            model = PatientGNNSAGE_v4(**model_params)
            
            # Configurar trainer
            logger_name = self._create_logger_name(fold, model_params, use_vaf)
            logger = TensorBoardLogger(log_dir, name=logger_name)
            
            trainer = L.Trainer(
                logger=logger,
                max_epochs=self.max_epochs,
                enable_progress_bar=True,
                callbacks=[PyTorchLightningPruningCallback(trial, monitor="val_c_index")]
            )   
            
            # Entrenar modelo
            trainer.fit(model, train_loader, val_loader)
            
            # Usar c_index (mayor es mejor, así que negativo para minimizar)
            c_index = trainer.callback_metrics.get("val_c_index", 0.5)
            fold_scores.append(-c_index.item())  # Negativo porque Optuna minimiza
            
        # Registrar scores individuales de cada fold
        for fold, score in enumerate(fold_scores):
            trial.set_user_attr(f'fold_{fold}_score', score)
        
        return np.mean(fold_scores)
    
    def _create_logger_name(self, fold, params, use_vaf):
        """Crear nombre descriptivo para el logger"""
        base_name = f"optuna_trial_fold_{fold}_gd_{params['hidden_gene_dim']}_hd_{params['hidden_dim']}_dr_{params['dropout']:.2f}_l2_{params['l2_reg']:.0e}"
        
        if use_vaf:
            base_name += f"_vd_{params['hidden_vaf']}"
            
        return base_name
    
    def get_base_model_params(self, trial_params):
        """Parámetros base comunes a ambos modelos"""
        return {
            'patient_feat_dim': 13,
            'gene_feat_dim': 16,
            'hidden_gene_dim': trial_params['hidden_gene_dim'],
            'hidden_dim': trial_params['hidden_dim'],
            'out_dim': 1,
            'dropout': trial_params['dropout'],
            'l2_reg': trial_params['l2_reg'],
            'learning_rate': trial_params['learning_rate']
        }

def objective_v4_boolean(trial):
    """Función objetivo para modelo Boolean v4"""
    # Definir el espacio de búsqueda (sin batch_size)
    trial_params = {
        'learning_rate': trial.suggest_float('learning_rate', 1e-5, 1e-2, log=True),
        'hidden_dim': trial.suggest_int('hidden_dim', 16, 128),
        'hidden_gene_dim': trial.suggest_int('hidden_gene_dim', 8, 64),
        'dropout': trial.suggest_float('dropout', 0.1, 0.5),
        'l2_reg': trial.suggest_float('l2_reg', 1e-6, 1e-3, log=True)
    }
    
    # Usar el optimizador común
    optimizer = GNNOptimizer()
    
    # Parámetros del modelo Boolean
    model_params = optimizer.get_base_model_params(trial_params)
    model_params['use_vaf'] = False
    
    # Ejecutar validación cruzada
    return optimizer.run_cross_validation(
        trial, 
        model_params, 
        use_vaf=False, 
        log_dir="lightning_logs/gnn_v4_boolean/"
    )

def objective_v4_vaf(trial):
    """Función objetivo para modelo VAF v4"""
    # Definir el espacio de búsqueda (sin batch_size, pero con hidden_vaf)
    trial_params = {
        'learning_rate': trial.suggest_float('learning_rate', 1e-5, 1e-2, log=True),
        'hidden_dim': trial.suggest_int('hidden_dim', 16, 128),
        'hidden_gene_dim': trial.suggest_int('hidden_gene_dim', 8, 64),
        'hidden_vaf': trial.suggest_categorical('hidden_vaf', [4, 6, 8, 12, 16]),
        'dropout': trial.suggest_float('dropout', 0.1, 0.5),
        'l2_reg': trial.suggest_float('l2_reg', 1e-6, 1e-3, log=True)
    }
    
    # Usar el optimizador común
    optimizer = GNNOptimizer()
    
    # Parámetros del modelo VAF
    model_params = optimizer.get_base_model_params(trial_params)
    model_params['use_vaf'] = True
    model_params['hidden_vaf'] = trial_params['hidden_vaf']
    
    # Ejecutar validación cruzada
    return optimizer.run_cross_validation(
        trial, 
        model_params, 
        use_vaf=True, 
        log_dir="lightning_logs/gnn_v4_vaf/"
    )

def save_results(study, model_type, output_dir):
    """Guardar resultados de la optimización"""
    results = {
        "model_type": model_type,
        "best_params": study.best_params,
        "best_c_index": -study.best_value,
        "best_value": study.best_value,
        "all_trials": [
            {
                "params": trial.params,
                "c_index": -trial.value,
                "value": trial.value,
                "fold_scores": {
                    f"fold_{i}": trial.user_attrs.get(f'fold_{i}_score', None)
                    for i in range(5)
                }
            }
            for trial in study.trials
        ]
    }
    
    with open(f"{output_dir}/optuna_results_gnn_v4_{model_type.lower()}.json", "w") as f:
        json.dump(results, f, indent=4)
    
    # Crear visualizaciones
    try:
        fig = optuna.visualization.plot_param_importances(study)
        fig.write_html(f"{output_dir}/param_importance.html")
        
        fig = optuna.visualization.plot_optimization_history(study)
        fig.write_html(f"{output_dir}/optimization_history.html")
        
        fig = optuna.visualization.plot_parallel_coordinate(study)
        fig.write_html(f"{output_dir}/parallel_coordinate.html")
    except Exception as e:
        print(f"Error creando visualizaciones: {e}")

# ===============================
# OPTIMIZACIÓN MODELO BOOLEAN V4
# ===============================

print("=== Iniciando optimización GNN v4 Boolean ===")

study_v4_bool = optuna.create_study(
    direction="minimize",  # Minimizamos porque usamos -c_index
    study_name="gnn_v4_optimization_boolean"
)
    
study_v4_bool.optimize(objective_v4_boolean, n_trials=100) 
    
print("Número de trials completados (Boolean):", len(study_v4_bool.trials))
print("Mejores hiperparámetros (Boolean):", study_v4_bool.best_params)
print("Mejor C-index (Boolean):", -study_v4_bool.best_value)

# Guardar resultados Boolean
save_results(study_v4_bool, "Boolean", "results/gnn/gnn_optimization_v4_boolean")

print("=== Optimización Boolean completada ===\n")

# ===============================
# OPTIMIZACIÓN MODELO VAF V4
# ===============================

print("=== Iniciando optimización GNN v4 VAF ===")

study_v4_vaf = optuna.create_study(
    direction="minimize",  # Minimizamos porque usamos -c_index
    study_name="gnn_v4_optimization_vaf"
)
    
study_v4_vaf.optimize(objective_v4_vaf, n_trials=100) 
    
print("Número de trials completados (VAF):", len(study_v4_vaf.trials))
print("Mejores hiperparámetros (VAF):", study_v4_vaf.best_params)
print("Mejor C-index (VAF):", -study_v4_vaf.best_value)

# Guardar resultados VAF
save_results(study_v4_vaf, "VAF", "results/gnn/gnn_optimization_v4_vaf")

print("=== Optimización VAF completada ===")

# ===============================
# RESUMEN FINAL
# ===============================

print("\n" + "="*50)
print("RESUMEN DE OPTIMIZACIÓN GNN v4")
print("="*50)
print(f"Batch size fijo: 32")
print(f"Boolean - Mejor C-index: {-study_v4_bool.best_value:.4f}")
print(f"Boolean - Mejores params: {study_v4_bool.best_params}")
print(f"\nVAF - Mejor C-index: {-study_v4_vaf.best_value:.4f}")
print(f"VAF - Mejores params: {study_v4_vaf.best_params}")
print("="*50)