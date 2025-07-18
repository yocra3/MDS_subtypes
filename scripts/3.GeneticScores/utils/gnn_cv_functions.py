"""
GNN Cross-Validation Helper Functions
====================================

Funciones auxiliares para realizar cross-validation con modelos GNN.
Reutiliza el código de entrenamiento existente de trainGNNmodel.py y lo adapta 
para trabajar con estructura de folds de CV.

Estructura de datos esperada:
- graph_data: Dict con claves 'fold_0', 'fold_1', ..., 'fold_9'
- Cada fold contiene una lista de objetos HeteroData
- Se deja un fold para test y el resto para train

Autores: Adaptado del pipeline MDS_subtypes
Fecha: 2024
"""

import os
import sys
import logging
import pickle
import yaml
from typing import Dict, List, Tuple, Any, Optional
from pathlib import Path

import torch
import numpy as np
import pandas as pd
import lightning as L
from torch_geometric.loader import DataLoader
from lightning.pytorch.loggers import TensorBoardLogger
from lightning.pytorch.callbacks import ModelCheckpoint

# Añadir path para imports de modelos
sys.path.append('./scripts/3.GeneticScores/pytorch_files/')
from GNNPatientSAGE import PatientGNNSAGE
from GNNPatientSAGE_v2 import PatientGNNSAGE_v2
from GNNPatientSAGE_v3 import PatientGNNSAGE_v3
from GNNPatientSAGE_v4 import PatientGNNSAGE_v4
from GNNPatientSAGEFocal import PatientGNNSAGEFocal

# Configurar logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class GNNCrossValidator:
    """
    Clase principal para realizar cross-validation con modelos GNN.
    """
    
    def __init__(self, 
                 graph_data_path: str,
                 config_path: str = "scripts/3.GeneticScores/configs/gnn_models_config.yaml",
                 results_dir: str = "results/gnn/cv_results",
                 lightning_logs_dir: str = "lightning_logs/gnn_cv"):
        """
        Inicializa el validador cruzado.
        
        Args:
            graph_data_path: Ruta al archivo .pt con los datos de grafos por fold
            config_path: Ruta al archivo de configuración YAML
            results_dir: Directorio para guardar resultados de CV
            lightning_logs_dir: Directorio para logs de Lightning
        """
        self.graph_data_path = graph_data_path
        self.config_path = config_path
        self.results_dir = Path(results_dir)
        self.lightning_logs_dir = Path(lightning_logs_dir)
        
        # Crear directorios si no existen
        self.results_dir.mkdir(parents=True, exist_ok=True)
        self.lightning_logs_dir.mkdir(parents=True, exist_ok=True)
        
        # Almacenar datos de grafos y configuración
        self.graph_data = None
        self.n_folds = None
        self.config = None
        self.patient_feat_dim = None
        self.gene_feat_dim = None
        
        # Cargar configuración
        self.load_config()
        
    def load_config(self):
        """Carga la configuración desde el archivo YAML."""
        if not os.path.exists(self.config_path):
            logger.warning(f"Archivo de configuración no encontrado: {self.config_path}")
            logger.info("Usando configuración por defecto")
            self.config = self._get_default_config()
        else:
            with open(self.config_path, 'r') as f:
                self.config = yaml.safe_load(f)
            logger.info(f"Configuración cargada desde: {self.config_path}")
    
    def _get_default_config(self) -> Dict[str, Any]:
        """Devuelve configuración por defecto si no se encuentra archivo YAML."""
        return {
            'models': {
                'v1_boolean': {
                    'model_type': 'v1',
                    'use_vaf': False,
                    'hyperparameters': {
                        'hidden_gene_dim': 32,
                        'hidden_dim': 16,
                        'out_dim': 1,
                        'learning_rate': 0.001,
                        'hidden_vaf': 8
                    },
                    'training': {
                        'max_epochs': 30,
                        'batch_size_train': 32,
                        'batch_size_test': 64
                    }
                }
            },
            'cross_validation': {
                'shuffle_train': True,
                'shuffle_test': False
            }
        }
        
    def load_graph_folds(self) -> Dict[str, List]:
        """
        Carga todos los folds de datos de grafos desde archivo.
        También extrae las dimensiones de entrada automáticamente.
        
        Returns:
            Dict con datos de grafos por fold
        """
        logger.info(f"Cargando datos de grafos desde: {self.graph_data_path}")
        
        if not os.path.exists(self.graph_data_path):
            raise FileNotFoundError(f"Archivo de grafos no encontrado: {self.graph_data_path}")
        
        self.graph_data = torch.load(self.graph_data_path)
        self.n_folds = len(self.graph_data)
        
        # Extraer dimensiones de entrada del primer grafo disponible
        first_fold_key = list(self.graph_data.keys())[0]
        first_graph = self.graph_data[first_fold_key][0]
        
        self.patient_feat_dim = first_graph['patient'].x.shape[1]
        self.gene_feat_dim = first_graph['gene'].x.shape[1]
        
        logger.info(f"Cargados {self.n_folds} folds")
        logger.info(f"Dimensiones detectadas - Patient: {self.patient_feat_dim}, Gene: {self.gene_feat_dim}")
        
        for fold_name, graphs in self.graph_data.items():
            logger.info(f"{fold_name}: {len(graphs)} grafos")
            
        return self.graph_data
    
    def create_fold_dataloaders(self, 
                               test_fold: int,
                               batch_size_train: int = 32,
                               batch_size_test: int = 64,
                               shuffle_train: bool = True) -> Tuple[DataLoader, DataLoader]:
        """
        Crea DataLoaders para entrenamiento y test basado en el fold de test especificado.
        
        Args:
            test_fold: Número del fold a usar como test (0-based)
            batch_size_train: Tamaño de batch para entrenamiento
            batch_size_test: Tamaño de batch para test
            shuffle_train: Si hacer shuffle de datos de entrenamiento
            
        Returns:
            Tuple de (train_loader, test_loader)
        """
        if self.graph_data is None:
            raise ValueError("Debe cargar los datos de grafos primero con load_graph_folds()")
        
        if test_fold >= self.n_folds:
            raise ValueError(f"test_fold ({test_fold}) debe ser menor que n_folds ({self.n_folds})")
        
        # Obtener datos de test
        test_fold_key = f'fold_{test_fold}'
        test_graphs = self.graph_data[test_fold_key]
        
        # Obtener datos de entrenamiento (todos los otros folds)
        train_graphs = []
        for fold_idx in range(self.n_folds):
            if fold_idx != test_fold:
                fold_key = f'fold_{fold_idx}'
                train_graphs.extend(self.graph_data[fold_key])
        
        # Crear DataLoaders
        train_loader = DataLoader(train_graphs, 
                                 batch_size=batch_size_train, 
                                 shuffle=shuffle_train)
        test_loader = DataLoader(test_graphs, 
                                batch_size=batch_size_test, 
                                shuffle=False)
        
        logger.info(f"Fold {test_fold} como test: {len(test_graphs)} grafos")
        logger.info(f"Entrenamiento: {len(train_graphs)} grafos")
        
        return train_loader, test_loader
    
    def define_checkpoint_callback(self, model_name: str, fold_idx: int) -> ModelCheckpoint:
        """
        Define callback para guardar el mejor modelo por fold.
        
        Args:
            model_name: Nombre del modelo
            fold_idx: Índice del fold
            
        Returns:
            ModelCheckpoint callback
        """
        checkpoint_dir = self.results_dir / "model_checkpoints" / f"fold_{fold_idx}"
        checkpoint_dir.mkdir(parents=True, exist_ok=True)

        return ModelCheckpoint(
            dirpath=str(checkpoint_dir),
            filename=f"{model_name}_fold_{fold_idx}",
            save_top_k=self.config['checkpoint_config'].get('save_top_k', 1),
            monitor=self.config['checkpoint_config'].get('monitor', 'val_c_index'),
            mode=self.config['checkpoint_config'].get('mode', 'max')
        )
    
    def train_single_fold(self,
                         model: L.LightningModule,
                         train_loader: DataLoader,
                         test_loader: DataLoader,
                         fold_idx: int,
                         model_name: str,
                         max_epochs: int = 30) -> Dict[str, Any]:
        """
        Entrena un modelo en un solo fold.
        
        Args:
            model: Modelo Lightning a entrenar
            train_loader: DataLoader de entrenamiento
            test_loader: DataLoader de test
            fold_idx: Índice del fold
            model_name: Nombre del modelo para logging
            max_epochs: Número máximo de épocas
            
        Returns:
            Dict con métricas del fold
        """
        # Logger para este fold
        logger_name = f"{model_name}_fold_{fold_idx}"
        logger_tb = TensorBoardLogger(self.lightning_logs_dir, name=logger_name)
        
        # Checkpoint callback
        checkpoint_callback = self.define_checkpoint_callback(model_name, fold_idx)
        
        # Trainer
        trainer = L.Trainer(
            max_epochs=max_epochs,
            logger=logger_tb,
            callbacks=[checkpoint_callback],
            enable_progress_bar=True, 
            enable_model_summary=False
        )
        
        # Entrenar
        logger.info(f"Entrenando {model_name} - Fold {fold_idx}")
        trainer.fit(model=model, 
                   train_dataloaders=train_loader, 
                   val_dataloaders=test_loader)
        
        # Extraer métricas finales
        # Extraer el mejor c_index del checkpoint callback (si está disponible)
        best_c_index = None
        if hasattr(checkpoint_callback, 'best_model_score') and checkpoint_callback.best_model_score is not None:
            best_c_index = checkpoint_callback.best_model_score.item()
        
        metrics = {
            'fold': fold_idx,
            'final_val_loss': trainer.logged_metrics.get('val_loss', None).item(),
            'final_val_c_index': trainer.logged_metrics.get('val_c_index', None).item(),
            'best_val_c_index': best_c_index,
            'checkpoint_path': checkpoint_callback.best_model_path
        }
        
        logger.info(f"Fold {fold_idx} completado. C-index val: {metrics['final_val_c_index']}")
        
        return metrics
    
    def create_model_instance(self, 
                             model_name: str,
                             model_config: Optional[Dict[str, Any]] = None) -> L.LightningModule:
        """
        Crea una instancia del modelo especificado usando configuración YAML.
        
        Args:
            model_name: Nombre del modelo en el archivo de configuración
            model_config: Configuración manual (opcional, sobrescribe YAML)
            
        Returns:
            Instancia del modelo Lightning
        """
        if self.patient_feat_dim is None or self.gene_feat_dim is None:
            raise ValueError("Debe cargar los datos de grafos primero para obtener dimensiones")

        # Usar configuración del YAML o manual
        if model_config is None:
            if model_name not in self.config['models']:
                raise ValueError(f"Modelo {model_name} no encontrado en configuración")
            else:
                config = self.config['models'][model_name]
        else:
            # Convertir configuración manual al formato esperado
            config = {
                'model_type': model_config.get('model_type', 'v1'),
                'use_vaf': model_config.get('use_vaf', False),
                'hyperparameters': model_config,
                'training': {}
            }
        
        model_type = config['model_type']
        use_vaf = config['use_vaf']
        hyperparams = config['hyperparameters']
        
        # Parámetros base con dimensiones automáticas
        base_params = {
            'patient_feat_dim': self.patient_feat_dim,
            'gene_feat_dim': self.gene_feat_dim,
            'hidden_gene_dim': hyperparams.get('hidden_gene_dim', 32),
            'hidden_dim': hyperparams.get('hidden_dim', 16),
            'out_dim': hyperparams.get('out_dim', 1),
            'learning_rate': hyperparams.get('learning_rate', 0.001),
            'use_vaf': use_vaf,
            'hidden_vaf': hyperparams.get('hidden_vaf', 8)
        }
        if model_type == 'v1':
            return PatientGNNSAGE(**base_params)
        
        elif model_type == 'v2':
            return PatientGNNSAGE_v2(**base_params)
        
        elif model_type == 'v3':
            return PatientGNNSAGE_v3(**base_params)
        
        elif model_type == 'v4':
            v4_params = base_params.copy()
            v4_params.update({
                'dropout': hyperparams.get('dropout', 0.3),
                'l2_reg': hyperparams.get('l2_reg', 1e-4)
            })
            return PatientGNNSAGE_v4(**v4_params)
        
        elif model_type == 'focal':
            focal_params = base_params.copy()
            focal_params.update({
                'focal_alpha': hyperparams.get('focal_alpha', 0.25),
                'focal_gamma': hyperparams.get('focal_gamma', 2.0)
            })
            return PatientGNNSAGEFocal(**focal_params)
        
        else:
            raise ValueError(f"Tipo de modelo no soportado: {model_type}")
    
    def run_cross_validation(self,
                            model_name: str,
                            model_config: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
        """
        Ejecuta cross-validation completo para un modelo usando configuración YAML.
        
        Args:
            model_name: Nombre del modelo en el archivo de configuración
            model_config: Configuración manual (opcional, sobrescribe YAML)
            
        Returns:
            Dict con resultados de CV
        """
        
        if self.graph_data is None:
            self.load_graph_folds()
        
        # Obtener configuración del modelo
        if model_config is None:
            if model_name not in self.config['models']:
                raise ValueError(f"Modelo {model_name} no encontrado en configuración")
            config = self.config['models'][model_name]
            training_config = config['training']
        else:
            config = {
                'model_type': model_config.get('model_type', 'v1'),
                'use_vaf': model_config.get('use_vaf', False),
                'hyperparameters': model_config
            }
            training_config = {
                'max_epochs': model_config.get('max_epochs', 30),
                'batch_size_train': model_config.get('batch_size_train', 32),
                'batch_size_test': model_config.get('batch_size_test', 64)
            }
        
        max_epochs = training_config.get('max_epochs', 30)
        batch_size_train = training_config.get('batch_size_train', 32)
        batch_size_test = training_config.get('batch_size_test', 64)
        
        logger.info(f"Iniciando CV para modelo {model_name}")
        logger.info(f"Configuración: {max_epochs} épocas, batch sizes: {batch_size_train}/{batch_size_test}")
        
        fold_results = []
        
        # Iterar sobre cada fold como test
        for test_fold in range(self.n_folds):
            logger.info(f"Procesando fold {test_fold}/{self.n_folds-1}")
            
            # Crear DataLoaders para este fold
            train_loader, test_loader = self.create_fold_dataloaders(
                test_fold=test_fold,
                batch_size_train=batch_size_train,
                batch_size_test=batch_size_test
            )
            
            # Crear instancia del modelo (nueva para cada fold)
            model = self.create_model_instance(model_name, model_config)
            
            # Entrenar en este fold
            fold_metrics = self.train_single_fold(
                model=model,
                train_loader=train_loader,
                test_loader=test_loader,
                fold_idx=test_fold,
                model_name=model_name,
                max_epochs=max_epochs
            )
            
            fold_results.append(fold_metrics)
        
        # Calcular estadísticas agregadas
        val_c_indices = [r['final_val_c_index'] for r in fold_results if r['final_val_c_index'] is not None]
        val_losses = [r['final_val_loss'] for r in fold_results if r['final_val_loss'] is not None]
        
        cv_results = {
            'model_name': model_name,
            'model_config': config,
            'fold_results': fold_results,
            'n_folds': self.n_folds,
            
            # Dimensiones detectadas
            'data_dimensions': {
                'patient_feat_dim': self.patient_feat_dim,
                'gene_feat_dim': self.gene_feat_dim
            },
            
            # Estadísticas agregadas
            'mean_val_c_index': np.mean(val_c_indices) if val_c_indices else None,
            'std_val_c_index': np.std(val_c_indices) if val_c_indices else None,
            'mean_val_loss': np.mean(val_losses) if val_losses else None,
            'std_val_loss': np.std(val_losses) if val_losses else None,
            
            # Información adicional
            'training_config': training_config
        }
        
        # Guardar resultados
        self.save_cv_results(cv_results, model_name)
        
        logger.info(f"CV completado para {model_name}")
        if cv_results['mean_val_c_index'] is not None:
            logger.info(f"C-index promedio: {cv_results['mean_val_c_index']:.4f} ± {cv_results['std_val_c_index']:.4f}")
        
        return cv_results
    
    def save_cv_results(self, cv_results: Dict[str, Any], model_name: str):
        """
        Guarda los resultados de CV en archivo.
        
        Args:
            cv_results: Resultados de CV
            model_name: Nombre del modelo
        """
        # Guardar como pickle
        results_file = self.results_dir / f"{model_name}_cv_results.pkl"
        with open(results_file, 'wb') as f:
            pickle.dump(cv_results, f)
        
        # Guardar resumen como CSV
        fold_df = pd.DataFrame(cv_results['fold_results'])
        summary_file = self.results_dir / f"{model_name}_cv_summary.csv"
        fold_df.to_csv(summary_file, index=False)
        
        logger.info(f"Resultados guardados en: {results_file}")
        logger.info(f"Resumen guardado en: {summary_file}")


# Funciones de conveniencia
def run_gnn_cv_pipeline(graph_data_path, config_path, results_dir, lightning_logs_dir, models_to_run=None):
 
    """
    Ejecuta el pipeline completo de CV para múltiples modelos GNN usando configuración YAML.
    
    Args:
        graph_data_path: Ruta al archivo de datos de grafos
        config_path: Ruta al archivo de configuración YAML
        results_dir: Directorio para resultados
        models_to_run: Lista de modelos a ejecutar (None = usar modelos por defecto)
        
    Returns:
        Dict con resultados para cada modelo
    """
    # Crear validador
    cv_validator = GNNCrossValidator(
        graph_data_path=graph_data_path,
        config_path=config_path,
        results_dir=results_dir,
        lightning_logs_dir=lightning_logs_dir
    )
    
    # Determinar modelos a ejecutar
    if models_to_run is None:
        # Usar algunos modelos por defecto si no se especifica
        models_to_run = list(cv_validator.config['models'].keys())   
    
    logger.info(f"Ejecutando CV para modelos: {models_to_run}")
    
    # Ejecutar CV para cada modelo
    all_results = {}

    for model_name in models_to_run:
        if model_name not in cv_validator.config['models']:
            logger.warning(f"Modelo {model_name} no encontrado en configuración, saltando")
            continue
        try:
            results = cv_validator.run_cross_validation(model_name)
            all_results[model_name] = results
            
        except Exception as e:
            logger.error(f"Error ejecutando CV para {model_name}: {str(e)}")
            continue
    
    # Guardar resumen comparativo
    save_comparative_summary(all_results, results_dir)
    
    return all_results


def save_comparative_summary(all_results: Dict[str, Dict[str, Any]], results_dir: str):
    """
    Guarda un resumen comparativo de todos los modelos.
    
    Args:
        all_results: Resultados de CV para todos los modelos
        results_dir: Directorio de resultados
    """
    summary_data = []
    
    for model_name, results in all_results.items():
        model_config = results.get('model_config', {})
        data_dims = results.get('data_dimensions', {})
        
        summary_data.append({
            'model': model_name,
            'model_type': model_config.get('model_type', 'unknown'),
            'use_vaf': model_config.get('use_vaf', False),
            'patient_feat_dim': data_dims.get('patient_feat_dim', 'N/A'),
            'gene_feat_dim': data_dims.get('gene_feat_dim', 'N/A'),
            'mean_c_index': results['mean_val_c_index'],
            'std_c_index': results['std_val_c_index'],
            'mean_loss': results['mean_val_loss'],
            'std_loss': results['std_val_loss'],
            'n_folds': results['n_folds']
        })
    
    summary_df = pd.DataFrame(summary_data)
    if len(summary_df) > 0 and 'mean_c_index' in summary_df.columns:
        # Ordenar por C-index solo si hay datos válidos
        valid_c_index = summary_df['mean_c_index'].notna()
        if valid_c_index.any():
            summary_df = summary_df.sort_values('mean_c_index', ascending=False)
    
    summary_file = Path(results_dir) / "comparative_cv_summary.csv"
    summary_df.to_csv(summary_file, index=False)
    
    logger.info(f"Resumen comparativo guardado en: {summary_file}")
    
    # Mostrar resumen en log
    logger.info("\n" + "="*80)
    logger.info("RESUMEN COMPARATIVO DE MODELOS GNN")
    logger.info("="*80)
    for _, row in summary_df.iterrows():
        c_idx = row['mean_c_index']
        c_std = row['std_c_index']
        if pd.notna(c_idx) and pd.notna(c_std):
            logger.info(f"{row['model']:20s} | C-index: {c_idx:.4f} ± {c_std:.4f} | Dims: {row['patient_feat_dim']}x{row['gene_feat_dim']}")
        else:
            logger.info(f"{row['model']:20s} | C-index: Error en cálculo")
    logger.info("="*80)


if __name__ == "__main__":
    """
    Ejemplo de uso del módulo.
    """
    # Ejemplo básico
    graph_data_path = "results/gnn/preprocess/graphs_genesEncScGPT_folds.pt"
    config_path = "scripts/3.GeneticScores/configs/gnn_models_config.yaml"
    
    # Ejecutar CV para modelos seleccionados
    results = run_gnn_cv_pipeline(
        graph_data_path=graph_data_path,
        config_path=config_path,
        models_to_run=['v1_boolean', 'v4_boolean_optim']  # Solo algunos modelos para prueba
    )
    
    print("Cross-validation completado!")
    for model_name, result in results.items():
        c_idx = result.get('mean_val_c_index', 'N/A')
        c_std = result.get('std_val_c_index', 'N/A')
        
        if c_idx != 'N/A' and c_std != 'N/A':
            print(f"{model_name}: C-index = {c_idx:.4f} ± {c_std:.4f}")
        else:
            print(f"{model_name}: Error en cálculo de C-index")
