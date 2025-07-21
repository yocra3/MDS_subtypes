"""
GNN gene function Helper Functions
====================================

Funciones auxiliares para entrenar modelos GNN a nivel de gene.

Estructura de datos esperada:
- graph_data: lista de objetos HeteroData, cada uno representando un paciente
- gene_dict: diccionario con qué pacientes tienes cada gen
- Se deja un gen para testear y se entrena con el resto de genes.

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

# Configurar logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class GNNGeneTrainer:
    """
    Clase principal para entrenar modelos GNN a nivel de gen.
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
        self.gene_dict = None
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

    def load_graph_genes(self) -> Tuple[List, Dict[str, List]]:
        """
        Carga todos los folds de datos de grafos desde archivo.
        También extrae las dimensiones de entrada automáticamente.
        
        Returns:
           Lista de grafos y diccionario de genes.
        """
        logger.info(f"Cargando datos de grafos desde: {self.graph_data_path}")
        
        if not os.path.exists(self.graph_data_path):
            raise FileNotFoundError(f"Archivo de grafos no encontrado: {self.graph_data_path}")

        self.graph_data, self.gene_dict = torch.load(self.graph_data_path)
        
        # Extraer dimensiones de entrada del primer grafo disponible
        first_graph = self.graph_data[0]
        
        self.patient_feat_dim = first_graph['patient'].x.shape[1]
        self.gene_feat_dim = first_graph['gene'].x.shape[1]
        
        logger.info(f"Dimensiones detectadas - Patient: {self.patient_feat_dim}, Gene: {self.gene_feat_dim}")

        return self.graph_data, self.gene_dict

       
    def create_gene_dataloaders(self, 
                               test_gene: str,
                               batch_size_train: int = 32,
                               batch_size_test: int = 64,
                               shuffle_train: bool = True) -> Tuple[DataLoader, DataLoader]:
        """
        Crea DataLoaders para entrenamiento y test basado en el gen especificado.
        
        Args:
            test_gene: Nombre del gen a usar como test
            batch_size_train: Tamaño de batch para entrenamiento
            batch_size_test: Tamaño de batch para test
            shuffle_train: Si hacer shuffle de datos de entrenamiento
            
        Returns:
            Tuple de (train_loader, test_loader)
        """
        if self.graph_data is None:
            raise ValueError("Debe cargar los datos de grafos primero con load_graph_genes()")
        
        if self.gene_dict is None:
            raise ValueError("Debe cargar el diccionario de genes primero con load_graph_genes()")

        if test_gene not in self.gene_dict:
            raise ValueError(f"test_gene ({test_gene}) no encontrado en los datos de grafos")

        # Obtener datos de test y train
        test_idx = self.gene_dict[test_gene]
        all_indices = set(range(len(self.graph_data)))
        test_indices = set(test_idx)
        train_indices = list(all_indices - test_indices)
        
        # Crear lista de grafos de entrenamiento
        train_graphs = [self.graph_data[i] for i in train_indices]
        test_graphs = [self.graph_data[i] for i in test_indices]

        # Crear DataLoaders
        train_loader = DataLoader(train_graphs, 
                                 batch_size=batch_size_train, 
                                 shuffle=shuffle_train)
        test_loader = DataLoader(test_graphs, 
                                batch_size=batch_size_test, 
                                shuffle=False)

        logger.info(f"Entrenamiento para gen {test_gene}: {len(train_graphs)} grafos train y {len(test_graphs)} grafos test")

        return train_loader, test_loader

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
        else:
            raise ValueError(f"Tipo de modelo no soportado: {model_type}")

    def define_checkpoint_callback(self, model_name: str, gene_name: str) -> ModelCheckpoint:
        """
        Define callback para guardar el mejor modelo por fold.
        
        Args:
            model_name: Nombre del modelo
            fold_idx: Índice del fold
            
        Returns:
            ModelCheckpoint callback
        """
        checkpoint_dir = self.results_dir / "model_checkpoints" / f"gene_{gene_name}"
        checkpoint_dir.mkdir(parents=True, exist_ok=True)

        return ModelCheckpoint(
            dirpath=str(checkpoint_dir),
            filename=f"{model_name}_gene_{gene_name}",
            save_top_k=self.config['checkpoint_config'].get('save_top_k', 1),
            monitor=self.config['checkpoint_config'].get('monitor', 'val_c_index'),
            mode=self.config['checkpoint_config'].get('mode', 'max')
        )
    
    def train_single_gene(self,
                         model: L.LightningModule,
                         train_loader: DataLoader,
                         test_loader: DataLoader,
                         gene_name: str,
                         model_name: str,
                         max_epochs: int = 30) -> Dict[str, Any]:
        """
        Entrena un modelo en un solo gen.
        
        Args:
            model: Modelo Lightning a entrenar
            train_loader: DataLoader de entrenamiento
            test_loader: DataLoader de test
            gene_name: Nombre del gen para logging
            model_name: Nombre del modelo para logging
            max_epochs: Número máximo de épocas
            
        Returns:
            Dict con métricas del gen
        """
        # Logger para este gen
        logger_name = f"{model_name}_gene_{gene_name}"
        logger_tb = TensorBoardLogger(self.lightning_logs_dir, name=logger_name)
        
        # Checkpoint callback
        checkpoint_callback = self.define_checkpoint_callback(model_name, gene_name)

        # Trainer
        trainer = L.Trainer(
            max_epochs=max_epochs,
            logger=logger_tb,
            callbacks=[checkpoint_callback],
            enable_progress_bar=True, 
            enable_model_summary=False
        )
        
        # Entrenar
        logger.info(f"Entrenando {model_name} - Gen {gene_name}")
        trainer.fit(model=model, 
                   train_dataloaders=train_loader, 
                   val_dataloaders=test_loader)
        
        # Extraer métricas finales
        # Extraer el mejor c_index del checkpoint callback (si está disponible)
        best_c_index = None
        if hasattr(checkpoint_callback, 'best_model_score') and checkpoint_callback.best_model_score is not None:
            best_c_index = checkpoint_callback.best_model_score.item()
        
        metrics = {
            'gene': gene_name,
            'final_val_loss': trainer.logged_metrics.get('val_loss', None).item(),
            'final_val_c_index': trainer.logged_metrics.get('val_c_index', None).item(),
            'best_val_c_index': best_c_index,
            'checkpoint_path': checkpoint_callback.best_model_path
        }

        logger.info(f"Gene {gene_name} completado. C-index val: {metrics['final_val_c_index']}")

        return metrics



    def run_gene_evaluation(self,
                            model_name: str,
                            min_patients: int = 10,
                            model_config: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
        """
        Entrena un modelo quitando un gen usando configuración YAML.
        
        Args:
            model_name: Nombre del modelo en el archivo de configuración
            model_config: Configuración manual (opcional, sobrescribe YAML)
            
        Returns:
            Dict con resultados de CV
        """
        
        if self.graph_data is None:
            self.load_graph_genes()
        
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
        
        gene_results = []

        sel_genes = [key for key, value in self.gene_dict.items() if len(value) > min_patients]

        logger.info(f"Número de genes seleccionados para {model_name}: {len(sel_genes)}")

        sel_genes = ['DNMT3A']
        # Iterar sobre cada fold como test
        for gene in sel_genes:
            logger.info(f"Procesando gen {gene}")

            test_size = len(self.gene_dict[gene])
            train_size = len(self.graph_data) - test_size

            while train_size % batch_size_train == 1:
                batch_size_train -= 1
            while test_size % batch_size_test == 1:
                batch_size_test -= 1
            logger.info(f"Batch sizes ajustados: {batch_size_train}/{batch_size_test} (train/test)")
            # Crear DataLoaders para este gen
            train_loader, test_loader = self.create_gene_dataloaders(
                test_gene=gene,
                batch_size_train=batch_size_train,
                batch_size_test=batch_size_test
            )
            
            # Crear instancia del modelo (nueva para cada fold)
            model = self.create_model_instance(model_name, model_config)
            
            # Entrenar en este fold
            gene_metrics = self.train_single_gene(
                model=model,
                train_loader=train_loader,
                test_loader=test_loader,
                gene_name=gene,
                model_name=model_name,
                max_epochs=max_epochs
            )

            gene_results.append(gene_metrics)
        
        gene_results = {
            'model_name': model_name,
            'model_config': config,
            'gene_results': gene_results,
            
            # Información adicional
            'training_config': training_config
        }
        
        # Guardar resultados
        self.save_gene_results(gene_results, model_name)
        
        logger.info(f"CV completado para {model_name}")

        return gene_results

    def save_gene_results(self, gene_results: Dict[str, Any], model_name: str):
        """
        Guarda los resultados de los genes en archivo.

        Args:
            gene_results: Resultados de los genes
            model_name: Nombre del modelo
        """
        # Guardar como pickle
        results_file = self.results_dir / f"{model_name}_gene_results.pkl"
        with open(results_file, 'wb') as f:
            pickle.dump(gene_results, f)

        # Guardar resumen como CSV
        fold_df = pd.DataFrame(gene_results['gene_results'])
        summary_file = self.results_dir / f"{model_name}_gene_summary.csv"
        fold_df.to_csv(summary_file, index=False)
        
        logger.info(f"Resultados guardados en: {results_file}")
        logger.info(f"Resumen guardado en: {summary_file}")

    
# Funciones de conveniencia
def run_gnn_gene_pipeline(graph_data_path, config_path, results_dir, lightning_logs_dir, models_to_run=None):
 
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
    gene_trainer = GNNGeneTrainer(
        graph_data_path=graph_data_path,
        config_path=config_path,
        results_dir=results_dir,
        lightning_logs_dir=lightning_logs_dir
    )
    
    # Determinar modelos a ejecutar
    if models_to_run is None:
        # Usar algunos modelos por defecto si no se especifica
        models_to_run = list(gene_trainer.config['models'].keys())   

    logger.info(f"Ejecutando Genes para modelos: {models_to_run}")

    # Ejecutar Genes para cada modelo
    all_results = {}

    for model_name in models_to_run:
        if model_name not in gene_trainer.config['models']:
            logger.warning(f"Modelo {model_name} no encontrado en configuración, saltando")
            continue
        try:
            results = gene_trainer.run_gene_evaluation(model_name)
            all_results[model_name] = results
            
        except Exception as e:
            logger.error(f"Error ejecutando Genes para {model_name}: {str(e)}")
            continue
    logger.info("Pipeline de GNN Genes completado")
    return all_results  