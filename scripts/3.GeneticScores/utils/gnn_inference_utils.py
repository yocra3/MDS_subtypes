import yaml
import json
import torch
import glob
import os
from pathlib import Path
from typing import Dict, Any, Tuple
import logging

logger = logging.getLogger(__name__)

def load_config(config_path: str) -> Dict[str, Any]:
    """
    Cargar configuración desde archivo YAML o JSON.
    
    Args:
        config_path (str): Ruta al archivo de configuración
        
    Returns:
        Dict[str, Any]: Configuración cargada
        
    Raises:
        ValueError: Si el formato no es soportado
        Exception: Si hay error leyendo el archivo
    """
    config_path = Path(config_path)
    
    try:
        with open(config_path, 'r', encoding='utf-8') as f:
            if config_path.suffix.lower() in ['.yaml', '.yml']:
                config = yaml.safe_load(f)
            elif config_path.suffix.lower() == '.json':
                config = json.load(f)
            else:
                raise ValueError(f"Formato de configuración no soportado: {config_path.suffix}")
        
        logger.info(f"Configuración cargada exitosamente desde {config_path}")
        return config
        
    except Exception as e:
        logger.error(f"Error cargando configuración desde {config_path}: {e}")
        raise

def get_graph_dimensions(graph_data: Dict[str, Any]) -> Tuple[int, int]:
    """
    Extrae las dimensiones de características de pacientes y genes del primer grafo.
    
    Args:
        graph_data (Dict[str, Any]): Datos de grafos cargados
        
    Returns:
        Tuple[int, int]: Dimensiones de características (paciente, gene)
    """
    first_fold_key = list(graph_data.keys())[0]
    first_graph = graph_data[first_fold_key][0]
    
    patient_feat_dim = first_graph['patient'].x.shape[1]
    gene_feat_dim = first_graph['gene'].x.shape[1]
    
    return patient_feat_dim, gene_feat_dim

def find_checkpoint(checkpoint_path: str, model_name: str, fold: int) -> str:
    """
    Busca y selecciona el checkpoint más reciente para un modelo y fold específicos.
    
    Args:
        checkpoint_path (str): Ruta base de los checkpoints
        model_name (str): Nombre del modelo
        fold (int): Número del fold
        
    Returns:
        str: Ruta al checkpoint seleccionado
        
    Raises:
        FileNotFoundError: Si no se encuentran checkpoints
    """
    base_checkpoint_path = f"{checkpoint_path}/fold_{fold}/{model_name}_fold_{fold}"
    checkpoint_pattern = f"{base_checkpoint_path}*.ckpt"
    
    checkpoint_files = glob.glob(checkpoint_pattern)
    
    if not checkpoint_files:
        raise FileNotFoundError(f"No se encontraron checkpoints para el patrón: {checkpoint_pattern}")
    
    # Seleccionar el archivo más reciente basado en la fecha de modificación
    checkpoint = max(checkpoint_files, key=os.path.getmtime)
    print(f"Using checkpoint: {checkpoint}")
    
    return checkpoint

def instantiate_gnn_model(model_name: str, config: Dict[str, Any], 
                         patient_feat_dim: int, gene_feat_dim: int, 
                         checkpoint_path: str):
    """
    Instancia un modelo GNN específico basado en la configuración.
    
    Args:
        model_name (str): Nombre del modelo
        config (Dict[str, Any]): Configuración del modelo
        patient_feat_dim (int): Dimensión de características de paciente
        gene_feat_dim (int): Dimensión de características de genes
        checkpoint_path (str): Ruta al checkpoint del modelo
        
    Returns:
        model: Instancia del modelo GNN
        
    Raises:
        ValueError: Si el tipo de modelo no es reconocido
    """
    # Importar clases de evaluación GNN
    from evalGNN_classes import EvalGNN, EvalGNN2, EvalGNN3, EvalGNN4
    
    model_config = config['models'].get(model_name, {})
    model_type = model_config['model_type']
    use_vaf = model_config['use_vaf']
    hyperparams = model_config['hyperparameters']
    
    # Parámetros base con dimensiones automáticas
    base_params = {
        'patient_feat_dim': patient_feat_dim,
        'gene_feat_dim': gene_feat_dim,
        'hidden_gene_dim': hyperparams.get('hidden_gene_dim', 32),
        'hidden_dim': hyperparams.get('hidden_dim', 16),
        'out_dim': hyperparams.get('out_dim', 1),
        'learning_rate': hyperparams.get('learning_rate', 0.001),
        'use_vaf': use_vaf,
        'hidden_vaf': hyperparams.get('hidden_vaf', 8)
    }
    
    # Agregar parámetros específicos para modelo v4
    if model_type == 'v4':
        base_params.update({
            'dropout': hyperparams.get('dropout', 0.3),
            'l2_reg': hyperparams.get('l2_reg', 1e-4)
        })
    
    # Instanciar el modelo basado en el tipo
    if model_type == "v1":
        model = EvalGNN(ck_path=checkpoint_path, **base_params)
    elif model_type == "v2":
        model = EvalGNN2(ck_path=checkpoint_path, **base_params)
    elif model_type == "v3":
        model = EvalGNN3(ck_path=checkpoint_path, **base_params)
    elif model_type == "v4":
        model = EvalGNN4(ck_path=checkpoint_path, **base_params)
    else:
        raise ValueError(f"Unknown model type: {model_type}")
    
    return model
