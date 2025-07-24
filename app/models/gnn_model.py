"""
Wrapper del modelo GNN para predicción de riesgo MDS.

Este módulo implementa la clase MDSGNNModel que actúa como interfaz unificada
para el modelo PyTorch Lightning GNN, integrando todas las funcionalidades
necesarias para la predicción de riesgo en la aplicación web.
"""

import logging
import time
import pickle
from pathlib import Path
from typing import Dict, Any, Optional
import yaml

import torch
import numpy as np
from sklearn.compose import ColumnTransformer

from models.base_model import BaseModel, PatientData, RiskPrediction, ValidationResult
from models.pytorch_models.gnn_models import EvalGNN
from utils.graph_builder import GraphBuilder, load_gene_embeddings_from_file, create_dummy_gene_embeddings

logger = logging.getLogger(__name__)


class MDSGNNModel(BaseModel):
    """
    Wrapper para el modelo GNN de predicción de riesgo MDS.
    
    Implementa la interfaz BaseModel para integrar el modelo PyTorch Lightning
    en la aplicación web, manejando toda la lógica de preprocesamiento,
    construcción de grafos y predicción de riesgo.
    """
    
    def __init__(self, config_path: str):
        """
        Inicializa el modelo GNN con configuración desde archivo YAML.
        
        Args:
            config_path: Ruta al archivo de configuración YAML
            
        Raises:
            FileNotFoundError: Si no se encuentran los archivos requeridos
            RuntimeError: Si hay errores en la carga del modelo
        """
        self.config = self._load_config(config_path)
        self.model_config = self.config['models']['gnn']
        self.data_config = self.config['data']
        self.validation_config = self.config['validation']
        self.risk_config = self.config['risk_categories']
        
        # Estado de inicialización
        self.model = None
        self.transformer = None
        self.gene_embeddings = None
        self.graph_builder = None
        
        # Inicializar componentes
        self._initialize_components()
        
        logger.info("MDSGNNModel inicializado correctamente")
    
    def _load_config(self, config_path: str) -> Dict[str, Any]:
        """Carga configuración desde archivo YAML."""
        try:
            with open(config_path, 'r', encoding='utf-8') as f:
                config = yaml.safe_load(f)
            logger.info(f"Configuración cargada desde {config_path}")
            return config
        except Exception as e:
            logger.error(f"Error cargando configuración: {e}")
            raise
    
    def _initialize_components(self):
        """
        Inicializa todos los componentes del modelo:
        - Modelo PyTorch Lightning
        - Transformer de features clínicas  
        - Gene embeddings
        - Graph builder
        """
        try:
            # 1. Cargar modelo GNN
            self._load_gnn_model()
            
            # 2. Cargar transformer clínico
            self._load_clinical_transformer()
            
            # 3. Cargar gene embeddings
            self._load_gene_embeddings()
            
            # 4. Inicializar graph builder
            self._initialize_graph_builder()
            
        except Exception as e:
            logger.error(f"Error inicializando componentes: {e}")
            raise RuntimeError(f"Fallo en inicialización del modelo: {e}")
    
    def _load_gnn_model(self):
        """Carga el modelo GNN desde checkpoint."""
        try:
            checkpoint_path = self.model_config['checkpoint_path']
            arch_config = self.model_config['architecture']
            
            # Verificar que existe el checkpoint
            if not Path(checkpoint_path).exists():
                raise FileNotFoundError(f"Checkpoint no encontrado: {checkpoint_path}")
            
            # Cargar modelo real
            self.model = EvalGNN(
                ck_path=checkpoint_path,
                patient_feat_dim=arch_config['patient_feat_dim'],
                gene_feat_dim=arch_config['gene_feat_dim'],
                hidden_gene_dim=arch_config['hidden_gene_dim'],
                hidden_dim=arch_config['hidden_dim'],
                learning_rate=self.model_config['hyperparameters']['learning_rate']
            )
            
            logger.info("Modelo GNN cargado correctamente")
            
        except Exception as e:
            logger.error(f"Error cargando modelo GNN: {e}")
            raise
    
    def _load_clinical_transformer(self):
        """Carga el transformer pre-entrenado para features clínicas."""
        try:
            transformer_path = self.data_config['clinical_transformer']['file_path']
            
            if not Path(transformer_path).exists():
                raise FileNotFoundError(f"Transformer no encontrado: {transformer_path}")
            
            # Cargar transformer real
            with open(transformer_path, 'rb') as f:
                self.transformer = pickle.load(f)
            
            logger.info("Clinical transformer cargado correctamente")
            
        except Exception as e:
            logger.error(f"Error cargando clinical transformer: {e}")
            raise
    
    def _load_gene_embeddings(self):
        """Carga los embeddings de genes."""
        try:
            embeddings_path = self.data_config['gene_embeddings']['file_path']
            supported_genes = self.data_config['gene_embeddings']['supported_genes']
            embedding_dim = self.data_config['gene_embeddings']['embedding_dim']
            
            if not Path(embeddings_path).exists():
                raise FileNotFoundError(f"Gene embeddings no encontrados: {embeddings_path}")
            
            # Cargar embeddings reales
            self.gene_embeddings = load_gene_embeddings_from_file(embeddings_path)
            
            logger.info(f"Gene embeddings cargados: {len(self.gene_embeddings)} genes")
            
        except Exception as e:
            logger.error(f"Error cargando gene embeddings: {e}")
            raise
    
    def _initialize_graph_builder(self):
        """Inicializa el constructor de grafos."""
        try:
            graph_config = {
                'genes': self.data_config['gene_embeddings']['supported_genes'],
                'embedding_dim': self.data_config['gene_embeddings']['embedding_dim']
            }
            
            self.graph_builder = GraphBuilder(
                config=graph_config,
                transformer=self.transformer
            )
            
            logger.info("Graph builder inicializado")
            
        except Exception as e:
            logger.error(f"Error inicializando graph builder: {e}")
            raise
    
    def predict(self, patient_data: PatientData) -> RiskPrediction:
        """
        Predice el riesgo para un paciente usando el modelo GNN.
        
        Args:
            patient_data: Datos del paciente
            
        Returns:
            RiskPrediction: Predicción completa de riesgo
            
        Raises:
            ValueError: Si los datos no son válidos
            RuntimeError: Si hay errores en la predicción
        """
        start_time = time.time()
        
        try:
            # 1. Validar datos de entrada
            validation = self.validate_input(patient_data)
            if not validation.is_valid:
                raise ValueError(f"Datos inválidos: {validation.errors}")
            
            # 2. Construir grafo heterogéneo
            patient_graph = self.graph_builder.build_patient_graph_from_data(
                patient_data, self.gene_embeddings
            )

            reference_score = 0.8384760022163391
            # 3. Realizar predicción
            raw_score = self._predict_with_model(patient_graph)
            score = (raw_score - reference_score) / np.log(2)
            # 4. Procesar resultados
            risk_category = self.get_risk_category(score)
            #risk_probability = self._score_to_probability(raw_score)
            #hazard_ratio = self._score_to_hazard_ratio(raw_score)
            #confidence = self._calculate_confidence(patient_data, raw_score)
            
            processing_time = time.time() - start_time
            
            return RiskPrediction(
                raw_score=float(score),
                risk_category=risk_category,
                model_version=self.get_model_info()['version'],
                processing_time=processing_time,
                warnings=validation.warnings
            )
            
        except Exception as e:
            logger.error(f"Error en predicción: {e}")
            raise RuntimeError(f"Error durante predicción: {e}")
    
    def _predict_with_model(self, patient_graph) -> float:
        """Ejecuta predicción con el modelo cargado."""
        try:
            if hasattr(self.model, 'predict'):
                # Modelo real
                prediction = self.model.predict(patient_graph)
                return float(prediction.squeeze())
                
        except Exception as e:
            logger.error(f"Error en predicción del modelo: {e}")
            raise
    
    def validate_input(self, patient_data: PatientData) -> ValidationResult:
        """
        Valida los datos de entrada del paciente.
        
        Args:
            patient_data: Datos del paciente a validar
            
        Returns:
            ValidationResult: Resultado detallado de la validación
        """
        errors = []
        warnings = []
        validated_fields = []
        
        try:
            # Validar variables clínicas numéricas
            clinical_ranges = self.validation_config['clinical_ranges']
            
            numeric_fields = {
                'AGE': patient_data.AGE,
                'BM_BLAST': patient_data.BM_BLAST,
                'WBC': patient_data.WBC,
                'ANC': patient_data.ANC,
                'HB': patient_data.HB,
                'PLT': patient_data.PLT,
                'MONOCYTES': patient_data.MONOCYTES
            }
            
            for field, value in numeric_fields.items():
                validated_fields.append(field)
                if field in clinical_ranges:
                    range_config = clinical_ranges[field]
                    if value < range_config['min'] or value > range_config['max']:
                        errors.append(
                            f"{field}: {value} fuera del rango válido "
                            f"({range_config['min']}-{range_config['max']} {range_config['unit']})"
                        )
            
            # Validar variables categóricas
            categorical_values = self.validation_config['categorical_values']
            
            if patient_data.SEX not in categorical_values['SEX']:
                errors.append(f"Sexo inválido: {patient_data.SEX}")
            validated_fields.append('SEX')
            
            if patient_data.CYTO_IPSSR not in categorical_values['CYTO_IPSSR']:
                errors.append(f"Riesgo citogenético inválido: {patient_data.CYTO_IPSSR}")
            validated_fields.append('CYTO_IPSSR')
            
            # Validar mutaciones
            supported_genes = set(self.data_config['gene_embeddings']['supported_genes'])
            unsupported_genes = [gene for gene in patient_data.mutations 
                               if gene not in supported_genes]
            
            if unsupported_genes:
                warnings.append(
                    f"Genes no soportados ignorados: {', '.join(unsupported_genes)}"
                )
            
            # Advertir si no hay mutaciones
            if not patient_data.mutations:
                warnings.append(self.validation_config['warnings']['missing_mutations']['message'])
            
            is_valid = len(errors) == 0

            return ValidationResult(
                is_valid=is_valid,
                errors=errors,
                warnings=warnings,
                validated_fields=validated_fields,
                missing_fields=[]
            )
            
        except Exception as e:
            logger.error(f"Error en validación: {e}")
            return ValidationResult(
                is_valid=False,
                errors=[f"Error en validación: {e}"],
                warnings=[],
                validated_fields=[],
                missing_fields=[]
            )
    
    def get_risk_category(self, raw_score: float) -> str:
        """
        Convierte score numérico a categoría de riesgo basada en configuración.
        
        Args:
            raw_score: Score numérico del modelo
            
        Returns:
            str: Categoría de riesgo
        """
        thresholds = self.risk_config['thresholds']
        
        if raw_score <= thresholds['very_low']:
            return 'very_low'
        elif raw_score <= thresholds['low']:
            return 'low'
        elif raw_score <= thresholds['moderate_low']:
            return 'moderate_low'
        elif raw_score <= thresholds['moderate_high']:
            return 'moderate_high'
        elif raw_score <= thresholds['high']:
            return 'high'
        else:
            return 'very_high'
    
    def get_risk_category_label(self, category: str) -> str:
        """Obtiene etiqueta en español para categoría de riesgo."""
        return self.risk_config['labels'].get(category, category)
    
    def get_risk_category_color(self, category: str) -> str:
        """Obtiene color CSS para categoría de riesgo."""
        return self.risk_config['colors'].get(category, '#666666')
    
    def _score_to_probability(self, raw_score: float) -> float:
        """Convierte score a probabilidad interpretable."""
        # Transformación sigmoide simple
        # En implementación real, usar curvas de calibración
        return 1.0 / (1.0 + np.exp(-raw_score))
    
    def _score_to_hazard_ratio(self, raw_score: float) -> float:
        """Convierte score a hazard ratio."""
        # Aproximación: HR = exp(score)
        # En implementación real, usar coeficientes calibrados
        return float(np.exp(raw_score))
    
    def _calculate_confidence(self, patient_data: PatientData, raw_score: float) -> float:
        """Calcula confianza en la predicción."""
        # Implementación simplificada
        # En realidad debería considerar incertidumbre del modelo
        base_confidence = 0.8
        
        # Reducir confianza si hay mutaciones no soportadas
        supported_genes = set(self.data_config['gene_embeddings']['supported_genes'])
        unsupported_ratio = len([g for g in patient_data.mutations 
                               if g not in supported_genes]) / max(len(patient_data.mutations), 1)
        
        confidence = base_confidence * (1 - 0.3 * unsupported_ratio)
        
        return max(0.5, min(1.0, confidence))
    
    def get_model_info(self) -> Dict[str, str]:
        """Retorna información del modelo."""
        return {
            'name': 'MDS GNN Risk Model',
            'type': 'graph_neural_network',
            'version': self.config['metadata']['model_version'],
            'framework': 'pytorch_lightning',
            'description': 'Graph Neural Network for MDS risk prediction'
        }
    
    # Métodos para desarrollo y testing
    def _create_mock_model(self):
        """Crea modelo mock para desarrollo."""
        class MockGNNModel:
            def predict(self, patient_graph):
                # Predicción mock basada en número de mutaciones
                n_mutations = patient_graph['gene'].x.shape[0]
                return np.array([0.1 * n_mutations - 0.5])
            
            def __call__(self, patient_graph):
                return float(self.predict(patient_graph))
        
        return MockGNNModel()
    
    def _create_mock_transformer(self):
        """Crea transformer mock para desarrollo."""
        class MockTransformer:
            def transform(self, X):
                # Transformación mock que normaliza valores
                return np.random.randn(X.shape[0], 20)  # 20 features
        
        return MockTransformer()


# Función de conveniencia para cargar modelo
def load_mds_gnn_model(config_path: str = "app/config.yaml") -> MDSGNNModel:
    """
    Función de conveniencia para cargar el modelo GNN.
    
    Args:
        config_path: Ruta al archivo de configuración
        
    Returns:
        MDSGNNModel: Modelo inicializado y listo para usar
    """
    return MDSGNNModel(config_path)
