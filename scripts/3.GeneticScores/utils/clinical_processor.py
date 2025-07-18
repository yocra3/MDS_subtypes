"""
Módulo para procesamiento de variables clínicas.
"""
import pandas as pd
import numpy as np
import logging
import json
from typing import Dict, List, Any, Tuple
from sklearn.preprocessing import StandardScaler, OneHotEncoder
from sklearn.compose import ColumnTransformer

logger = logging.getLogger(__name__)

class ClinicalProcessor:
    """Clase para procesar variables clínicas."""
    
    def __init__(self, config: Dict[str, Any]):
        """Inicializa el procesador clínico."""
        self.config = config
        self.var_processing = config['variable_processing']
        self.fitted_transformer = None
        
        logger.info("ClinicalProcessor inicializado")


    def prepare_clinical_features(self, df_clinical: pd.DataFrame, selected_genes: List[str]) -> Tuple[pd.DataFrame, Dict[str, Any], ColumnTransformer]:
        """
        Preparar características clínicas con filtro obligatorio de casos completos.
        """
        logger.info("Iniciando procesamiento de características clínicas...")
        
        # FILTRO OBLIGATORIO: Seleccionar solo casos completos basado en configuración
        df_complete, filtering_metadata = self._filter_complete_cases_from_config(df_clinical)
        

        # Procesamiento existente de variables clínicas
       # df_processed = self._apply_existing_processing(df_complete, selected_genes)

        transformer = self._define_processing_pipeline(df_complete, selected_genes)

        # Actualizar metadata
        metadata = self._generate_processing_metadata(df_complete, df_clinical, selected_genes)
        metadata['complete_cases_filtering'] = filtering_metadata
        
        return df_complete, metadata, transformer

    def _define_processing_pipeline(self, df_clinical: pd.DataFrame, selected_genes: List[str]) -> ColumnTransformer:
        """Aplica el procesamiento existente de variables clínicas."""
        logger.info("Aplicando procesamiento de variables clínicas...")
              
        # Extraer variables por tipo
        continuous_vars = self.var_processing.get('continuous_variables', [])
        categorical_vars = self.var_processing.get('categorical_variables', [])
        binary_vars = self.var_processing.get('binary_variables', [])
        binary_vars = binary_vars + selected_genes
        transformers = []

        # Procesar variables continuas con StandardScaler
        if continuous_vars:
            scaler = StandardScaler()
            cont_scaler = scaler.fit(df_clinical[continuous_vars])
            transformers.append(('continuous', scaler, continuous_vars))
        
        # Procesar variables categóricas con OneHotEncoder
        if categorical_vars:
            encoder = OneHotEncoder(
                drop='first',  # Evitar multicolinealidad
                sparse_output=False,
                handle_unknown='ignore'
            )
            cat_encoder = encoder.fit(df_clinical[categorical_vars])
            transformers.append(('categorical', encoder, categorical_vars))
        # Procesar variables binarias
        if binary_vars:
            available_binary = [var for var in binary_vars if var in df_clinical.columns]
            if available_binary:
                transformers.append(('binary', 'passthrough', available_binary))
        
        # Crear ColumnTransformer unificado
        if transformers:
            unified_pipeline = ColumnTransformer(
                transformers=transformers,
                remainder='drop'  # Mantener otras columnas
            )
            unified_pipeline.fit(df_clinical)
            # Guardar pipeline unificada        
            logger.info("Pipeline unificada creada con transformers individuales")

        return unified_pipeline



    def _apply_existing_processing(self, df_clinical: pd.DataFrame, selected_genes: List[str]) -> pd.DataFrame:
        """Aplica el procesamiento existente de variables clínicas."""
        logger.info("Aplicando procesamiento de variables clínicas...")
        
        df_processed = df_clinical.copy()
        
        # Extraer variables por tipo
        continuous_vars = self.var_processing.get('continuous_variables', [])
        categorical_vars = self.var_processing.get('categorical_variables', [])
        binary_vars = self.var_processing.get('binary_variables', [])
        
        # Procesar variables continuas con StandardScaler
        if continuous_vars:
            df_processed = self._process_continuous_variables(df_processed, continuous_vars)
        
        # Procesar variables categóricas con OneHotEncoder
        if categorical_vars:
            df_processed = self._process_categorical_variables(df_processed, categorical_vars)
        
        self._create_unified_pipeline(df_clinical, continuous_vars, categorical_vars, binary_vars)

        return df_processed

    def _filter_complete_cases_from_config(self, df: pd.DataFrame) -> Tuple[pd.DataFrame, Dict[str, Any]]:
        """
        Filtro obligatorio de casos completos usando variables de configuración.
        """
        initial_count = len(df)
        logger.info(f"Aplicando filtro obligatorio de casos completos...")
        logger.info(f"Muestras iniciales: {initial_count}")
        
        # Extraer variables de la configuración
        model_variables = []
        
        # 1. Variables continuas desde configuración
        continuous_vars = self.var_processing.get('continuous_variables', [])
        model_variables.extend(continuous_vars)
        
        # 2. Variables categóricas desde configuración
        categorical_vars = self.var_processing.get('categorical_variables', [])
        model_variables.extend(categorical_vars)
        
        # 3. Variables binarias desde configuración
        binary_vars = self.var_processing.get('binary_variables', [])
        model_variables.extend(binary_vars)

        # 4. Genes seleccionados desde configuración
        selected_genes = self.var_processing.get('gene_selection', [])
        if isinstance(selected_genes, list):
            model_variables.extend(selected_genes)

        logger.info(f"Variables críticas identificadas desde configuración: {len(model_variables)} variables")

        # Verificar que las variables existen en el DataFrame
        available_vars = [var for var in model_variables if var in df.columns]
        missing_vars = [var for var in model_variables if var not in df.columns]
        
        if missing_vars:
            logger.warning(f"Variables esperadas no encontradas en datos: {missing_vars}")
        
        # Analizar patrones de valores faltantes antes del filtro
        missing_analysis = self._analyze_missing_patterns(df[available_vars])

        # Aplicar filtro de casos completos
        df_complete = df.dropna(subset=available_vars)
        final_count = len(df_complete)
        removed_count = initial_count - final_count
        
        # Metadata del filtrado
        filtering_metadata = {
            'initial_samples': initial_count,
            'final_samples': final_count,
            'removed_samples': removed_count,
            'removal_rate': removed_count / initial_count if initial_count > 0 else 0,
            'model_variables_used': available_vars,
            'missing_variables': missing_vars,
            'missing_analysis': missing_analysis,
            'config_source': self._get_config_source_info()
        }
               
        logger.info(f"Casos completos: {initial_count} → {final_count} muestras")
        logger.info(f"Muestras eliminadas: {removed_count} ({filtering_metadata['removal_rate']:.1%})")
        
        if removed_count > initial_count * 0.5:  # Advertencia si se elimina >50%
            logger.warning(f"⚠ Se eliminó una proporción alta de muestras: {filtering_metadata['removal_rate']:.1%}")
            logger.warning("Considerar revisar la calidad de los datos de entrada o la configuración")
        
        return df_complete, filtering_metadata

   
    def _get_config_source_info(self) -> Dict[str, Any]:
        """
        Obtener información sobre la fuente de configuración para metadata.
        """
        gene_selection = self.var_processing.get('gene_selection', [])
        genes_count = (
            len(gene_selection) if isinstance(gene_selection, list) else
            len(gene_selection.get('selected_genes', []))
        )
        
        return {
            'continuous_variables_count': len(self.var_processing.get('continuous_variables', [])),
            'categorical_variables_count': len(self.var_processing.get('categorical_variables', [])),
            'binary_variables_count': len(self.var_processing.get('binary_variables', [])),
            'selected_genes_count': genes_count,
            'config_sections_used': list(self.var_processing.keys())
        }

    def _analyze_missing_patterns(self, df: pd.DataFrame) -> Dict[str, Any]:
        """
        Analizar patrones de valores faltantes antes del filtrado.
        """
        analysis = {}
        
        # Análisis por variable
        missing_by_var = df.isnull().sum()
        missing_percentages = (missing_by_var / len(df) * 100).round(2)
        
        analysis['by_variable'] = {
            'counts': missing_by_var.to_dict(),
            'percentages': missing_percentages.to_dict()
        }
        
        # Variables con valores faltantes
        vars_with_missing = missing_by_var[missing_by_var > 0]
        analysis['variables_with_missing'] = vars_with_missing.to_dict()
        
        # Análisis por muestra
        missing_by_sample = df.isnull().sum(axis=1)
        analysis['by_sample'] = {
            'mean_missing_per_sample': missing_by_sample.mean(),
            'max_missing_per_sample': missing_by_sample.max(),
            'samples_with_no_missing': (missing_by_sample == 0).sum(),
            'samples_with_any_missing': (missing_by_sample > 0).sum()
        }
        
        # Top variables con más valores faltantes
        top_missing = missing_percentages.sort_values(ascending=False).head(10)
        analysis['top_missing_variables'] = top_missing.to_dict()
        
        # Log resumen de patrones
        if len(vars_with_missing) > 0:
            logger.info(f"Variables con valores faltantes:")
            for var, count in vars_with_missing.head(10).items():
                percentage = missing_percentages[var]
                logger.info(f"  - {var}: {count} muestras ({percentage}%)")
            if len(vars_with_missing) > 10:
                logger.info(f"  ... y {len(vars_with_missing) - 10} variables más")
        else:
            logger.info("No se encontraron valores faltantes en las variables de modelo")
        
        return analysis

    def _generate_processing_metadata(self, df_initial: pd.DataFrame, df_final: pd.DataFrame, selected_genes: List[str]) -> Dict[str, Any]:
        """Genera metadata del procesamiento completo."""
        # Extraer variables por tipo
        continuous_vars = self.var_processing.get('continuous_variables', [])
        categorical_vars = self.var_processing.get('categorical_variables', [])
        binary_vars = self.var_processing.get('binary_variables', [])

        # Obtener lista final de features
        all_binary_vars = binary_vars + selected_genes
        feature_columns = self._get_feature_columns(df_final, continuous_vars, categorical_vars, all_binary_vars)

        metadata = {
            'feature_columns': feature_columns,
            'selected_genes': selected_genes,
            'processing_summary': {
                'initial_samples': len(df_initial),
                'final_samples': len(df_final),
                'feature_count': len(feature_columns),
                'continuous_variables': continuous_vars,
                'categorical_variables': categorical_vars,
                'binary_variables': binary_vars
            }
        }
        
        logger.info(f"Features preparadas: {len(feature_columns)} variables")
        return metadata
    
    def _process_continuous_variables(self, df: pd.DataFrame, continuous_vars: List[str]) -> pd.DataFrame:
        """Procesa variables continuas con StandardScaler."""
        available_vars = [var for var in continuous_vars if var in df.columns]
        
        if not available_vars:
            return df
        
        logger.info(f"Escalando {len(available_vars)} variables continuas con StandardScaler")
        
        # Aplicar StandardScaler
        scaler = StandardScaler()
        df_scaled = df.copy()
        df_scaled[available_vars] = scaler.fit_transform(df[available_vars])
        
        # Guardar scaler
        self.fitted_transformers['continuous_scaler'] = scaler
        
        return df_scaled
    
    def _process_categorical_variables(self, df: pd.DataFrame, categorical_vars: List[str]) -> pd.DataFrame:
        """Procesa variables categóricas con OneHotEncoder de sklearn."""
        available_vars = [var for var in categorical_vars if var in df.columns]
        
        if not available_vars:
            return df
        
        logger.info(f"Codificando {len(available_vars)} variables categóricas con OneHotEncoder")
        
        df_encoded = df.copy()
        
        # Crear y entrenar OneHotEncoder
        encoder = OneHotEncoder(
            drop='first',  # Evitar multicolinealidad
            sparse_output=False,
            handle_unknown='ignore'
        )
        
        # Entrenar el encoder con todas las variables categóricas
        encoded_data = encoder.fit_transform(df[available_vars])
        
        # Obtener nombres de las características
        feature_names = encoder.get_feature_names_out(available_vars)
        
        # Crear DataFrame con los datos codificados
        encoded_df = pd.DataFrame(encoded_data, columns=feature_names, index=df.index)
        
        # Eliminar variables originales y añadir las codificadas
        df_encoded = df_encoded.drop(columns=available_vars)
        df_encoded = pd.concat([df_encoded, encoded_df], axis=1)
        
        # Guardar el encoder
        self.fitted_transformers['categorical_encoder'] = encoder
        
        return df_encoded
    
    def _get_feature_columns(self, df: pd.DataFrame, continuous_vars: List[str], 
                            categorical_vars: List[str], binary_vars: List[str]) -> List[str]:
        """Obtiene la lista final de columnas de features."""
        feature_cols = []
        
        # Variables continuas
        feature_cols.extend([var for var in continuous_vars if var in df.columns])
        # Variables categóricas (one-hot encoding genera nuevas columnas)
        for var in categorical_vars:
            if var not in df.columns:
                # Buscar columnas con prefijo del one-hot encoding
                onehot_cols = [col for col in df.columns if col.startswith(f"{var}_")]
                feature_cols.extend(onehot_cols)
        
        # Variables binarias
        feature_cols.extend([var for var in binary_vars if var in df.columns])
        
        return feature_cols
    
    def extract_survival_data(self, df: pd.DataFrame) -> Tuple[np.ndarray, np.ndarray]:
        """Extrae variables de supervivencia."""
        survival_vars = self.var_processing['survival_variables']
        time_var = survival_vars['time_variable']
        status_var = survival_vars['status_variable']
        
        if time_var not in df.columns:
            raise ValueError(f"Variable de tiempo '{time_var}' no encontrada")
        
        if status_var not in df.columns:
            raise ValueError(f"Variable de estatus '{status_var}' no encontrada")
        
        times = df[time_var].values
        events = df[status_var].values
        
        logger.info(f"Datos de supervivencia extraídos: {len(times)} observaciones")
        logger.info(f"Eventos observados: {events.sum()}/{len(events)} ({events.mean():.2%})")
        
        return times, events
    
    def extract_cv_folds(self, df: pd.DataFrame) -> np.ndarray:
        """Extrae información de folds de cross-validation."""
        cv_config = self.config.get('cross_validation', {})
        fold_column = cv_config.get('fold_column', 'fold')
        
        if fold_column not in df.columns:
            raise ValueError(f"Columna de folds '{fold_column}' no encontrada")
        
        folds = df[fold_column].values
        
        logger.info(f"Folds de CV extraídos: {len(np.unique(folds))} folds")
        return folds
    
    def get_feature_matrix(self, df: pd.DataFrame, feature_columns: List[str]) -> np.ndarray:
        """Extrae matriz de features."""
        missing_cols = [col for col in feature_columns if col not in df.columns]
        if missing_cols:
            raise ValueError(f"Columnas de features no encontradas: {missing_cols}")
        
        X = df[feature_columns].values
        
        logger.info(f"Matriz de features extraída: shape {X.shape}")
        return X

    def _create_unified_pipeline(self, df_original: pd.DataFrame, 
        continuous_vars: List[str], 
        categorical_vars: List[str], 
        binary_vars: List[str]) -> None:
        """Crear pipeline unificada usando los transformers ya entrenados."""
        
        transformers = []
        
        # Agregar transformer continuo si existe
        if continuous_vars and 'continuous_scaler' in self.fitted_transformers:
            available_continuous = [var for var in continuous_vars if var in df_original.columns]
            if available_continuous:
                transformers.append(('continuous', self.fitted_transformers['continuous_scaler'], available_continuous))
        
        # Agregar transformer categórico si existe
        if categorical_vars and 'categorical_encoder' in self.fitted_transformers:
            available_categorical = [var for var in categorical_vars if var in df_original.columns]
            if available_categorical:
                transformers.append(('categorical', self.fitted_transformers['categorical_encoder'], available_categorical))
        
        if binary_vars:
            available_binary = [var for var in binary_vars if var in df_original.columns]
            if available_binary:
                transformers.append(('binary', 'passthrough', available_binary))
    
        # Crear ColumnTransformer unificado
        if transformers:
            unified_pipeline = ColumnTransformer(
                transformers=transformers,
                remainder='drop'  # Mantener otras columnas
            )
            
            # Guardar pipeline unificada
            self.fitted_transformers['unified_pipeline'] = unified_pipeline
            
            logger.info("Pipeline unificada creada con transformers individuales")