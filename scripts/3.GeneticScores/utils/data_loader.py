"""
Módulo para cargar datos desde diferentes fuentes.
"""
import pandas as pd
import numpy as np
import logging
from pathlib import Path
from typing import Dict, Optional, List, Any

logger = logging.getLogger(__name__)

class DataLoader:
    """Clase para cargar datos desde diferentes fuentes."""
    
    def __init__(self, config: Dict[str, Any]):
        """Inicializa el DataLoader."""
        self.config = config
        self.data_paths = config['data_paths']
        self.base_path = Path(self.data_paths['base_path'])
        
        logger.info("DataLoader inicializado")
    
    def load_clinical_data(self) -> pd.DataFrame:
        """Carga datos clínicos y mutaciones binarizadas."""
        clinical_file = self.base_path / self.data_paths['clinical_mutations_file']
        
        if not clinical_file.exists():
            raise FileNotFoundError(f"Archivo clínico no encontrado: {clinical_file}")
        
        logger.info(f"Cargando datos clínicos desde: {clinical_file}")
        df_clinical = pd.read_csv(clinical_file, sep='\t', low_memory=False)
        
        logger.info(f"Datos clínicos cargados: {df_clinical.shape[0]} pacientes, {df_clinical.shape[1]} variables")
        
        # Validar estructura del DataFrame
        self._validate_clinical_data(df_clinical)
        
        return df_clinical
    
    def load_mutation_vaf_data(self) -> Optional[pd.DataFrame]:
        """Carga datos de mutaciones con VAF."""
        mutations_file = self.data_paths.get('mutations_vaf_file')
        
        if mutations_file is None:
            logger.info("Archivo de mutaciones VAF no configurado")
            return None
        
        mutations_path = self.base_path / mutations_file
        
        if not mutations_path.exists():
            logger.warning(f"Archivo de mutaciones VAF no encontrado: {mutations_path}")
            return None
        
        logger.info(f"Cargando datos de mutaciones VAF desde: {mutations_path}")
        df_mutations = pd.read_csv(mutations_path, sep='\t')
        
        logger.info(f"Datos de mutaciones cargados: {df_mutations.shape[0]} mutaciones")
        
        df_mutations = df_mutations.dropna(subset=['ID', 'GENE', 'VAF'])
        # Validar estructura del DataFrame
        self._validate_mutation_data(df_mutations)
        
        return df_mutations
    
    def load_gene_embeddings(self) -> Optional[Dict[str, np.ndarray]]:
        """Carga embeddings de genes."""
        embeddings_file = self.data_paths.get('gene_embeddings_file')
        
        if embeddings_file is None:
            logger.info("Archivo de embeddings de genes no configurado")
            return None
        
        embeddings_path = self.base_path / embeddings_file
        
        if not embeddings_path.exists():
            logger.warning(f"Archivo de embeddings no encontrado: {embeddings_path}")
            return None
        
        logger.info(f"Cargando embeddings de genes desde: {embeddings_path}")
        
        try:
            # Cargar CSV con genes en filas y embeddings en columnas
            embeddings_df = pd.read_csv(embeddings_path, sep="\t", index_col=0)
            
            logger.info(f"Embeddings cargados: {embeddings_df.shape[0]} genes, {embeddings_df.shape[1]} dimensiones")
            
            # Convertir DataFrame a diccionario {gene_name: embedding_array}
            gene_embeddings = {}
            for gene_name in embeddings_df.index:
                gene_embeddings[gene_name] = embeddings_df.loc[gene_name].values.astype(np.float32)
            
            # Añadir embedding dummy para "Gene_0" (gene ficticio)
            embedding_dim = embeddings_df.shape[1]
            gene_embeddings["Gene_0"] = np.zeros(embedding_dim, dtype=np.float32)
            
            logger.info(f"Embeddings procesados: {len(gene_embeddings)} genes (incluyendo Gene_0 dummy)")
            
            return gene_embeddings
            
        except Exception as e:
            logger.error(f"Error cargando embeddings: {e}")
            return None
    
    def _validate_clinical_data(self, df: pd.DataFrame) -> None:
        """Valida la estructura del DataFrame clínico."""
        required_columns = []
        
        # Variables de supervivencia
        survival_vars = self.config['variable_processing']['survival_variables']
        required_columns.extend([survival_vars['time_variable'], survival_vars['status_variable']])
        
        # Variables procesadas
        var_processing = self.config['variable_processing']
        required_columns.extend(var_processing.get('continuous_variables', []))
        required_columns.extend(var_processing.get('categorical_variables', []))
        required_columns.extend(var_processing.get('binary_variables', []))
        
        # Fold de cross-validation
        cv_config = self.config.get('cross_validation', {})
        fold_column = cv_config.get('fold_column', 'fold')
        required_columns.append(fold_column)
        
        # Verificar columnas requeridas
        missing_columns = [col for col in required_columns if col not in df.columns]
        
        if missing_columns:
            raise ValueError(f"Columnas requeridas faltantes en datos clínicos: {missing_columns}")
        
        logger.info("Validación de datos clínicos: OK")
    
    def _validate_mutation_data(self, df: pd.DataFrame) -> None:
        """Valida la estructura del DataFrame de mutaciones."""
        expected_columns = ['ID', 'GENE', 'VAF']
        
        missing_columns = [col for col in expected_columns if col not in df.columns]
        
        if missing_columns:
            raise ValueError(f"Columnas requeridas faltantes en datos de mutaciones: {missing_columns}")
        
        # Validar que VAF esté en rango válido
        if not df['VAF'].between(0, 100).all():
            logger.warning("Algunos valores de VAF están fuera del rango [0, 100]")
        
        logger.info("Validación de datos de mutaciones: OK")
