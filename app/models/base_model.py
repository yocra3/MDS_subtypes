"""
Clase base abstracta para todos los modelos de predicción de riesgo MDS.

Este módulo define la interfaz común que deben implementar todos los modelos
de predicción, así como las estructuras de datos compartidas.
"""

from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Dict, List, Optional


@dataclass
class PatientData:
    """
    Estructura de datos del paciente para entrada a los modelos.
    
    Contiene todas las variables clínicas, citogenéticas y mutacionales
    necesarias para la predicción de riesgo MDS.
    """
    # Variables clínicas continuas
    AGE: float
    BM_BLAST: float  # % blastos en médula ósea (0-100)
    WBC: float       # White blood cells (×10⁹/L)
    ANC: float       # Absolute neutrophil count (×10⁹/L)
    HB: float        # Hemoglobina (g/dL)
    PLT: float       # Plaquetas (×10⁹/L)
    MONOCYTES: float # Monocitos (×10⁹/L)

    # Variables categóricas
    SEX: str  # 'M' o 'F'
    CYTO_IPSSR: str  # 'Very Good', 'Good', 'Intermediate', 'Poor', 'Very Poor'
    complex: str    # Cariotipo complejo (≥3 alteraciones)

    # Variables binarias citogenéticas
    plus8: bool      # Trisomía 8
    del7: bool       # Deleción cromosoma 7
    del20q: bool     # Deleción 20q
    del7q: bool      # Deleción 7q
    delY: bool       # Deleción cromosoma Y
    del5q: bool      # Deleción 5q
    
    # Lista de mutaciones del paciente
    mutations: List[str]  # gen -> VAF (0-100)


@dataclass
class RiskPrediction:
    """
    Resultado de predicción de riesgo de un modelo.
    
    Contiene el score numérico, categoría de riesgo y métricas adicionales
    para interpretación clínica.
    """
    raw_score: float         # Score numérico del modelo
    risk_category: str       # 'Very Low', 'Low', 'Moderate Low', 'Moderate High', 'High', 'Very High'
    #risk_probability: float  # Probabilidad interpretable (0-1)
    #hazard_ratio: float      # Hazard ratio vs paciente promedio
    model_version: str       # Versión/nombre del modelo usado
    #confidence: float        # Confianza en la predicción (0-1)
    
    # Metadata opcional
    processing_time: Optional[float] = None  # Tiempo de procesamiento en segundos
    warnings: Optional[List[str]] = None     # Advertencias sobre la predicción


@dataclass
class ValidationResult:
    """
    Resultado de validación de datos de entrada.
    
    Indica si los datos son válidos y proporciona detalles sobre
    errores y advertencias encontradas.
    """
    is_valid: bool           # True si los datos pasaron todas las validaciones
    errors: List[str]        # Lista de errores que impiden la predicción
    warnings: List[str]      # Lista de advertencias que no impiden la predicción
    
    # Metadata de validación
    validated_fields: Optional[List[str]] = None  # Campos que fueron validados
    missing_fields: Optional[List[str]] = None    # Campos requeridos faltantes


class BaseModel(ABC):
    """
    Clase base abstracta para todos los modelos de predicción de riesgo MDS.
    
    Define la interfaz común que deben implementar todos los modelos:
    - Validación de datos de entrada
    - Predicción de riesgo
    - Categorización de resultados
    
    Todos los modelos específicos (GNN, Legacy, etc.) deben heredar de esta clase.
    """
    
    @abstractmethod
    def predict(self, patient_data: PatientData) -> RiskPrediction:
        """
        Predice el riesgo para un paciente dado.
        
        Args:
            patient_data: Datos del paciente a evaluar
            
        Returns:
            RiskPrediction: Predicción completa de riesgo
            
        Raises:
            ValueError: Si los datos de entrada no son válidos
            RuntimeError: Si hay errores durante la predicción
        """
        pass
    
    @abstractmethod
    def validate_input(self, patient_data: PatientData) -> ValidationResult:
        """
        Valida los datos de entrada del paciente.
        
        Args:
            patient_data: Datos del paciente a validar
            
        Returns:
            ValidationResult: Resultado detallado de la validación
        """
        pass
    
    def get_risk_category(self, raw_score: float) -> str:
        """
        Convierte un score numérico a categoría de riesgo IPSS-M.
        
        Implementación por defecto basada en thresholds estándar IPSS-M.
        Los modelos específicos pueden sobrescribir este método.
        
        Args:
            raw_score: Score numérico del modelo
            
        Returns:
            str: Categoría de riesgo
        """
        if raw_score <= -1.5:
            return 'Very Low'
        elif raw_score <= -0.5:
            return 'Low'
        elif raw_score <= 0:
            return 'Moderate Low'
        elif raw_score <= 0.5:
            return 'Moderate High'
        elif raw_score <= 1.5:
            return 'High'
        else:
            return 'Very High'
    
    def get_model_info(self) -> Dict[str, str]:
        """
        Retorna información básica sobre el modelo.
        
        Método opcional que los modelos pueden sobrescribir para
        proporcionar metadata adicional.
        
        Returns:
            Dict con información del modelo (nombre, versión, etc.)
        """
        return {
            'name': self.__class__.__name__,
            'type': 'base',
            'version': 'unknown'
        }


# Funciones de utilidad para trabajar con las estructuras de datos

def create_empty_patient_data() -> PatientData:
    """
    Crea una instancia de PatientData con valores por defecto.
    
    Útil para inicialización de formularios o testing.
    
    Returns:
        PatientData con valores por defecto
    """
    return PatientData(
        AGE=65.0,
        BM_BLAST=5.0,
        WBC=3.0,
        ANC=1.5,
        HB=10.0,
        PLT=150.0,
        MONOCYTES=0.3,
        SEX='M',
        CYTO_IPSSR='Intermediate',
        plus8=False,
        del7=False,
        del20q=False,
        del7q=False,
        delY=False,
        del5q=False,
        complex="complex",
        mutations={}
    )


def create_validation_result(is_valid: bool = True, 
                           errors: Optional[List[str]] = None,
                           warnings: Optional[List[str]] = None) -> ValidationResult:
    """
    Función helper para crear ValidationResult.
    
    Args:
        is_valid: Si la validación pasó
        errors: Lista de errores (opcional)
        warnings: Lista de advertencias (opcional)
        
    Returns:
        ValidationResult configurado
    """
    return ValidationResult(
        is_valid=is_valid,
        errors=errors or [],
        warnings=warnings or []
    )
