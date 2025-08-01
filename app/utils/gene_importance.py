# -*- coding: utf-8 -*-
"""
Orquestador para evaluación de importancia de genes usando valores SHAP.

Este módulo coordina todo el proceso de evaluación de importancia de genes:
- Generación de grafos para todas las combinaciones
- Ejecución de inferencia en lote
- Cálculo de valores SHAP
- Formateo de resultados para visualización
"""

import logging
import time
from typing import Dict, List, Tuple, Optional, Any
from itertools import combinations

from models.shap_calculator import SHAPCalculator
from models.base_model import PatientData

logger = logging.getLogger(__name__)


class GeneImportanceEvaluator:
    """
    Evaluador de importancia de genes que orquesta todo el proceso.
    
    Coordina la generación de grafos, inferencia batch y cálculo SHAP
    para determinar la importancia de cada gen en las predicciones.
    """
    
    def __init__(self, gnn_model, max_genes: int = 8):
        """
        Inicializa el evaluador de importancia.
        
        Args:
            gnn_model: Modelo GNN para ejecutar inferencias
            max_genes: Número máximo de genes a analizar (complejidad exponencial)
        """
        self.gnn_model = gnn_model
        self.max_genes = max_genes
        self.shap_calculator = SHAPCalculator()
        
        logger.info(f"GeneImportanceEvaluator inicializado con max_genes={max_genes}")

    def evaluate_gene_importance(self, patient_data: PatientData, gene_embeddings: Dict[str, Any]) -> Dict[str, Any]:
        """
        Función principal para evaluar importancia de genes.
        
        Args:
            patient_data: Datos del paciente
            
        Returns:
            Dict[str, Any]: Resultados completos de importancia de genes
        """
        start_time = time.time()
        
        try:
            # 1. Validar que hay genes para analizar
            patient_genes = patient_data.mutations
            if not patient_genes:
                return self._empty_result("No hay mutaciones para analizar")
            
            # 2. Limitar número de genes si es necesario
            if len(patient_genes) > self.max_genes:
                logger.warning(f"Limitando análisis a {self.max_genes} genes de {len(patient_genes)} disponibles")
                patient_genes = patient_genes[:self.max_genes]
            
            
            # 3. Evaluar complejidad computacional
            complexity_info = self.shap_calculator.get_computation_complexity(len(patient_genes))
            if complexity_info['complexity_level'] == 'muy alta':
                return self._empty_result(f"Complejidad muy alta ({complexity_info['n_combinations']} combinaciones)")
            
            # 4. Generar PatientData para todas las combinaciones
            patient_dict = {}
            input_genes = [gene for gene in patient_genes if gene != 'Gene_0']
            logger.info(f"Evaluando importancia para {len(input_genes)} genes: {input_genes}")

            combinations_list = self.shap_calculator.generate_all_combinations(input_genes)
            patient_dict['base'] = self._create_patient_data_subset(patient_data, [])
            for combo_key, gene_combo in combinations_list:
                patient_dict[combo_key] = self._create_patient_data_subset(patient_data, gene_combo)

            # 5. Ejecutar inferencia 
            scores_dict = {}
            for combo_key, modified_patient_data in patient_dict.items():
                try:
                    # Ejecutar predicción individual
                    score = self.gnn_model.predict(modified_patient_data).raw_score
                    scores_dict[combo_key] = float(score)
                except Exception as e:
                    logger.warning(f"Error en inferencia para {combo_key}: {e}")
                    scores_dict[combo_key] = 0.0

            # 6. Calcular valores SHAP
            shap_values = self.shap_calculator.compute_shap_values(scores_dict, input_genes)
            
            formatted_results = {}
            formatted_results['base'] = scores_dict['base']
            formatted_results['shap_values'] = shap_values
            formatted_results['original_genes'] = input_genes
            formatted_results['complexity_info'] = complexity_info
            # 8. Agregar información de tiempo
            total_time = time.time() - start_time
            formatted_results['computation_time'] = total_time
            formatted_results['success'] = True
            
            logger.info(f"Evaluación completada en {total_time:.2f}s")
            return formatted_results
            
        except Exception as e:
            logger.error(f"Error en evaluación de importancia: {e}")
            return self._empty_result(f"Error en cálculo: {str(e)}")

    def _create_patient_data_subset(self, patient_data: PatientData, gene_subset: List[str]) -> PatientData:
        """
        Crear nueva instancia de PatientData con subset de genes.
        
        Args:
            patient_data: Datos originales del paciente
            gene_subset: Subset de genes a incluir
            
        Returns:
            PatientData: Nueva instancia con solo los genes especificados
        """
        # Crear copia de los datos del paciente
        modified_data = PatientData(
            AGE=patient_data.AGE,
            BM_BLAST=patient_data.BM_BLAST,
            WBC=patient_data.WBC,
            ANC=patient_data.ANC,
            HB=patient_data.HB,
            PLT=patient_data.PLT,
            MONOCYTES=patient_data.MONOCYTES,
            SEX=patient_data.SEX,
            CYTO_IPSSR=patient_data.CYTO_IPSSR,
            plus8=patient_data.plus8,
            del7=patient_data.del7,
            del20q=patient_data.del20q,
            del7q=patient_data.del7q,
            delY=patient_data.delY,
            del5q=patient_data.del5q,
            complex=patient_data.complex,
            mutations=gene_subset  # Solo los genes de este subset
        )
        
        return modified_data

      
    def _empty_result(self, message: str) -> Dict[str, Any]:
        """
        Generar resultado vacío con mensaje de error.
        
        Args:
            message: Mensaje explicativo
            
        Returns:
            Dict[str, Any]: Resultado vacío estructurado
        """
        return {
            'success': False,
            'message': message,
            'genes': [],
            'shap_values': [],
            'colors': [],
            'summary': message,
            'top_risk_genes': [],
            'top_protective_genes': [],
            'total_genes_analyzed': 0,
            'computation_time': 0.0,
            'complexity_info': None
        }
    
    def get_max_recommended_genes(self) -> int:
        """
        Retorna el número máximo recomendado de genes para análisis.
        
        Returns:
            int: Número máximo recomendado
        """
        return self.max_genes
    
    def estimate_analysis_time(self, n_genes: int) -> Dict[str, Any]:
        """
        Estima el tiempo de análisis para n genes.
        
        Args:
            n_genes: Número de genes
            
        Returns:
            Dict[str, Any]: Información de tiempo estimado
        """
        complexity_info = self.shap_calculator.get_computation_complexity(n_genes)
        
        # Agregar estimaciones más detalladas
        graph_generation_time = n_genes * 0.1  # ~0.1s per gene
        inference_time = (2 ** n_genes) * 0.05  # ~0.05s per combination
        shap_calculation_time = complexity_info['estimated_time_seconds']
        
        total_estimated_time = graph_generation_time + inference_time + shap_calculation_time
        
        return {
            'n_genes': n_genes,
            'graph_generation_time': graph_generation_time,
            'inference_time': inference_time,
            'shap_calculation_time': shap_calculation_time,
            'total_estimated_time': total_estimated_time,
            'complexity_level': complexity_info['complexity_level'],
            'recommendation': complexity_info['recommendation']
        }


def create_gene_importance_evaluator(gnn_model, max_genes: int = 8) -> GeneImportanceEvaluator:
    """
    Factory function para crear evaluador de importancia de genes.
    
    Args:
        gnn_model: Modelo GNN
        max_genes: Número máximo de genes a analizar
        
    Returns:
        GeneImportanceEvaluator: Evaluador inicializado
    """
    return GeneImportanceEvaluator(gnn_model, max_genes)


# Función de demostración
def demo_gene_importance():
    """Función de demostración del evaluador de importancia."""
    from models.base_model import PatientData
    
    # Crear datos de paciente de ejemplo
    patient_data = PatientData(
        AGE=65.0,
        BM_BLAST=8.5,
        WBC=3.2,
        ANC=1.8,
        HB=9.1,
        PLT=120.0,
        MONOCYTES=0.4,
        SEX='M',
        CYTO_IPSSR='Int',
        plus8=False,
        del7=False,
        del20q=False,
        del7q=True,
        delY=False,
        del5q=False,
        complex='non-complex',
        mutations=['TP53', 'ASXL1', 'SF3B1']
    )
    
    # Crear evaluador (con modelo mock)
    class MockGNNModel:
        def predict(self, patient_data):
            # Simulación simple
            from models.base_model import RiskPrediction
            score = 0.5 + len(patient_data.mutations) * 3
            return RiskPrediction(
                raw_score=score,
                risk_category='moderate',
                hazard_ratio=1.5,
                model_version='demo',
                warnings=[]
            )
    
    mock_model = MockGNNModel()
    evaluator = GeneImportanceEvaluator(mock_model, max_genes=3)
    
    # Ejecutar análisis (esto fallará con el mock, pero demuestra la interfaz)
    print("=== Demo Gene Importance Evaluator ===")
    print(f"Paciente con mutaciones: {patient_data.mutations}")
    
    time_estimate = evaluator.estimate_analysis_time(len(patient_data.mutations))
    print(f"Tiempo estimado: {time_estimate['total_estimated_time']:.2f}s")
    print(f"Complejidad: {time_estimate['complexity_level']}")
    evaluation_result = evaluator.evaluate_gene_importance(patient_data)
    print(f"Resultado de evaluación: {evaluation_result}")
    return evaluator


if __name__ == "__main__":
    # Ejecutar demo si se ejecuta directamente
    demo_gene_importance()
