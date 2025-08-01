# -*- coding: utf-8 -*-
"""
Calculadora de valores SHAP para importancia de genes en MDS.

Este módulo implementa el cálculo de valores SHAP (SHapley Additive exPlanations)
para determinar la importancia individual de cada gen en las predicciones del modelo GNN.
"""

import logging
import numpy as np
import pandas as pd
from itertools import combinations
from typing import Dict, List, Tuple, Optional, Any
from math import factorial

logger = logging.getLogger(__name__)


class SHAPCalculator:
    """
    Calculadora de valores SHAP para análisis de importancia de genes.
    
    Implementa el algoritmo SHAP para determinar la contribución marginal
    de cada gen en las predicciones de riesgo MDS.
    """
    
    def __init__(self):
        """Inicializa el calculador SHAP."""
        self.scores = {}
        self.genes = []
        self.base_score = 0.0
        
    def compute_shap_values(self, scores_dict: Dict[str, float], genes: List[str]) -> Dict[str, float]:
        """
        Calcula valores SHAP para un conjunto de genes y sus combinaciones.
        
        Args:
            scores_dict: Diccionario con scores para cada combinación de genes
                        Formato: {'base': score, 'GENE1': score, 'GENE1_GENE2': score, ...}
            genes: Lista de genes individuales presentes
            
        Returns:
            Dict[str, float]: Valores SHAP para cada gen
        """
        self.genes = genes
        self.base_score = scores_dict.get('base', 0.0)
        
        if not genes:
            logger.warning("No hay genes para calcular valores SHAP")
            return {}
        
        logger.info(f"Calculando valores SHAP para {len(genes)} genes")
        
        if not self.validate_scores_dict(scores_dict, genes):
            logger.error("Diccionario de scores no válido, abortando cálculo SHAP")
            return {}

        self.scores = self._correct_scores_dict(scores_dict)
        shap_values = {}
        n = len(genes)
        
        for gene in genes:
            shap_value = self._calculate_shap_for_gene(gene, n)
            shap_values[gene] = shap_value
            
        logger.info(f"Valores SHAP calculados: {shap_values}")
        return shap_values
    
    def _calculate_shap_for_gene(self, target_gene: str, n: int) -> float:
        """
        Calcula el valor SHAP para un gen específico.
        
        Args:
            target_gene: Gen para el cual calcular SHAP
            n: Número total de genes
            
        Returns:
            float: Valor SHAP para el gen
        """
        other_genes = [g for g in self.genes if g != target_gene]
        shap_value = 0.0
        
        # Caso 1: Conjunto vacío (S = ∅)
        weight_empty = factorial(0) * factorial(n - 1) / factorial(n)  # = 1/n
        f_empty = self.base_score
        f_with_gene = self.scores.get(target_gene, f_empty)
        shap_value += weight_empty * (f_with_gene - f_empty)
        
        # Caso 2: Todos los subconjuntos posibles de otros genes
        for s in range(1, len(other_genes) + 1):
            # Generar todas las combinaciones de tamaño s
            for subset in combinations(other_genes, s):
                subset_list = list(subset)
                
                # Calcular peso de Shapley
                weight = factorial(s) * factorial(n - s - 1) / factorial(n)
                
                # Score sin el gen objetivo
                subset_key = self._get_combination_key(subset_list)
                f_without = self.scores.get(subset_key)
                
                # Score con el gen objetivo
                subset_with_gene = subset_list + [target_gene]
                subset_with_key = self._get_combination_key(subset_with_gene)
                f_with = self.scores.get(subset_with_key)

                # Contribución marginal ponderada
                marginal_contribution = weight * (f_with - f_without)
                shap_value += marginal_contribution
        
        return shap_value
    
    def _get_combination_key(self, gene_list: List[str]) -> str:
        """
        Genera la clave para una combinación de genes.
        
        Args:
            gene_list: Lista de genes
            
        Returns:
            str: Clave ordenada para la combinación
        """
        if not gene_list:
            return 'base'
        return '_'.join(sorted(gene_list))

    def _correct_scores_dict(self, scores_dict: Dict[str, float]) -> Dict[str, float]:
        """
        Corrige el diccionario de scores para asegurar que todas las combinaciones
        necesarias están presentes.
        
        Args:
            scores_dict: Diccionario de scores original
            
        Returns:
            Dict[str, float]: Diccionario corregido con todas las combinaciones
        """
        corrected_scores = {}

        for key, item in scores_dict.items():

            if key == 'base':
                corrected_scores['base'] = scores_dict['base']
            
            else:
                ## Divide key en genes individuales
                genes = key.split('_')
                if len(genes) == 1:
                    corrected_scores[genes[0]] = item
                else:
                    # Genera clave combinada
                    combo_key = '_'.join(sorted(genes))
                    corrected_scores[combo_key] = item
        return corrected_scores
    
    def generate_all_combinations(self, genes: List[str]) -> List[Tuple[str, List[str]]]:
        """
        Genera todas las combinaciones posibles de genes.
        
        Args:
            genes: Lista de genes disponibles
            
        Returns:
            List[Tuple[str, List[str]]]: Lista de (clave, combinación) para todas las combinaciones
        """
        combinations_list = []
        
        # Agregar caso base (sin genes)
        combinations_list.append(('base', []))
        
        # Generar todas las combinaciones de 1 a n genes
        for r in range(1, len(genes) + 1):
            for combo in combinations(genes, r):
                combo_list = list(combo)
                key = self._get_combination_key(combo_list)
                combinations_list.append((key, combo_list))
        
        logger.info(f"Generadas {len(combinations_list)} combinaciones para {len(genes)} genes")
        return combinations_list
    
  
    
    def validate_scores_dict(self, scores_dict: Dict[str, float], genes: List[str]) -> bool:
        """
        Valida que el diccionario de scores tenga todas las combinaciones necesarias.
        
        Args:
            scores_dict: Diccionario de scores
            genes: Lista de genes
            
        Returns:
            bool: True si el diccionario es válido
        """
        # Verificar que existe el score base
        if 'base' not in scores_dict:
            logger.error("Falta score base en el diccionario")
            return False
        
        # Verificar que existen scores para genes individuales
        for gene in genes:
            if gene not in scores_dict:
                logger.warning(f"Falta score individual para gen: {gene}")
        
        # Contar combinaciones esperadas vs disponibles
        expected_combinations = 2 ** len(genes)  # Incluye caso base
        available_combinations = len(scores_dict)
        
        if available_combinations < expected_combinations:
            logger.warning(f"Combinaciones disponibles ({available_combinations}) < esperadas ({expected_combinations})")
            logger.info("Se procederá con las combinaciones disponibles")
        
        return True
    
    def estimate_computation_time(self, n_genes: int) -> float:
        """
        Estima el tiempo de computación para n genes.
        
        Args:
            n_genes: Número de genes
            
        Returns:
            float: Tiempo estimado en segundos
        """
        # Estimación basada en complejidad exponencial
        # Tiempo base por combinación: ~0.01 segundos
        n_combinations = 2 ** n_genes
        estimated_time = n_combinations * 0.01
        
        return estimated_time
    
    def get_computation_complexity(self, n_genes: int) -> Dict[str, Any]:
        """
        Retorna información sobre la complejidad computacional.
        
        Args:
            n_genes: Número de genes
            
        Returns:
            Dict[str, Any]: Información de complejidad
        """
        n_combinations = 2 ** n_genes
        estimated_time = self.estimate_computation_time(n_genes)
        
        complexity_level = "baja"
        if n_genes > 10:
            complexity_level = "muy alta"
        elif n_genes > 7:
            complexity_level = "alta"
        elif n_genes > 4:
            complexity_level = "moderada"
        
        return {
            'n_genes': n_genes,
            'n_combinations': n_combinations,
            'estimated_time_seconds': estimated_time,
            'complexity_level': complexity_level,
            'recommendation': "Proceder" if n_genes <= 8 else "Considerar limitar genes"
        }


def create_shap_calculator() -> SHAPCalculator:
    """
    Factory function para crear una instancia del calculador SHAP.
    
    Returns:
        SHAPCalculator: Instancia inicializada
    """
    return SHAPCalculator()


# Funciones de utilidad
def demo_shap_calculation():
    """Función de demostración del cálculo SHAP."""
    # Datos de ejemplo
    scores = {
        'base': 0.0,
        'TP53': 0.5,
        'ASXL1': 0.3,
        'SF3B1': -0.2,
        'TP53_ASXL1': 0.7,
        'TP53_SF3B1': 0.4,
        'ASXL1_SF3B1': 0.2,
        'TP53_ASXL1_SF3B1': 0.8
    }
    
    genes = ['TP53', 'ASXL1', 'SF3B1']
    
    calculator = SHAPCalculator()
    shap_values = calculator.compute_shap_values(scores, genes)
    
    print("=== Demo SHAP Calculator ===")
    print(f"Genes analizados: {genes}")
    print(f"Valores SHAP: {shap_values}")

    return shap_values


if __name__ == "__main__":
    # Ejecutar demo si se ejecuta directamente
    demo_shap_calculation()
