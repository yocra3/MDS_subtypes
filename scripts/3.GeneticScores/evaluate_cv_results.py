#!/usr/bin/env python3
"""
Script para evaluar y comparar resultados de cross-validation de modelos GNN
===========================================================================

Este script carga, analiza y compara los resultados de cross-validation de múltiples
modelos GNN, generando métricas agregadas, visualizaciones comparativas y reportes
estadísticos.

Uso:
    python evaluate_cv_results.py --results_dir results/gnn/basic_model [opciones]

Autor: Sistema de Evaluación de Modelos GNN
Fecha: 2025
"""

import argparse
import sys
import os
import pickle
import json
import warnings
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any, Union
from dataclasses import dataclass, asdict
import logging

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import yaml

# Configurar matplotlib para evitar problemas con displays
plt.switch_backend('Agg')
warnings.filterwarnings('ignore', category=UserWarning)

# Configurar logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout),
    ]
)
logger = logging.getLogger(__name__)

@dataclass
class ModelMetrics:
    """Clase para almacenar métricas de un modelo."""
    model_name: str
    mean_c_index: float
    std_c_index: float
    median_c_index: float
    min_c_index: float
    max_c_index: float
    q25_c_index: float
    q75_c_index: float
    mean_loss: float
    std_loss: float
    n_folds: int
    fold_c_indices: List[float]
    fold_losses: List[float]
    model_config: Dict[str, Any]
    
    def to_dict(self) -> Dict[str, Any]:
        """Convierte a diccionario para serialización."""
        return asdict(self)

class CVResultsEvaluator:
    """
    Evaluador de resultados de cross-validation para modelos GNN.
    """
    
    def __init__(self, results_dir: str, output_dir: str = None):
        """
        Inicializa el evaluador.
        
        Args:
            results_dir: Directorio con los resultados de CV
            output_dir: Directorio de salida (por defecto: {results_dir}/evaluation)
        """
        self.results_dir = Path(results_dir)
        self.output_dir = Path(output_dir) if output_dir else self.results_dir / "evaluation"
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Crear subdirectorios
        (self.output_dir / "plots").mkdir(exist_ok=True)
        (self.output_dir / "reports").mkdir(exist_ok=True)
        (self.output_dir / "data").mkdir(exist_ok=True)
        
        # Almacenar resultados cargados
        self.models_data: Dict[str, Dict[str, Any]] = {}
        self.models_metrics: Dict[str, ModelMetrics] = {}
        self.external_benchmarks: Dict[str, Dict[str, Any]] = {}
        
        logger.info(f"Evaluador inicializado")
        logger.info(f"Directorio de resultados: {self.results_dir}")
        logger.info(f"Directorio de salida: {self.output_dir}")
    
    def load_cv_results(self, patterns: List[str] = None) -> Dict[str, Dict[str, Any]]:
        """
        Carga los resultados de CV desde archivos pickle y CSV.
        
        Args:
            patterns: Lista de patrones para filtrar modelos (ej: ['v1_', 'v2_'])
        
        Returns:
            Diccionario con los resultados cargados
        """
        logger.info("Cargando resultados de CV...")
        
        # Buscar archivos de resultados
        pkl_files = list(self.results_dir.glob("*_cv_results.pkl"))
        
        logger.info(f"Encontrados {len(pkl_files)} archivos .pkl")
        
        if not pkl_files:
            raise FileNotFoundError(f"No se encontraron archivos de resultados en {self.results_dir}")
        
        # Cargar desde archivos pickle (preferido)
        for pkl_file in pkl_files:
            model_name = pkl_file.stem.replace("_cv_results", "")
            
            # Filtrar por patrones si se especifican
            if patterns and not any(pattern in model_name for pattern in patterns):
                continue
                
            try:
                with open(pkl_file, 'rb') as f:
                    data = pickle.load(f)
                self.models_data[model_name] = data
                logger.info(f"Cargado {model_name} desde {pkl_file.name}")
            except Exception as e:
                logger.warning(f"Error cargando {pkl_file}: {e}")
        
       
        logger.info(f"Cargados {len(self.models_data)} modelos en total")
        return self.models_data
    
    def load_external_benchmarks(self, benchmark_paths: Dict[str, str]) -> Dict[str, Dict[str, Any]]:
        """
        Carga benchmarks externos para comparación.
        
        Args:
            benchmark_paths: Diccionario con nombres y rutas de benchmarks
        
        Returns:
            Diccionario con benchmarks cargados
        """
        logger.info("Cargando benchmarks externos...")
        
        for name, path in benchmark_paths.items():
            try:
                if not Path(path).exists():
                    logger.warning(f"Benchmark {name} no encontrado en {path}")
                    continue
                
                if path.endswith('.csv') or path.endswith('.tsv'):
                    sep = '\t' if path.endswith('.tsv') else ','
                    df = pd.read_csv(path, sep=sep)
                    
                    # Buscar columnas relevantes
                    c_index_col = None
                    fold_col = None
                    
                    for col in df.columns:
                        if 'c_index' in col.lower() or 'cindex' in col.lower():
                            c_index_col = col
                        if 'fold' in col.lower():
                            fold_col = col
                    
                    if c_index_col is None:
                        logger.warning(f"No se encontró columna de C-index en {path}")
                        continue
                    
                    # Procesar datos
                    if fold_col and 'GLOBAL' in df[fold_col].values:
                        # Formato con folds + global
                        fold_data = df[df[fold_col] != 'GLOBAL']
                        global_data = df[df[fold_col] == 'GLOBAL']
                        
                        c_indices = fold_data[c_index_col].dropna().tolist()
                        global_c_index = global_data[c_index_col].iloc[0] if len(global_data) > 0 else np.mean(c_indices)
                    else:
                        # Formato simple
                        c_indices = df[c_index_col].dropna().tolist()
                        global_c_index = np.mean(c_indices)
                    
                    benchmark_data = {
                        'name': name,
                        'c_indices': c_indices,
                        'mean_c_index': global_c_index,
                        'std_c_index': np.std(c_indices) if len(c_indices) > 1 else 0,
                        'median_c_index': np.median(c_indices) if c_indices else np.nan,
                        'min_c_index': np.min(c_indices) if c_indices else np.nan,
                        'max_c_index': np.max(c_indices) if c_indices else np.nan,
                        'q25_c_index': np.percentile(c_indices, 25) if c_indices else np.nan,
                        'q75_c_index': np.percentile(c_indices, 75) if c_indices else np.nan,
                        'iqr_c_index': np.percentile(c_indices, 75) - np.percentile(c_indices, 25) if c_indices else np.nan,
                        'range_c_index': np.max(c_indices) - np.min(c_indices) if c_indices else np.nan,
                        'n_folds': len(c_indices),
                        'source_file': path
                    }
                    
                    self.external_benchmarks[name] = benchmark_data
                    logger.info(f"Cargado benchmark {name}: C-index = {global_c_index:.4f}")
                
            except Exception as e:
                logger.warning(f"Error cargando benchmark {name}: {e}")
        
        return self.external_benchmarks
    
    def calculate_metrics(self) -> Dict[str, ModelMetrics]:
        """
        Calcula métricas agregadas para todos los modelos.
        
        Returns:
            Diccionario con métricas por modelo
        """
        logger.info("Calculando métricas agregadas...")
        
        for model_name, data in self.models_data.items():
            try:
                # Extraer C-indices y losses de los folds
                fold_results = data.get('fold_results', [])
                
                c_indices = []
                losses = []
                
                for fold_result in fold_results:
                    if isinstance(fold_result, dict):
                        c_idx = fold_result.get('c_index', fold_result.get('best_val_c_index', np.nan))
                        loss = fold_result.get('loss', fold_result.get('final_val_loss', np.nan))
                    else:
                        # Si fold_result es una lista o tupla
                        c_idx = fold_result[1] if len(fold_result) > 1 else np.nan
                        loss = fold_result[0] if len(fold_result) > 0 else np.nan
                    
                    if not pd.isna(c_idx):
                        c_indices.append(float(c_idx))
                    if not pd.isna(loss):
                        losses.append(float(loss))
                
                # Calcular estadísticas
                if c_indices:
                    mean_c_index = np.mean(c_indices)
                    std_c_index = np.std(c_indices)
                    median_c_index = np.median(c_indices)
                    min_c_index = np.min(c_indices)
                    max_c_index = np.max(c_indices)
                    q25_c_index = np.percentile(c_indices, 25)
                    q75_c_index = np.percentile(c_indices, 75)
                else:
                    mean_c_index = data.get('mean_val_c_index', np.nan)
                    std_c_index = data.get('std_val_c_index', np.nan)
                    median_c_index = np.nan
                    min_c_index = np.nan
                    max_c_index = np.nan
                    q25_c_index = np.nan
                    q75_c_index = np.nan
                
                if losses:
                    mean_loss = np.mean(losses)
                    std_loss = np.std(losses)
                else:
                    mean_loss = data.get('mean_val_loss', np.nan)
                    std_loss = data.get('std_val_loss', np.nan)
                
                # Crear objeto ModelMetrics
                metrics = ModelMetrics(
                    model_name=model_name,
                    mean_c_index=mean_c_index,
                    std_c_index=std_c_index,
                    median_c_index=median_c_index,
                    min_c_index=min_c_index,
                    max_c_index=max_c_index,
                    q25_c_index=q25_c_index,
                    q75_c_index=q75_c_index,
                    mean_loss=mean_loss,
                    std_loss=std_loss,
                    n_folds=len(c_indices) if c_indices else data.get('n_folds', 0),
                    fold_c_indices=c_indices,
                    fold_losses=losses,
                    model_config=data.get('model_config', {})
                )
                
                self.models_metrics[model_name] = metrics
                logger.info(f"Métricas calculadas para {model_name}: C-index = {mean_c_index:.4f} ± {std_c_index:.4f}")
                
            except Exception as e:
                logger.warning(f"Error calculando métricas para {model_name}: {e}")
        
        return self.models_metrics
    
    def generate_summary_table(self) -> pd.DataFrame:
        """
        Genera tabla resumen con métricas de todos los modelos.
        
        Returns:
            DataFrame con métricas agregadas
        """
        logger.info("Generando tabla resumen...")
        
        summary_data = []
        
        # Agregar modelos GNN
        for model_name, metrics in self.models_metrics.items():
            summary_data.append({
                'Model': model_name,
                'Type': 'GNN',
                'Mean_C_Index': metrics.mean_c_index,
                'Std_C_Index': metrics.std_c_index,
                'Median_C_Index': metrics.median_c_index,
                'Min_C_Index': metrics.min_c_index,
                'Max_C_Index': metrics.max_c_index,
                'Q25_C_Index': metrics.q25_c_index,
                'Q75_C_Index': metrics.q75_c_index,
                'IQR_C_Index': metrics.q75_c_index - metrics.q25_c_index if not pd.isna(metrics.q75_c_index) else np.nan,
                'Range_C_Index': metrics.max_c_index - metrics.min_c_index if not pd.isna(metrics.max_c_index) else np.nan,
                'Mean_Loss': metrics.mean_loss,
                'Std_Loss': metrics.std_loss,
                'N_Folds': metrics.n_folds,
                'Use_VAF': metrics.model_config.get('use_vaf', False),
                'Model_Type': metrics.model_config.get('model_type', 'Unknown')
            })
        
        # Agregar benchmarks externos
        for name, benchmark in self.external_benchmarks.items():
            c_indices = benchmark['c_indices']
            summary_data.append({
                'Model': name,
                'Type': 'Benchmark',
                'Mean_C_Index': benchmark['mean_c_index'],
                'Std_C_Index': benchmark['std_c_index'],
                'Median_C_Index':  benchmark['median_c_index'],
                'Min_C_Index': benchmark['min_c_index'],
                'Max_C_Index': benchmark['max_c_index'],
                'Q25_C_Index': benchmark['q25_c_index'],
                'Q75_C_Index': benchmark['q75_c_index'],
                'IQR_C_Index': benchmark['iqr_c_index'],
                'Range_C_Index': np.max(c_indices) - np.min(c_indices) if c_indices else np.nan,
                'Mean_Loss': np.nan,
                'Std_Loss': np.nan,
                'N_Folds': benchmark['n_folds'],
                'Use_VAF': np.nan,
                'Model_Type': 'External'
            })
        
        df = pd.DataFrame(summary_data)
        
        # Ordenar por C-index medio
        df = df.sort_values('Mean_C_Index', ascending=False)
        
        # Guardar tabla
        output_path = self.output_dir / "data" / "summary_table.csv"
        df.to_csv(output_path, index=False)
        logger.info(f"Tabla resumen guardada en {output_path}")
        
        return df
    
    def create_comparison_plots(self) -> Dict[str, str]:
        """
        Crea visualizaciones comparativas de los resultados.
        
        Returns:
            Diccionario con rutas de los plots generados
        """
        logger.info("Generando visualizaciones comparativas...")
        
        plot_paths = {}
        
        # Preparar datos para plotting
        all_models = list(self.models_metrics.keys()) + list(self.external_benchmarks.keys())
        
        if not all_models:
            logger.warning("No hay modelos para visualizar")
            return plot_paths
        
        # Configurar estilo
        plt.style.use('seaborn-v0_8')
        sns.set_palette("husl")
        
        # 1. Boxplot de C-index por modelo
        self._create_boxplot_comparison(plot_paths)
        
        # 2. Barplot con barras de error
        self._create_barplot_comparison(plot_paths)
        
        # 3. Violin plot
        self._create_violin_plot(plot_paths)
        
        # 4. Scatter plot de C-index vs Loss
        self._create_scatter_plot(plot_paths)
              
        logger.info(f"Generadas {len(plot_paths)} visualizaciones")
        return plot_paths
    
    def _create_boxplot_comparison(self, plot_paths: Dict[str, str]):
        """Crea boxplot comparativo de C-index."""
        try:
            plt.figure(figsize=(12, 8))
            
            # Preparar datos
            plot_data = []
            for model_name, metrics in self.models_metrics.items():
                for c_idx in metrics.fold_c_indices:
                    plot_data.append({
                        'Model': model_name,
                        'C_Index': c_idx,
                        'Type': 'GNN'
                    })
            
            # Agregar benchmarks
            for name, benchmark in self.external_benchmarks.items():
                for c_idx in benchmark['c_indices']:
                    plot_data.append({
                        'Model': name,
                        'C_Index': c_idx,
                        'Type': 'Benchmark'
                    })
            
            if not plot_data:
                logger.warning("No hay datos para boxplot")
                return
            
            df_plot = pd.DataFrame(plot_data)
            
            # Crear boxplot
            ax = sns.boxplot(data=df_plot, x='Model', y='C_Index', hue='Type')
            plt.title('Comparación de C-index por Modelo (Cross-Validation)', fontsize=14, fontweight='bold')
            plt.xlabel('Modelo', fontsize=12)
            plt.ylabel('C-index', fontsize=12)
            plt.xticks(rotation=45, ha='right')
            plt.legend(title='Tipo', bbox_to_anchor=(1.05, 1), loc='upper left')
            plt.tight_layout()
            
            # Guardar
            output_path = self.output_dir / "plots" / "boxplot_comparison.png"
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close()
            
            plot_paths['boxplot'] = str(output_path)
            logger.info(f"Boxplot guardado en {output_path}")
            
        except Exception as e:
            logger.error(f"Error creando boxplot: {e}")
    
    def _create_barplot_comparison(self, plot_paths: Dict[str, str]):
        """Crea barplot con barras de error."""
        try:
            plt.figure(figsize=(12, 8))
            
            # Preparar datos
            models = []
            means = []
            stds = []
            types = []
            
            for model_name, metrics in self.models_metrics.items():
                models.append(model_name)
                means.append(metrics.mean_c_index)
                stds.append(metrics.std_c_index)
                types.append('GNN')
            
            for name, benchmark in self.external_benchmarks.items():
                models.append(name)
                means.append(benchmark['mean_c_index'])
                stds.append(benchmark['std_c_index'])
                types.append('Benchmark')
            
            if not models:
                logger.warning("No hay datos para barplot")
                return
            
            # Crear DataFrame
            df_plot = pd.DataFrame({
                'Model': models,
                'Mean_C_Index': means,
                'Std_C_Index': stds,
                'Type': types
            })
            
            # Ordenar por C-index
            df_plot = df_plot.sort_values('Mean_C_Index', ascending=False)
            
            # Crear barplot
            ax = sns.barplot(data=df_plot, x='Model', y='Mean_C_Index', hue='Type')
            
            # Agregar barras de error
            for i, (_, row) in enumerate(df_plot.iterrows()):
                ax.errorbar(i, row['Mean_C_Index'], yerr=row['Std_C_Index'], 
                           fmt='none', color='black', capsize=3)
            
            plt.title('Comparación de C-index Promedio por Modelo', fontsize=14, fontweight='bold')
            plt.xlabel('Modelo', fontsize=12)
            plt.ylabel('C-index Promedio', fontsize=12)
            plt.xticks(rotation=45, ha='right')
            plt.legend(title='Tipo', bbox_to_anchor=(1.05, 1), loc='upper left')
            plt.tight_layout()
            
            # Guardar
            output_path = self.output_dir / "plots" / "barplot_comparison.png"
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close()
            
            plot_paths['barplot'] = str(output_path)
            logger.info(f"Barplot guardado en {output_path}")
            
        except Exception as e:
            logger.error(f"Error creando barplot: {e}")
    
    def _create_violin_plot(self, plot_paths: Dict[str, str]):
        """Crea violin plot de distribuciones."""
        try:
            plt.figure(figsize=(14, 8))
            
            # Preparar datos
            plot_data = []
            for model_name, metrics in self.models_metrics.items():
                for c_idx in metrics.fold_c_indices:
                    plot_data.append({
                        'Model': model_name,
                        'C_Index': c_idx,
                        'Type': 'GNN'
                    })
            
            for name, benchmark in self.external_benchmarks.items():
                for c_idx in benchmark['c_indices']:
                    plot_data.append({
                        'Model': name,
                        'C_Index': c_idx,
                        'Type': 'Benchmark'
                    })
            
            if not plot_data:
                logger.warning("No hay datos para violin plot")
                return
            
            df_plot = pd.DataFrame(plot_data)
            
            # Crear violin plot
            ax = sns.violinplot(data=df_plot, x='Model', y='C_Index', hue='Type')
            plt.title('Distribución de C-index por Modelo', fontsize=14, fontweight='bold')
            plt.xlabel('Modelo', fontsize=12)
            plt.ylabel('C-index', fontsize=12)
            plt.xticks(rotation=45, ha='right')
            plt.legend(title='Tipo', bbox_to_anchor=(1.05, 1), loc='upper left')
            plt.tight_layout()
            
            # Guardar
            output_path = self.output_dir / "plots" / "violin_plot.png"
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close()
            
            plot_paths['violin'] = str(output_path)
            logger.info(f"Violin plot guardado en {output_path}")
            
        except Exception as e:
            logger.error(f"Error creando violin plot: {e}")
    
    def _create_scatter_plot(self, plot_paths: Dict[str, str]):
        """Crea scatter plot de C-index vs Loss."""
        try:
            plt.figure(figsize=(10, 8))
            
            # Preparar datos (solo GNN models que tienen loss)
            models = []
            c_indices = []
            losses = []
            
            for model_name, metrics in self.models_metrics.items():
                if not pd.isna(metrics.mean_loss):
                    models.append(model_name)
                    c_indices.append(metrics.mean_c_index)
                    losses.append(metrics.mean_loss)
            
            if not models:
                logger.warning("No hay datos para scatter plot")
                return
            
            # Crear scatter plot
            plt.scatter(losses, c_indices, alpha=0.7, s=100)
            
            # Agregar etiquetas
            for i, model in enumerate(models):
                plt.annotate(model, (losses[i], c_indices[i]), 
                           xytext=(5, 5), textcoords='offset points', 
                           fontsize=10, alpha=0.8)
            
            plt.title('C-index vs Loss Promedio', fontsize=14, fontweight='bold')
            plt.xlabel('Loss Promedio', fontsize=12)
            plt.ylabel('C-index Promedio', fontsize=12)
            plt.grid(True, alpha=0.3)
            plt.tight_layout()
            
            # Guardar
            output_path = self.output_dir / "plots" / "scatter_cindex_loss.png"
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close()
            
            plot_paths['scatter'] = str(output_path)
            logger.info(f"Scatter plot guardado en {output_path}")
            
        except Exception as e:
            logger.error(f"Error creando scatter plot: {e}")
    
       
    def perform_statistical_tests(self) -> Dict[str, Dict[str, Any]]:
        """
        Realiza tests estadísticos entre modelos.
        
        Returns:
            Diccionario con resultados de tests estadísticos
        """
        logger.info("Realizando tests estadísticos...")
        
        test_results = {}
        
        # Preparar datos para tests
        models_with_data = {}
        for model_name, metrics in self.models_metrics.items():
            if len(metrics.fold_c_indices) > 1:
                models_with_data[model_name] = metrics.fold_c_indices
        
        # Agregar benchmarks
        for name, benchmark in self.external_benchmarks.items():
            if len(benchmark['c_indices']) > 1:
                models_with_data[name] = benchmark['c_indices']
        
        if len(models_with_data) < 2:
            logger.warning("Insuficientes modelos para tests estadísticos")
            return test_results
        
        model_names = list(models_with_data.keys())
        
        # 1. Test de normalidad (Shapiro-Wilk)
        normality_tests = {}
        for model_name, data in models_with_data.items():
            if len(data) >= 3:  # Mínimo para Shapiro-Wilk
                stat, p_value = stats.shapiro(data)
                normality_tests[model_name] = {
                    'statistic': stat,
                    'p_value': p_value,
                    'is_normal': p_value > 0.05
                }
        
        test_results['normality'] = normality_tests
        
        # 2. Tests pareados (t-test o Mann-Whitney)
        pairwise_tests = {}
        for i in range(len(model_names)):
            for j in range(i + 1, len(model_names)):
                model1, model2 = model_names[i], model_names[j]
                data1, data2 = models_with_data[model1], models_with_data[model2]
                
                # Determinar si usar test paramétrico o no paramétrico
                is_normal1 = normality_tests.get(model1, {}).get('is_normal', False)
                is_normal2 = normality_tests.get(model2, {}).get('is_normal', False)
                
                if is_normal1 and is_normal2:
                    # Test t de Student
                    stat, p_value = stats.ttest_ind(data1, data2)
                    test_name = 't-test'
                else:
                    # Test Mann-Whitney U
                    stat, p_value = stats.mannwhitneyu(data1, data2, alternative='two-sided')
                    test_name = 'Mann-Whitney U'
                
                pairwise_tests[f"{model1}_vs_{model2}"] = {
                    'test_name': test_name,
                    'statistic': stat,
                    'p_value': p_value,
                    'significant': p_value < 0.05,
                    'effect_size': abs(np.mean(data1) - np.mean(data2)) / np.sqrt((np.var(data1) + np.var(data2)) / 2)
                }
        
        test_results['pairwise'] = pairwise_tests
        
        # 3. ANOVA (si hay más de 2 modelos)
        if len(models_with_data) > 2:
            all_data = list(models_with_data.values())
            f_stat, p_value = stats.f_oneway(*all_data)
            
            test_results['anova'] = {
                'f_statistic': f_stat,
                'p_value': p_value,
                'significant': p_value < 0.05,
                'models_compared': list(models_with_data.keys())
            }
        
        # Guardar resultados
        output_path = self.output_dir / "data" / "statistical_tests.json"
        with open(output_path, 'w') as f:
            json.dump(test_results, f, indent=2, default=str)
        
        logger.info(f"Resultados de tests estadísticos guardados en {output_path}")
        return test_results
    
    def generate_report(self, include_plots: bool = True) -> str:
        """
        Genera reporte completo en formato texto.
        
        Args:
            include_plots: Si incluir referencias a plots
        
        Returns:
            Ruta del archivo de reporte generado
        """
        logger.info("Generando reporte completo...")
        
        report_lines = []
        
        # Encabezado
        report_lines.extend([
            "=" * 80,
            "REPORTE DE EVALUACIÓN DE MODELOS GNN - CROSS-VALIDATION",
            "=" * 80,
            f"Fecha: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}",
            f"Directorio de resultados: {self.results_dir}",
            f"Directorio de salida: {self.output_dir}",
            "",
            "RESUMEN EJECUTIVO",
            "-" * 20,
        ])
        
        # Información general
        n_models = len(self.models_metrics)
        n_benchmarks = len(self.external_benchmarks)
        
        report_lines.append(f"• Modelos GNN evaluados: {n_models}")
        report_lines.append(f"• Benchmarks externos: {n_benchmarks}")
        report_lines.append("")
        
        # Mejor modelo
        if self.models_metrics:
            best_model = max(self.models_metrics.items(), key=lambda x: x[1].mean_c_index)
            report_lines.extend([
                f"• Mejor modelo GNN: {best_model[0]}",
                f"  - C-index promedio: {best_model[1].mean_c_index:.4f} ± {best_model[1].std_c_index:.4f}",
                f"  - Rango C-index: [{best_model[1].min_c_index:.4f}, {best_model[1].max_c_index:.4f}]",
                f"  - Número de folds: {best_model[1].n_folds}",
                ""
            ])
        
        # Tabla de resultados
        report_lines.extend([
            "RESULTADOS DETALLADOS",
            "-" * 25,
            ""
        ])
        
        # Modelos GNN
        if self.models_metrics:
            report_lines.append("Modelos GNN:")
            report_lines.append("-" * 15)
            
            for model_name, metrics in sorted(self.models_metrics.items(), 
                                            key=lambda x: x[1].mean_c_index, reverse=True):
                report_lines.extend([
                    f"• {model_name}:",
                    f"  - C-index promedio: {metrics.mean_c_index:.4f} ± {metrics.std_c_index:.4f}",
                    f"  - C-index mediano: {metrics.median_c_index:.4f}",
                    f"  - Rango C-index: [{metrics.min_c_index:.4f}, {metrics.max_c_index:.4f}]",
                    f"  - IQR C-index: [{metrics.q25_c_index:.4f}, {metrics.q75_c_index:.4f}]",
                    f"  - Loss promedio: {metrics.mean_loss:.4f} ± {metrics.std_loss:.4f}" if not pd.isna(metrics.mean_loss) else "  - Loss: N/A",
                    f"  - Número de folds: {metrics.n_folds}",
                    f"  - Configuración: {metrics.model_config}",
                    ""
                ])
        
        # Benchmarks externos
        if self.external_benchmarks:
            report_lines.append("Benchmarks Externos:")
            report_lines.append("-" * 20)
            
            for name, benchmark in sorted(self.external_benchmarks.items(), 
                                        key=lambda x: x[1]['mean_c_index'], reverse=True):
                report_lines.extend([
                    f"• {name}:",
                    f"  - C-index promedio: {benchmark['mean_c_index']:.4f} ± {benchmark['std_c_index']:.4f}",
                    f"  - C-index mediano: {benchmark['median_c_index']:.4f}",
                    f"  - Rango C-index: [{benchmark['min_c_index']:.4f}, {benchmark['max_c_index']:.4f}]",
                    f"  - IQR C-index: [{benchmark['q25_c_index']:.4f}, {benchmark['q75_c_index']:.4f}]",
                    f"  - Número de folds: {benchmark['n_folds']}",
                    f"  - Archivo fuente: {benchmark['source_file']}",
                    ""
                ])
        
        # Tests estadísticos
        if hasattr(self, 'statistical_tests'):
            report_lines.extend([
                "TESTS ESTADÍSTICOS",
                "-" * 20,
                ""
            ])
            
            # Tests de normalidad
            if 'normality' in self.statistical_tests:
                report_lines.append("Tests de Normalidad (Shapiro-Wilk):")
                for model, test in self.statistical_tests['normality'].items():
                    status = "Normal" if test['is_normal'] else "No Normal"
                    report_lines.append(f"• {model}: {status} (p={test['p_value']:.4f})")
                report_lines.append("")
            
            # Tests pareados
            if 'pairwise' in self.statistical_tests:
                report_lines.append("Tests Pareados:")
                for comparison, test in self.statistical_tests['pairwise'].items():
                    status = "Significativo" if test['significant'] else "No Significativo"
                    report_lines.append(f"• {comparison}: {status} ({test['test_name']}, p={test['p_value']:.4f})")
                report_lines.append("")
            
            # ANOVA
            if 'anova' in self.statistical_tests:
                anova = self.statistical_tests['anova']
                status = "Significativo" if anova['significant'] else "No Significativo"
                report_lines.extend([
                    f"ANOVA: {status} (F={anova['f_statistic']:.4f}, p={anova['p_value']:.4f})",
                    f"Modelos comparados: {', '.join(anova['models_compared'])}",
                    ""
                ])
        
        # Referencias a archivos generados
        report_lines.extend([
            "ARCHIVOS GENERADOS",
            "-" * 18,
            "",
            "Datos:",
            f"• Tabla resumen: {self.output_dir / 'data' / 'summary_table.csv'}",
            f"• Tests estadísticos: {self.output_dir / 'data' / 'statistical_tests.json'}",
            ""
        ])
        
        if include_plots:
            report_lines.extend([
                "Visualizaciones:",
                f"• Boxplot: {self.output_dir / 'plots' / 'boxplot_comparison.png'}",
                f"• Barplot: {self.output_dir / 'plots' / 'barplot_comparison.png'}",
                f"• Violin plot: {self.output_dir / 'plots' / 'violin_plot.png'}",
                f"• Scatter plot: {self.output_dir / 'plots' / 'scatter_cindex_loss.png'}",
                f"• Heatmap: {self.output_dir / 'plots' / 'heatmap_metrics.png'}",
                f"• Distribuciones: {self.output_dir / 'plots' / 'distribution_cindex.png'}",
                ""
            ])
        
        # Recomendaciones
        report_lines.extend([
            "RECOMENDACIONES",
            "-" * 15,
            "",
            "Basado en los resultados del análisis:",
            ""
        ])
        
        if self.models_metrics:
            best_model = max(self.models_metrics.items(), key=lambda x: x[1].mean_c_index)
            report_lines.extend([
                f"1. El modelo {best_model[0]} muestra el mejor rendimiento promedio",
                f"   con C-index = {best_model[1].mean_c_index:.4f}",
                ""
            ])
            
            # Variabilidad
            most_stable = min(self.models_metrics.items(), key=lambda x: x[1].std_c_index)
            report_lines.extend([
                f"2. El modelo {most_stable[0]} muestra la mayor estabilidad",
                f"   con desviación estándar = {most_stable[1].std_c_index:.4f}",
                ""
            ])
        
        report_lines.extend([
            "3. Considere validar los mejores modelos en un conjunto de test independiente",
            "4. Evalúe la significancia estadística de las diferencias entre modelos",
            "5. Considere el balance entre rendimiento y estabilidad para la selección final",
            "",
            "=" * 80
        ])
        
        # Guardar reporte
        output_path = self.output_dir / "reports" / "evaluation_report.txt"
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write('\n'.join(report_lines))
        
        logger.info(f"Reporte completo guardado en {output_path}")
        return str(output_path)
    
    def run_complete_evaluation(self, model_patterns: List[str] = None, 
                              external_benchmarks: Dict[str, str] = None) -> Dict[str, Any]:
        """
        Ejecuta evaluación completa de modelos.
        
        Args:
            model_patterns: Patrones para filtrar modelos
            external_benchmarks: Diccionario con benchmarks externos
        
        Returns:
            Diccionario con resumen de resultados
        """
        logger.info("Iniciando evaluación completa...")
        
        results = {}
        
        try:
            # 1. Cargar datos
            self.load_cv_results(model_patterns)
            results['models_loaded'] = len(self.models_data)
            
            # 2. Cargar benchmarks externos
            if external_benchmarks:
                self.load_external_benchmarks(external_benchmarks)
                results['benchmarks_loaded'] = len(self.external_benchmarks)
            
            # 3. Calcular métricas
            self.calculate_metrics()
            results['metrics_calculated'] = len(self.models_metrics)
            
            # 4. Generar tabla resumen
            summary_table = self.generate_summary_table()
            results['summary_table_shape'] = summary_table.shape
            
            # 5. Crear visualizaciones
            plot_paths = self.create_comparison_plots()
            results['plots_created'] = len(plot_paths)
            
            # 6. Tests estadísticos
            self.statistical_tests = self.perform_statistical_tests()
            results['statistical_tests_performed'] = len(self.statistical_tests)
            
            # 7. Generar reporte
            report_path = self.generate_report()
            results['report_generated'] = report_path
            
            # 8. Guardar configuración completa
            config_data = {
                'evaluation_config': {
                    'results_dir': str(self.results_dir),
                    'output_dir': str(self.output_dir),
                    'model_patterns': model_patterns,
                    'external_benchmarks': external_benchmarks,
                    'execution_time': pd.Timestamp.now().isoformat()
                },
                'results_summary': results
            }
            
            config_path = self.output_dir / "data" / "evaluation_config.json"
            with open(config_path, 'w') as f:
                json.dump(config_data, f, indent=2, default=str)
            
            logger.info("Evaluación completa terminada exitosamente")
            return results
            
        except Exception as e:
            logger.error(f"Error durante la evaluación completa: {e}")
            raise

def main():
    """Función principal del script."""
    parser = argparse.ArgumentParser(
        description="Evalúa y compara resultados de cross-validation de modelos GNN",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Ejemplos de uso:

1. Evaluación básica:
   python evaluate_cv_results.py --results_dir results/gnn/basic_model

2. Evaluación con filtros:
   python evaluate_cv_results.py --results_dir results/gnn/basic_model --patterns v1_ v2_

3. Evaluación con benchmarks externos:
   python evaluate_cv_results.py --results_dir results/gnn/basic_model \\
       --external_benchmarks IPSSM:results/gnn/IPSSM_cv_folds/cv_summary.tsv

4. Evaluación completa con salida personalizada:
   python evaluate_cv_results.py --results_dir results/gnn/basic_model \\
       --output_dir results/gnn/evaluation_2025 \\
       --patterns v1_ v2_ v3_ \\
       --external_benchmarks IPSSM:results/gnn/IPSSM_cv_folds/cv_summary.tsv \\
       --verbose
        """
    )
    
    parser.add_argument(
        '--results_dir',
        type=str,
        required=True,
        help='Directorio con resultados de CV (archivos *_cv_results.pkl y *_cv_summary.csv)'
    )
    
    parser.add_argument(
        '--output_dir',
        type=str,
        default=None,
        help='Directorio de salida (por defecto: {results_dir}/evaluation)'
    )
    
    parser.add_argument(
        '--patterns',
        type=str,
        nargs='*',
        default=None,
        help='Patrones para filtrar modelos (ej: v1_ v2_boolean)'
    )
    
    parser.add_argument(
        '--external_benchmarks',
        type=str,
        nargs='*',
        default=None,
        help='Benchmarks externos en formato nombre:ruta (ej: IPSSM:path/to/file.csv)'
    )
    
    parser.add_argument(
        '--verbose',
        action='store_true',
        help='Activar logging detallado'
    )
    
    parser.add_argument(
        '--no_plots',
        action='store_true',
        help='No generar visualizaciones (solo datos y reportes)'
    )
    
    args = parser.parse_args()
    
    # Configurar logging
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Validar argumentos
    if not Path(args.results_dir).exists():
        logger.error(f"Directorio de resultados no encontrado: {args.results_dir}")
        sys.exit(1)
    
    # Procesar benchmarks externos
    external_benchmarks = {}
    if args.external_benchmarks:
        for benchmark in args.external_benchmarks:
            if ':' not in benchmark:
                logger.error(f"Formato inválido para benchmark: {benchmark}")
                logger.error("Use formato: nombre:ruta")
                sys.exit(1)
            
            name, path = benchmark.split(':', 1)
            external_benchmarks[name] = path
    
    try:
        # Crear evaluador
        evaluator = CVResultsEvaluator(
            results_dir=args.results_dir,
            output_dir=args.output_dir
        )
        
        # Ejecutar evaluación completa
        results = evaluator.run_complete_evaluation(
            model_patterns=args.patterns,
            external_benchmarks=external_benchmarks if external_benchmarks else None
        )
        
        # Mostrar resumen
        logger.info("\n" + "="*60)
        logger.info("RESUMEN DE EVALUACIÓN")
        logger.info("="*60)
        logger.info(f"Modelos cargados: {results['models_loaded']}")
        if 'benchmarks_loaded' in results:
            logger.info(f"Benchmarks cargados: {results['benchmarks_loaded']}")
        logger.info(f"Métricas calculadas: {results['metrics_calculated']}")
        logger.info(f"Visualizaciones creadas: {results['plots_created']}")
        logger.info(f"Tests estadísticos: {results['statistical_tests_performed']}")
        logger.info(f"Reporte generado: {results['report_generated']}")
        logger.info(f"Directorio de salida: {evaluator.output_dir}")
        logger.info("="*60)
        
        # Mostrar mejores modelos
        if evaluator.models_metrics:
            logger.info("\nMEJORES MODELOS (por C-index):")
            sorted_models = sorted(evaluator.models_metrics.items(), 
                                 key=lambda x: x[1].mean_c_index, reverse=True)
            for i, (name, metrics) in enumerate(sorted_models[:5]):
                logger.info(f"{i+1}. {name}: {metrics.mean_c_index:.4f} ± {metrics.std_c_index:.4f}")
        
        logger.info("\nEvaluación completada exitosamente!")
        
    except Exception as e:
        logger.error(f"Error durante la evaluación: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()