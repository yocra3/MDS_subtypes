#!/usr/bin/env python3
"""
Ejemplo avanzado de uso del pipeline unificado MDS.

Este ejemplo demuestra:
1. Configuración compleja con múltiples experimentos
2. Análisis comparativo de diferentes configuraciones
3. Validación exhaustiva y generación de reportes
4. Uso de embeddings personalizados
"""

import yaml
import subprocess
import pandas as pd
import numpy as np
from pathlib import Path
import sys
import logging
from typing import Dict, List, Any
import json
from datetime import datetime

# Configurar logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class AdvancedMDSAnalysis:
    """Clase para análisis avanzado con múltiples configuraciones."""
    
    def __init__(self, base_output_dir: str = "results/advanced_analysis"):
        """Inicializa el análisis avanzado."""
        self.base_output_dir = Path(base_output_dir)
        self.base_output_dir.mkdir(parents=True, exist_ok=True)
        
        self.experiments = {}
        self.results = {}
        
        logger.info(f"🔬 Inicializando análisis avanzado MDS")
        logger.info(f"📁 Directorio base: {self.base_output_dir}")
    
    def create_experiment_configs(self) -> Dict[str, Dict[str, Any]]:
        """Crea múltiples configuraciones experimentales."""
        
        logger.info("⚙️ Creando configuraciones experimentales...")
        
        # Configuración base común
        base_config = {
            'data_sources': {
                'clinical_data': 'data/IPSSMol/df_clinical.tsv',
                'mutation_vaf_data': 'data/IPSSMol/maf.tsv',
                'gene_embeddings': 'data/scGPT/gene_embeddings.csv'
            },
            'cross_validation': {
                'fold_column': 'fold'
            },
            'data_generation': {
                'output_formats': ['table', 'graph']
            }
        }
        
        # Experimento 1: Genes de alto impacto (top 10)
        exp1_config = base_config.copy()
        exp1_config.update({
            'variable_processing': {
                'gene_selection': [
                    'TP53', 'ASXL1', 'RUNX1', 'SRSF2', 'SF3B1',
                    'TET2', 'DNMT3A', 'IDH2', 'CBL', 'EZH2'
                ],
                'continuous_variables': ['AGE', 'HB', 'PLT', 'ANC', 'BM_BLAST'],
                'categorical_variables': ['SEX'],
                'binary_variables': ['AML_TRANSF'],
                'survival_variables': {
                    'time_variable': 'LFS_YEARS',
                    'status_variable': 'LFS_STATUS'
                }
            },
            'output': {'file_prefix': 'top_genes'}
        })
        exp1_config['data_generation'].update({
            'table_format': {'output_filename': 'top_genes_table'},
            'graph_format': {'output_filename': 'top_genes_graph'}
        })
        
        # Experimento 2: Genes extendidos (top 20)
        exp2_config = base_config.copy()
        exp2_config.update({
            'variable_processing': {
                'gene_selection': [
                    'TP53', 'ASXL1', 'RUNX1', 'SRSF2', 'SF3B1',
                    'TET2', 'DNMT3A', 'IDH2', 'CBL', 'EZH2',
                    'U2AF1', 'ZRSR2', 'STAG2', 'RAD21', 'SMC3',
                    'BCOR', 'BCORL1', 'PHF6', 'WT1', 'NPM1'
                ],
                'continuous_variables': ['AGE', 'HB', 'PLT', 'ANC', 'BM_BLAST', 'LDH'],
                'categorical_variables': ['SEX', 'WHO_SUBTYPE'],
                'binary_variables': ['AML_TRANSF', 'PRIOR_MPN'],
                'survival_variables': {
                    'time_variable': 'LFS_YEARS',
                    'status_variable': 'LFS_STATUS'
                }
            },
            'output': {'file_prefix': 'extended_genes'}
        })
        exp2_config['data_generation'].update({
            'table_format': {'output_filename': 'extended_genes_table'},
            'graph_format': {'output_filename': 'extended_genes_graph'}
        })
        
        # Experimento 3: Solo variables clínicas (sin genes)
        exp3_config = base_config.copy()
        exp3_config.update({
            'variable_processing': {
                'gene_selection': [],  # Sin genes
                'continuous_variables': [
                    'AGE', 'HB', 'PLT', 'ANC', 'BM_BLAST', 'LDH', 
                    'FERRITIN', 'B2_MICRO'
                ],
                'categorical_variables': ['SEX', 'WHO_SUBTYPE', 'CYTOGENETICS'],
                'binary_variables': ['AML_TRANSF', 'PRIOR_MPN', 'PRIOR_CHEMO'],
                'survival_variables': {
                    'time_variable': 'LFS_YEARS',
                    'status_variable': 'LFS_STATUS'
                }
            },
            'output': {'file_prefix': 'clinical_only'}
        })
        exp3_config['data_generation'].update({
            'output_formats': ['table'],  # Solo tabla para análisis clínico
            'table_format': {'output_filename': 'clinical_only_table'}
        })
        
        experiments = {
            'top_genes': exp1_config,
            'extended_genes': exp2_config,
            'clinical_only': exp3_config
        }
        
        logger.info(f"✓ Creadas {len(experiments)} configuraciones experimentales")
        return experiments
    
    def run_experiment(self, exp_name: str, config: Dict[str, Any]) -> Dict[str, Any]:
        """Ejecuta un experimento individual."""
        
        logger.info(f"🧪 Ejecutando experimento: {exp_name}")
        
        # Crear directorio específico
        exp_dir = self.base_output_dir / exp_name
        exp_dir.mkdir(exist_ok=True)
        
        # Guardar configuración
        config_file = exp_dir / f"{exp_name}_config.yaml"
        with open(config_file, 'w', encoding='utf-8') as f:
            yaml.dump(config, f, default_flow_style=False, allow_unicode=True)
        
        # Ejecutar pipeline
        cmd = [
            sys.executable, 'prepare_data_unified.py',
            '--config', str(config_file),
            '--output', str(exp_dir),
            '--verbose'
        ]
        
        try:
            start_time = datetime.now()
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
            end_time = datetime.now()
            
            execution_time = (end_time - start_time).total_seconds()
            
            if result.returncode == 0:
                logger.info(f"✅ Experimento {exp_name} completado en {execution_time:.1f}s")
                
                # Validar resultados
                validation_results = self._validate_experiment(exp_name, exp_dir)
                
                return {
                    'status': 'success',
                    'execution_time': execution_time,
                    'output_dir': str(exp_dir),
                    'validation': validation_results,
                    'stdout': result.stdout,
                    'stderr': result.stderr
                }
            else:
                logger.error(f"❌ Experimento {exp_name} falló")
                logger.error(f"Error: {result.stderr}")
                
                return {
                    'status': 'failed',
                    'execution_time': execution_time,
                    'error': result.stderr,
                    'stdout': result.stdout
                }
        
        except subprocess.TimeoutExpired:
            logger.error(f"❌ Experimento {exp_name} timeout (>10 min)")
            return {'status': 'timeout'}
        
        except Exception as e:
            logger.error(f"❌ Error inesperado en {exp_name}: {e}")
            return {'status': 'error', 'error': str(e)}
    
    def _validate_experiment(self, exp_name: str, exp_dir: Path) -> Dict[str, Any]:
        """Valida los resultados de un experimento."""
        
        logger.info(f"🔍 Validando experimento {exp_name}...")
        
        validation_results = {}
        
        # 1. Validación de calidad
        cmd_quality = [
            sys.executable, 'tests/validate_step3c.py',
            '--data-dir', str(exp_dir)
        ]
        
        try:
            result = subprocess.run(cmd_quality, capture_output=True, text=True, timeout=120)
            validation_results['quality'] = {
                'status': 'success' if result.returncode == 0 else 'failed',
                'output': result.stdout,
                'errors': result.stderr
            }
        except Exception as e:
            validation_results['quality'] = {'status': 'error', 'error': str(e)}
        
        # 2. Validación de compatibilidad (solo si hay archivos de grafo)
        if any(f.suffix == '.pt' for f in exp_dir.glob('*.pt')):
            cmd_compat = [
                sys.executable, 'tests/validate_step3d.py',
                '--data-dir', str(exp_dir)
            ]
            
            try:
                result = subprocess.run(cmd_compat, capture_output=True, text=True, timeout=120)
                validation_results['compatibility'] = {
                    'status': 'success' if result.returncode == 0 else 'failed',
                    'output': result.stdout,
                    'errors': result.stderr
                }
            except Exception as e:
                validation_results['compatibility'] = {'status': 'error', 'error': str(e)}
        
        return validation_results
    
    def analyze_file_sizes(self) -> Dict[str, Any]:
        """Analiza los tamaños de archivos generados."""
        
        logger.info("📊 Analizando tamaños de archivos...")
        
        file_analysis = {}
        
        for exp_name in self.results:
            if self.results[exp_name]['status'] == 'success':
                exp_dir = Path(self.results[exp_name]['output_dir'])
                
                files_info = {}
                for file_path in exp_dir.glob('*'):
                    if file_path.is_file() and file_path.suffix in ['.pkl', '.pt', '.csv', '.txt']:
                        size_mb = file_path.stat().st_size / (1024 * 1024)
                        files_info[file_path.name] = {
                            'size_mb': round(size_mb, 2),
                            'type': file_path.suffix
                        }
                
                file_analysis[exp_name] = files_info
        
        return file_analysis
    
    def generate_comparison_report(self) -> str:
        """Genera un reporte comparativo de todos los experimentos."""
        
        logger.info("📋 Generando reporte comparativo...")
        
        report_lines = [
            "="*80,
            "REPORTE COMPARATIVO - ANÁLISIS AVANZADO MDS",
            "="*80,
            f"Fecha: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
            f"Experimentos ejecutados: {len(self.results)}",
            ""
        ]
        
        # Resumen de experimentos
        report_lines.append("📊 RESUMEN DE EXPERIMENTOS:")
        for exp_name, result in self.results.items():
            status_icon = "✅" if result['status'] == 'success' else "❌"
            exec_time = result.get('execution_time', 0)
            
            report_lines.append(f"  {status_icon} {exp_name}: {result['status']} ({exec_time:.1f}s)")
        
        report_lines.append("")
        
        # Análisis de archivos
        file_analysis = self.analyze_file_sizes()
        
        if file_analysis:
            report_lines.append("📁 TAMAÑOS DE ARCHIVOS GENERADOS:")
            for exp_name, files in file_analysis.items():
                report_lines.append(f"  📂 {exp_name}:")
                for filename, info in files.items():
                    report_lines.append(f"    - {filename}: {info['size_mb']} MB")
            
            report_lines.append("")
        
        # Validaciones
        report_lines.append("🔬 RESULTADOS DE VALIDACIÓN:")
        for exp_name, result in self.results.items():
            if result['status'] == 'success' and 'validation' in result:
                validation = result['validation']
                
                report_lines.append(f"  📋 {exp_name}:")
                
                if 'quality' in validation:
                    quality_status = validation['quality']['status']
                    quality_icon = "✅" if quality_status == 'success' else "❌"
                    report_lines.append(f"    {quality_icon} Calidad: {quality_status}")
                
                if 'compatibility' in validation:
                    compat_status = validation['compatibility']['status']
                    compat_icon = "✅" if compat_status == 'success' else "❌"
                    report_lines.append(f"    {compat_icon} Compatibilidad: {compat_status}")
        
        report_lines.extend([
            "",
            "="*80,
            "RECOMENDACIONES:",
            ""
        ])
        
        # Generar recomendaciones
        successful_experiments = [name for name, result in self.results.items() 
                                if result['status'] == 'success']
        
        if successful_experiments:
            report_lines.append(f"✅ {len(successful_experiments)} experimentos exitosos")
            report_lines.append("💡 Revisar archivos generados para análisis downstream")
            
            if len(successful_experiments) > 1:
                report_lines.append("🔄 Considerar análisis comparativo de performance")
        else:
            report_lines.append("❌ No hay experimentos exitosos")
            report_lines.append("🔧 Revisar configuraciones y datos de entrada")
        
        report_lines.extend([
            "",
            "="*80
        ])
        
        return "\n".join(report_lines)
    
    def run_all_experiments(self) -> bool:
        """Ejecuta todos los experimentos configurados."""
        
        logger.info("🚀 Iniciando análisis avanzado con múltiples experimentos...")
        
        # Crear configuraciones
        self.experiments = self.create_experiment_configs()
        
        # Ejecutar cada experimento
        for exp_name, config in self.experiments.items():
            self.results[exp_name] = self.run_experiment(exp_name, config)
        
        # Generar reporte
        report = self.generate_comparison_report()
        
        # Guardar reporte
        report_file = self.base_output_dir / 'comparative_analysis_report.txt'
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write(report)
        
        # Guardar resultados detallados
        results_file = self.base_output_dir / 'detailed_results.json'
        with open(results_file, 'w', encoding='utf-8') as f:
            json.dump(self.results, f, indent=2, default=str)
        
        logger.info(f"📋 Reporte guardado en: {report_file}")
        logger.info(f"📊 Resultados detallados en: {results_file}")
        
        # Mostrar reporte
        print("\n" + report)
        
        # Verificar éxito general
        successful_count = sum(1 for result in self.results.values() 
                             if result['status'] == 'success')
        
        logger.info(f"🎯 Experimentos exitosos: {successful_count}/{len(self.results)}")
        
        return successful_count > 0

def main():
    """Función principal del ejemplo avanzado."""
    
    # Cambiar al directorio correcto
    script_dir = Path(__file__).parent.parent
    original_dir = Path.cwd()
    
    try:
        import os
        os.chdir(script_dir)
        
        logger.info(f"📁 Directorio de trabajo: {Path.cwd()}")
        
        # Crear y ejecutar análisis
        analysis = AdvancedMDSAnalysis()
        success = analysis.run_all_experiments()
        
        if success:
            logger.info("\n🎉 ¡Análisis avanzado completado exitosamente!")
            logger.info(f"📁 Resultados disponibles en: {analysis.base_output_dir}")
            logger.info("📊 Revisa el reporte comparativo para análisis detallado")
            return 0
        else:
            logger.error("\n💥 El análisis avanzado tuvo problemas")
            logger.info("🔧 Revisa los logs y reportes para más detalles")
            return 1
    
    finally:
        os.chdir(original_dir)

if __name__ == "__main__":
    exit(main())
