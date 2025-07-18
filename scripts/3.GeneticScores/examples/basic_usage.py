#!/usr/bin/env python3
"""
Ejemplo básico de uso del pipeline unificado MDS.

Este ejemplo demuestra:
1. Configuración simple de variables
2. Ejecución del pipeline
3. Validación básica de resultados
"""

import yaml
import subprocess
from pathlib import Path
import sys
import logging

# Configurar logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def create_basic_config():
    """Crea una configuración básica para análisis MDS."""
    
    config = {
        'data_sources': {
            'clinical_data': 'data/IPSSMol/df_clinical.tsv',
            'mutation_vaf_data': 'data/IPSSMol/maf.tsv',
            'gene_embeddings': 'data/scGPT/gene_embeddings.csv'
        },
        
        'variable_processing': {
            # Selección de genes más importantes en MDS
            'gene_selection': [
                'TP53', 'ASXL1', 'RUNX1', 'SRSF2', 'SF3B1',
                'TET2', 'DNMT3A', 'IDH2', 'CBL', 'EZH2'
            ],
            
            # Variables clínicas continuas
            'continuous_variables': [
                'AGE',          # Edad
                'HB',           # Hemoglobina
                'PLT',          # Plaquetas
                'ANC',          # Neutrófilos
                'BM_BLAST'      # Blastos en médula ósea
            ],
            
            # Variables categóricas
            'categorical_variables': [
                'SEX'           # Sexo (M/F)
            ],
            
            # Variables binarias
            'binary_variables': [
                'AML_TRANSF'    # Transformación a AML
            ],
            
            # Variables de supervivencia
            'survival_variables': {
                'time_variable': 'LFS_YEARS',    # Tiempo libre de leucemia
                'status_variable': 'LFS_STATUS'  # Estado (0=censurado, 1=evento)
            }
        },
        
        'cross_validation': {
            'fold_column': 'fold'  # Columna con información de folds
        },
        
        'data_generation': {
            'output_formats': ['table', 'graph'],  # Generar ambos formatos
            
            'table_format': {
                'output_filename': 'basic_table'
            },
            
            'graph_format': {
                'output_filename': 'basic_graph'
            }
        },
        
        'output': {
            'file_prefix': 'basic_example'  # Prefijo para archivos de salida
        }
    }
    
    return config

def run_basic_pipeline():
    """Ejecuta el pipeline con configuración básica."""
    
    logger.info("🚀 Iniciando ejemplo básico del pipeline MDS")
    
    # 1. Crear configuración
    logger.info("📋 Creando configuración básica...")
    config = create_basic_config()
    
    # Guardar configuración
    config_file = Path('configs/basic_example.yaml')
    config_file.parent.mkdir(exist_ok=True)
    
    with open(config_file, 'w', encoding='utf-8') as f:
        yaml.dump(config, f, default_flow_style=False, allow_unicode=True)
    
    logger.info(f"✓ Configuración guardada en: {config_file}")
    
    # 2. Verificar archivos de entrada
    logger.info("📂 Verificando archivos de entrada...")
    base_path = Path('../../')  # Ajustar según la estructura
    
    required_files = [
        base_path / config['data_sources']['clinical_data'],
        base_path / config['data_sources']['mutation_vaf_data']
    ]
    
    missing_files = []
    for file_path in required_files:
        if not file_path.exists():
            missing_files.append(str(file_path))
        else:
            logger.info(f"✓ Archivo encontrado: {file_path}")
    
    if missing_files:
        logger.error(f"❌ Archivos faltantes: {missing_files}")
        logger.info("💡 Asegúrate de que los datos estén en la ubicación correcta")
        return False
    
    # 3. Ejecutar pipeline
    logger.info("⚙️ Ejecutando pipeline...")
    output_dir = 'results/basic_example'
    
    cmd = [
        sys.executable, 'prepare_data_unified.py',
        '--config', str(config_file),
        '--output', output_dir,
        '--verbose'
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        
        if result.returncode == 0:
            logger.info("✅ Pipeline ejecutado exitosamente")
            logger.info(f"📁 Resultados disponibles en: {output_dir}")
        else:
            logger.error("❌ Error ejecutando pipeline")
            logger.error(f"STDERR: {result.stderr}")
            return False
            
    except subprocess.TimeoutExpired:
        logger.error("❌ Pipeline timeout (>5 min)")
        return False
    except Exception as e:
        logger.error(f"❌ Error inesperado: {e}")
        return False
    
    # 4. Verificar resultados
    logger.info("🔍 Verificando resultados...")
    output_path = Path(output_dir)
    
    expected_files = [
        output_path / 'basic_example_table.pkl',
        output_path / 'basic_example_graph.pt'
    ]
    
    for expected_file in expected_files:
        if expected_file.exists():
            size_mb = expected_file.stat().st_size / (1024 * 1024)
            logger.info(f"✓ Generado: {expected_file.name} ({size_mb:.2f} MB)")
        else:
            logger.warning(f"⚠ No encontrado: {expected_file.name}")
    
    # 5. Mostrar resumen
    logger.info("\n" + "="*60)
    logger.info("📊 RESUMEN DEL EJEMPLO BÁSICO")
    logger.info("="*60)
    logger.info(f"Genes seleccionados: {len(config['variable_processing']['gene_selection'])}")
    logger.info(f"Variables clínicas: {len(config['variable_processing']['continuous_variables']) + len(config['variable_processing']['categorical_variables']) + len(config['variable_processing']['binary_variables'])}")
    logger.info(f"Formatos generados: {', '.join(config['data_generation']['output_formats'])}")
    logger.info(f"Directorio de salida: {output_dir}")
    logger.info("="*60)
    
    return True

def validate_results():
    """Ejecuta validación básica de los resultados."""
    
    logger.info("🔬 Ejecutando validación básica...")
    
    output_dir = 'results/basic_example'
    
    # Validación de calidad
    cmd_quality = [
        sys.executable, 'tests/validate_step3c.py',
        '--data-dir', output_dir
    ]
    
    try:
        result = subprocess.run(cmd_quality, capture_output=True, text=True, timeout=60)
        
        if result.returncode == 0:
            logger.info("✅ Validación de calidad exitosa")
        else:
            logger.warning("⚠ Problemas en validación de calidad")
            logger.warning(f"Detalles: {result.stderr}")
    
    except Exception as e:
        logger.warning(f"⚠ Error en validación: {e}")

def main():
    """Función principal del ejemplo básico."""
    
    # Cambiar al directorio correcto
    script_dir = Path(__file__).parent.parent
    original_dir = Path.cwd()
    
    try:
        import os
        os.chdir(script_dir)
        
        logger.info(f"📁 Directorio de trabajo: {Path.cwd()}")
        
        # Ejecutar pipeline
        success = run_basic_pipeline()
        
        if success:
            # Validar resultados
            validate_results()
            
            logger.info("\n🎉 ¡Ejemplo básico completado exitosamente!")
            logger.info("💡 Revisa los archivos generados en results/basic_example/")
            logger.info("📖 Para más ejemplos, consulta el directorio examples/")
            
            return 0
        else:
            logger.error("\n💥 El ejemplo básico falló")
            logger.info("🔧 Revisa los logs anteriores para diagnosticar el problema")
            return 1
    
    finally:
        os.chdir(original_dir)

if __name__ == "__main__":
    exit(main())
