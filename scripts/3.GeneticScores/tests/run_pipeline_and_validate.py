#!/usr/bin/env python3
"""
Script automatizado para ejecutar el pipeline completo y validar la calidad de datos.
Combina la ejecución del pipeline (Paso 3A) con la validación de calidad (Paso 3C).
"""
import sys
import os
from pathlib import Path
import logging
import subprocess
import argparse
import yaml
import json

# Configurar logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def run_pipeline_and_validate(config_file: str, output_dir: str, verbose: bool = False) -> bool:
    """
    Ejecuta el pipeline y valida la calidad de datos.
    
    Args:
        config_file: Archivo de configuración del pipeline
        output_dir: Directorio de salida
        verbose: Activar logging detallado
        
    Returns:
        True si todo fue exitoso
    """
    logger.info("🚀 Iniciando ejecución automatizada del pipeline y validación...")
    
    try:
        # Preparar directorios
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        script_dir = Path(__file__).parent.parent
        pipeline_script = script_dir / "prepare_data_unified.py"
        validator_script = Path(__file__).parent / "validate_step3c.py"
        
        # Paso 1: Ejecutar pipeline
        logger.info("📊 Paso 1: Ejecutando pipeline de preparación de datos...")
        
        pipeline_cmd = [
            sys.executable,
            str(pipeline_script),
            "--config", config_file,
            "--output", output_dir
        ]
        
        if verbose:
            pipeline_cmd.append("--verbose")
        
        logger.info(f"Comando pipeline: {' '.join(pipeline_cmd)}")
        
        result = subprocess.run(
            pipeline_cmd,
            capture_output=True,
            text=True,
            timeout=300  # 5 minutos timeout
        )
        
        if result.returncode != 0:
            logger.error("❌ Pipeline falló")
            logger.error(f"STDOUT: {result.stdout}")
            logger.error(f"STDERR: {result.stderr}")
            return False
        
        logger.info("✅ Pipeline completado exitosamente")
        
        # Paso 2: Validar calidad de datos
        logger.info("🔍 Paso 2: Validando calidad de datos generados...")
        
        validator_cmd = [
            sys.executable,
            str(validator_script),
            "--data-dir", output_dir
        ]
        
        if verbose:
            validator_cmd.append("--verbose")
        
        logger.info(f"Comando validador: {' '.join(validator_cmd)}")
        
        result = subprocess.run(
            validator_cmd,
            capture_output=True,
            text=True,
            timeout=120  # 2 minutos timeout
        )
        
        # Mostrar output del validador
        if result.stdout:
            print("\n" + "="*60)
            print("REPORTE DE VALIDACIÓN DE CALIDAD")
            print("="*60)
            print(result.stdout)
        
        if result.returncode != 0:
            logger.error("❌ Validación de calidad falló")
            if result.stderr:
                logger.error(f"STDERR: {result.stderr}")
            return False
        
        logger.info("✅ Validación de calidad completada exitosamente")
        
        # Paso 3: Resumen final
        logger.info("📋 Paso 3: Generando resumen final...")
        
        summary = generate_execution_summary(config_file, output_dir)
        
        summary_file = output_path / "execution_summary.json"
        with open(summary_file, 'w', encoding='utf-8') as f:
            json.dump(summary, f, indent=2, ensure_ascii=False, default=str)
        
        logger.info(f"📄 Resumen guardado en: {summary_file}")
        
        logger.info("🎉 Ejecución automatizada completada exitosamente!")
        return True
        
    except subprocess.TimeoutExpired:
        logger.error("❌ Timeout en la ejecución")
        return False
    except Exception as e:
        logger.error(f"❌ Error en ejecución automatizada: {e}")
        import traceback
        logger.error(f"Traceback: {traceback.format_exc()}")
        return False


def generate_execution_summary(config_file: str, output_dir: str) -> dict:
    """Genera resumen de la ejecución."""
    summary = {
        'execution_timestamp': str(pd.Timestamp.now()),
        'config_file': config_file,
        'output_directory': output_dir,
        'pipeline_status': 'completed',
        'validation_status': 'completed',
        'generated_files': [],
        'file_sizes': {}
    }
    
    try:
        output_path = Path(output_dir)
        
        # Listar archivos generados
        for file_path in output_path.iterdir():
            if file_path.is_file():
                summary['generated_files'].append(file_path.name)
                summary['file_sizes'][file_path.name] = file_path.stat().st_size
        
        # Cargar configuración usada
        with open(config_file, 'r', encoding='utf-8') as f:
            if config_file.endswith('.yaml') or config_file.endswith('.yml'):
                config_data = yaml.safe_load(f)
            else:
                config_data = json.load(f)
        
        summary['config_summary'] = {
            'output_formats': config_data.get('data_generation', {}).get('output_formats', []),
            'gene_selection': len(config_data.get('variable_processing', {}).get('gene_selection', [])),
            'clinical_variables': len(config_data.get('variable_processing', {}).get('continuous_variables', [])) + 
                                len(config_data.get('variable_processing', {}).get('categorical_variables', [])) +
                                len(config_data.get('variable_processing', {}).get('binary_variables', []))
        }
        
    except Exception as e:
        logger.warning(f"⚠ Error generando resumen: {e}")
        summary['error'] = str(e)
    
    return summary


def main():
    """Función principal."""
    parser = argparse.ArgumentParser(
        description="Ejecutor Automatizado del Pipeline + Validación de Calidad",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Ejemplos de uso:
  python run_pipeline_and_validate.py --config configs/real_data_config.yaml --output results/test_run
  python run_pipeline_and_validate.py --config configs/real_data_config.yaml --output results/production --verbose
        """
    )
    
    parser.add_argument(
        '--config', '-c',
        type=str,
        required=True,
        help='Archivo de configuración para el pipeline'
    )
    
    parser.add_argument(
        '--output', '-o',
        type=str,
        required=True,
        help='Directorio de salida para los datos generados'
    )
    
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Activar logging detallado'
    )
    
    args = parser.parse_args()
    
    # Configurar nivel de logging
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    try:
        success = run_pipeline_and_validate(
            config_file=args.config,
            output_dir=args.output,
            verbose=args.verbose
        )
        
        if success:
            logger.info("🎉 Proceso completo exitoso!")
            return 0
        else:
            logger.error("💥 El proceso falló")
            return 1
            
    except Exception as e:
        logger.error(f"Error en proceso principal: {e}")
        return 1


if __name__ == "__main__":
    import pandas as pd  # Importar aquí para evitar error si no está disponible
    exit(main())
