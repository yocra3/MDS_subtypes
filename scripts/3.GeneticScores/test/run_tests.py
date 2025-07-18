#!/usr/bin/env python3
"""
Script para ejecutar todos los tests del pipeline unificado.

Autor: Sistema automatizado
Fecha: 2025-07-04
"""

import unittest
import sys
from pathlib import Path
import argparse

# Añadir directorios al path
current_dir = Path(__file__).parent
sys.path.append(str(current_dir.parent))
sys.path.append(str(current_dir.parent / "utils"))


def run_all_tests(verbose=False):
    """
    Ejecutar todos los tests disponibles.
    
    Args:
        verbose: Si mostrar salida detallada
    """
    print("=== Ejecutando Tests del Pipeline Unificado MDS ===")
    
    # Descubrir y ejecutar tests
    loader = unittest.TestLoader()
    start_dir = str(current_dir)
    
    # Cargar todos los tests
    suite = loader.discover(start_dir, pattern='test_*.py')
    
    # Ejecutar tests
    runner = unittest.TextTestRunner(
        verbosity=2 if verbose else 1,
        stream=sys.stdout,
        descriptions=True,
        failfast=False
    )
    
    result = runner.run(suite)
    
    # Mostrar resumen
    print(f"\n=== Resumen de Tests ===")
    print(f"Tests ejecutados: {result.testsRun}")
    print(f"Errores: {len(result.errors)}")
    print(f"Fallos: {len(result.failures)}")
    print(f"Omitidos: {len(result.skipped) if hasattr(result, 'skipped') else 0}")
    
    if result.errors:
        print(f"\nErrores encontrados:")
        for test, error in result.errors:
            error_first_line = error.split('\n')[0]
            print(f"  - {test}: {error_first_line}")
    
    if result.failures:
        print(f"\nFallos encontrados:")
        for test, failure in result.failures:
            failure_first_line = failure.split('\n')[0]
            print(f"  - {test}: {failure_first_line}")
    
    # Determinar éxito
    success = len(result.errors) == 0 and len(result.failures) == 0
    print(f"\nResultado: {'✅ ÉXITO' if success else '❌ FALLO'}")
    
    return success


def run_specific_test(test_name, verbose=False):
    """
    Ejecutar un test específico.
    
    Args:
        test_name: Nombre del archivo de test (sin .py)
        verbose: Si mostrar salida detallada
    """
    print(f"=== Ejecutando Test Específico: {test_name} ===")
    
    try:
        # Cargar módulo de test específico
        loader = unittest.TestLoader()
        suite = loader.loadTestsFromName(test_name)
        
        # Ejecutar test
        runner = unittest.TextTestRunner(
            verbosity=2 if verbose else 1,
            stream=sys.stdout,
            descriptions=True
        )
        
        result = runner.run(suite)
        
        success = len(result.errors) == 0 and len(result.failures) == 0
        print(f"\nResultado: {'✅ ÉXITO' if success else '❌ FALLO'}")
        
        return success
        
    except Exception as e:
        print(f"Error ejecutando test {test_name}: {e}")
        return False


def main():
    """Función principal."""
    parser = argparse.ArgumentParser(
        description="Ejecutar tests del pipeline unificado MDS",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Ejemplos de uso:
  python run_tests.py                           # Ejecutar todos los tests
  python run_tests.py --verbose                 # Ejecutar con salida detallada
  python run_tests.py --test test_data_loader   # Ejecutar test específico
        """
    )
    
    parser.add_argument(
        '--test', '-t',
        type=str,
        help='Ejecutar test específico (nombre del archivo sin .py)'
    )
    
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Mostrar salida detallada de los tests'
    )
    
    parser.add_argument(
        '--list', '-l',
        action='store_true',
        help='Listar tests disponibles'
    )
    
    args = parser.parse_args()
    
    if args.list:
        print("Tests disponibles:")
        test_files = list(current_dir.glob('test_*.py'))
        for test_file in test_files:
            print(f"  - {test_file.stem}")
        return 0
    
    try:
        if args.test:
            success = run_specific_test(args.test, args.verbose)
        else:
            success = run_all_tests(args.verbose)
        
        return 0 if success else 1
        
    except KeyboardInterrupt:
        print("\n\nTests interrumpidos por el usuario.")
        return 1
    except Exception as e:
        print(f"\nError ejecutando tests: {e}")
        return 1


if __name__ == "__main__":
    exit(main())
