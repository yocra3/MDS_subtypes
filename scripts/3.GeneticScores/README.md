# Pipeline Unificado de Preparación de Datos MDS

Este directorio contiene el pipeline unificado para la preparación de datos de síndrome mielodisplásico (MDS) que permite generar datos compatibles con modelos de supervivencia tanto de redes neuronales (NN) como de redes de grafos neuronales (GNN).

## 📋 Tabla de Contenidos

- [Descripción General](#descripción-general)
- [Características](#características)
- [Estructura del Proyecto](#estructura-del-proyecto)
- [Instalación y Configuración](#instalación-y-configuración)
- [Guía de Uso](#guía-de-uso)
- [Ejemplos Prácticos](#ejemplos-prácticos)
- [Validación y Pruebas](#validación-y-pruebas)
- [Troubleshooting](#troubleshooting)

## 🎯 Descripción General

El pipeline unificado permite:

1. **Preparación flexible de datos**: Selección configurable de variables clínicas y genéticas
2. **Múltiples formatos de salida**: Generación de datos para modelos de tabla (NN) y grafo (GNN)
3. **Cross-validation integrada**: Uso de folds pre-generados para validación cruzada
4. **Validación automatizada**: Verificación de calidad de datos y compatibilidad con modelos
5. **Configuración YAML**: Control total mediante archivos de configuración

## ✨ Características

- ✅ **Modular y extensible**: Arquitectura basada en componentes reutilizables
- ✅ **Docker compatible**: Ejecución en contenedores para reproducibilidad
- ✅ **Logging detallado**: Seguimiento completo del procesamiento
- ✅ **Validación exhaustiva**: Verificación automática de calidad y compatibilidad
- ✅ **Flexible**: Soporte para diferentes tipos de variables y embeddings

## 📁 Estructura del Proyecto

```
scripts/3.GeneticScores/
├── prepare_data_unified.py          # Script principal del pipeline
├── utils/                           # Módulos utilitarios
│   ├── data_loader.py              # Carga de datos (clínicos, mutaciones, embeddings)
│   ├── clinical_processor.py       # Procesamiento de variables clínicas
│   └── data_generators.py          # Generación de formatos tabla/grafo
├── configs/                         # Archivos de configuración
│   ├── real_data_config.yaml       # Configuración para datos reales
│   ├── test_config.yaml           # Configuración para pruebas
│   └── minimal_config.yaml        # Configuración mínima
├── tests/                          # Scripts de validación
│   ├── validate_step2b.py          # Validación de integración
│   ├── validate_step3c.py          # Validación de calidad de datos
│   ├── validate_step3d.py          # Validación de compatibilidad con modelos
│   └── run_pipeline_and_validate.py # Ejecución automatizada
├── examples/                       # Ejemplos de uso
│   ├── basic_usage.py              # Ejemplo básico
│   ├── advanced_usage.py           # Ejemplo avanzado
│   └── custom_config_example.py    # Configuración personalizada
└── docs/                          # Documentación adicional
    ├── configuration_guide.md      # Guía de configuración
    ├── validation_guide.md         # Guía de validación
    └── docker_usage.md            # Uso con Docker
```

## 🚀 Instalación y Configuración

### Requisitos del Sistema

- Python 3.8+
- Docker (recomendado para reproducibilidad)
- Dependencias de Python (ver `requirements.txt`)

### Instalación Local

```bash
# Clonar repositorio (si no lo tienes)
git clone <repo-url>
cd MDS_subtypes/scripts/3.GeneticScores

# Instalar dependencias
pip install -r requirements.txt

# Verificar instalación
python prepare_data_unified.py --help
```

### Uso con Docker (Recomendado)

```bash
# Construir imagen (si no existe)
docker build -t mds_subtypes_python:1.6 .

# Verificar imagen
docker run --rm mds_subtypes_python:1.6 python --version
```

## 📖 Guía de Uso

### 1. Configuración Básica

Crea un archivo de configuración YAML:

```yaml
# configs/mi_config.yaml
data_sources:
  clinical_data: "data/IPSSMol/df_clinical.tsv"
  mutation_vaf_data: "data/IPSSMol/maf.tsv"
  gene_embeddings: "data/scGPT/gene_embeddings.csv"

variable_processing:
  gene_selection:
    - "TP53"
    - "ASXL1"
    - "RUNX1"
  
  continuous_variables:
    - "AGE"
    - "HB"
    - "PLT"
  
  categorical_variables:
    - "SEX"
  
  binary_variables:
    - "AML_TRANSF"
  
  survival_variables:
    time_variable: "LFS_YEARS"
    status_variable: "LFS_STATUS"

cross_validation:
  fold_column: "fold"

data_generation:
  output_formats: ["table", "graph"]
  table_format:
    output_filename: "mi_experimento_table"
  graph_format:
    output_filename: "mi_experimento_graph"

output:
  file_prefix: "mi_experimento"
```

### 2. Ejecución Básica

```bash
# Ejecución local
python prepare_data_unified.py \
    --config configs/mi_config.yaml \
    --output results/mi_experimento

# Ejecución con Docker
docker run --rm \
    -v /ruta/completa/MDS_subtypes:/workspace \
    mds_subtypes_python:1.6 \
    python /workspace/scripts/3.GeneticScores/prepare_data_unified.py \
    --config /workspace/scripts/3.GeneticScores/configs/mi_config.yaml \
    --output /workspace/results/mi_experimento
```

### 3. Validación de Resultados

```bash
# Validación de calidad
python tests/validate_step3c.py \
    --data-dir results/mi_experimento

# Validación de compatibilidad con modelos
python tests/validate_step3d.py \
    --data-dir results/mi_experimento

# Validación completa automatizada
python tests/run_pipeline_and_validate.py \
    --config configs/mi_config.yaml \
    --output results/mi_experimento_validado
```

## 💡 Ejemplos Prácticos

### Ejemplo 1: Análisis Básico con Genes Específicos

```python
# examples/basic_analysis.py
import yaml
from pathlib import Path
import subprocess

# Configuración básica
config = {
    'data_sources': {
        'clinical_data': 'data/IPSSMol/df_clinical.tsv',
        'mutation_vaf_data': 'data/IPSSMol/maf.tsv',
        'gene_embeddings': 'data/scGPT/gene_embeddings.csv'
    },
    'variable_processing': {
        'gene_selection': ['TP53', 'ASXL1', 'RUNX1', 'SRSF2', 'SF3B1'],
        'continuous_variables': ['AGE', 'HB', 'PLT', 'ANC'],
        'categorical_variables': ['SEX'],
        'binary_variables': ['AML_TRANSF'],
        'survival_variables': {
            'time_variable': 'LFS_YEARS',
            'status_variable': 'LFS_STATUS'
        }
    },
    'cross_validation': {'fold_column': 'fold'},
    'data_generation': {
        'output_formats': ['table', 'graph'],
        'table_format': {'output_filename': 'basic_table'},
        'graph_format': {'output_filename': 'basic_graph'}
    },
    'output': {'file_prefix': 'basic_analysis'}
}

# Guardar configuración
config_file = Path('configs/basic_analysis.yaml')
with open(config_file, 'w') as f:
    yaml.dump(config, f)

# Ejecutar pipeline
result = subprocess.run([
    'python', 'prepare_data_unified.py',
    '--config', str(config_file),
    '--output', 'results/basic_analysis',
    '--verbose'
], capture_output=True, text=True)

print("STDOUT:", result.stdout)
if result.stderr:
    print("STDERR:", result.stderr)
```

### Ejemplo 2: Análisis Completo con Validación

```python
# examples/complete_analysis.py
import subprocess
from pathlib import Path

def run_complete_analysis():
    """Ejecuta análisis completo con validación."""
    
    # 1. Ejecutar pipeline
    print("🚀 Ejecutando pipeline de datos...")
    pipeline_result = subprocess.run([
        'python', 'prepare_data_unified.py',
        '--config', 'configs/real_data_config.yaml',
        '--output', 'results/complete_analysis'
    ], capture_output=True, text=True)
    
    if pipeline_result.returncode != 0:
        print("❌ Error en pipeline:", pipeline_result.stderr)
        return False
    
    # 2. Validar calidad
    print("🔍 Validando calidad de datos...")
    quality_result = subprocess.run([
        'python', 'tests/validate_step3c.py',
        '--data-dir', 'results/complete_analysis'
    ], capture_output=True, text=True)
    
    # 3. Validar compatibilidad
    print("🧪 Validando compatibilidad con modelos...")
    compat_result = subprocess.run([
        'python', 'tests/validate_step3d.py',
        '--data-dir', 'results/complete_analysis'
    ], capture_output=True, text=True)
    
    # 4. Resumen
    print("\n" + "="*50)
    print("RESUMEN DEL ANÁLISIS COMPLETO")
    print("="*50)
    print(f"Pipeline: {'✅ OK' if pipeline_result.returncode == 0 else '❌ ERROR'}")
    print(f"Calidad: {'✅ OK' if quality_result.returncode == 0 else '❌ ERROR'}")
    print(f"Compatibilidad: {'✅ OK' if compat_result.returncode == 0 else '❌ ERROR'}")
    
    return all(r.returncode == 0 for r in [pipeline_result, quality_result, compat_result])

if __name__ == "__main__":
    success = run_complete_analysis()
    exit(0 if success else 1)
```

### Ejemplo 3: Análisis con Docker

```bash
#!/bin/bash
# examples/docker_analysis.sh

# Variables
IMAGE="mds_subtypes_python:1.6"
WORKSPACE="/workspace"
CONFIG="configs/real_data_config.yaml"
OUTPUT="results/docker_analysis"

echo "🐳 Ejecutando análisis completo con Docker..."

# 1. Pipeline principal
echo "📊 Ejecutando pipeline..."
docker run --rm \
    -v "$(pwd)/../..:/workspace" \
    $IMAGE \
    python $WORKSPACE/scripts/3.GeneticScores/prepare_data_unified.py \
    --config $WORKSPACE/scripts/3.GeneticScores/$CONFIG \
    --output $WORKSPACE/$OUTPUT \
    --verbose

# 2. Validación de calidad
echo "🔍 Validando calidad..."
docker run --rm \
    -v "$(pwd)/../..:/workspace" \
    $IMAGE \
    python $WORKSPACE/scripts/3.GeneticScores/tests/validate_step3c.py \
    --data-dir $WORKSPACE/$OUTPUT

# 3. Validación de compatibilidad
echo "🧪 Validando compatibilidad..."
docker run --rm \
    -v "$(pwd)/../..:/workspace" \
    $IMAGE \
    python $WORKSPACE/scripts/3.GeneticScores/tests/validate_step3d.py \
    --data-dir $WORKSPACE/$OUTPUT

echo "✅ Análisis completado. Revisa $OUTPUT para los resultados."
```

## 🔬 Validación y Pruebas

### Tipos de Validación

1. **Validación de Integración** (`validate_step2b.py`)
   - Verifica que todos los módulos funcionan correctamente
   - Prueba el pipeline completo con datos de prueba

2. **Validación de Calidad** (`validate_step3c.py`)
   - Verifica dimensiones, tipos de datos y estructura
   - Valida datos de supervivencia y folds
   - Detecta valores faltantes o inconsistencias

3. **Validación de Compatibilidad** (`validate_step3d.py`)
   - Prueba compatibilidad con modelos NN y GNN existentes
   - Verifica forward pass de modelos
   - Valida estructuras de datos específicas

### Comandos de Validación

```bash
# Validación rápida (solo integración)
python tests/validate_step2b.py

# Validación de calidad de datos específicos
python tests/validate_step3c.py \
    --data-dir results/mi_experimento \
    --table-file mi_tabla.pkl \
    --graph-file mi_grafo.pt

# Validación completa automatizada
python tests/run_pipeline_and_validate.py \
    --config configs/mi_config.yaml \
    --output results/validacion_completa \
    --run-quality-validation \
    --run-model-validation
```

## 🔧 Troubleshooting

### Problemas Comunes

#### Error: "Archivo no encontrado"
```bash
# Verificar rutas en configuración
python prepare_data_unified.py --config configs/mi_config.yaml --dry-run

# Verificar estructura de archivos
ls -la data/IPSSMol/
```

#### Error: "Columnas faltantes"
```python
# Verificar columnas disponibles
import pandas as pd
df = pd.read_csv('data/IPSSMol/df_clinical.tsv', sep='\t')
print("Columnas disponibles:", df.columns.tolist())
```

#### Error de memoria con grafos grandes
```yaml
# Reducir tamaño en configuración
data_generation:
  graph_format:
    patient_feature_type: "clinical_only"  # Usar menos features
```

#### Incompatibilidad con modelos
```bash
# Verificar dimensiones esperadas
python tests/validate_step3d.py --data-dir results/mi_experimento --verbose

# Ajustar configuración según reporte de validación
```

### Logs y Depuración

```bash
# Activar logging detallado
python prepare_data_unified.py --config configs/mi_config.yaml --output results/debug --verbose

# Revisar logs
tail -f prepare_data_unified.log

# Verificar archivos generados
ls -la results/debug/
cat results/debug/quality_validation_report.txt
```

### Soporte Docker

```bash
# Verificar imagen Docker
docker images | grep mds_subtypes

# Probar container interactivo
docker run -it --rm -v $(pwd):/workspace mds_subtypes_python:1.6 bash

# Dentro del container:
cd /workspace/scripts/3.GeneticScores
python prepare_data_unified.py --help
```

## 📞 Contacto y Contribuciones

Para reportar bugs, solicitar features o contribuir al proyecto:

1. Crear issue en el repositorio
2. Documentar el problema con ejemplos reproducibles
3. Incluir logs relevantes y configuración utilizada

---

**Nota**: Esta documentación cubre la versión actual del pipeline. Para actualizaciones, consulta el changelog y la documentación técnica en `/docs/`.
