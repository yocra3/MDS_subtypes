# Pipeline Unificado de Preparación de Datos para Modelos MDS

Este sistema unifica la preparación de datos para modelos de supervivencia en Síndrome Mielodisplásico (MDS), permitiendo generar datos tanto para modelos de redes neuronales tradicionales (CoxPH, Random Forest) como para modelos de redes neuronales de grafos (GNN).

## 🎯 Objetivos

- **Unificación**: Un solo pipeline para todos los tipos de modelos
- **Flexibilidad**: Configuración completa vía archivos YAML/JSON
- **Reproducibilidad**: Uso de folds de cross-validation pre-generados
- **Modularidad**: Componentes independientes y testeables
- **Mantenibilidad**: Código limpio y bien documentado

## 📁 Estructura del Sistema

```
scripts/3.GeneticScores/
├── prepare_data_unified.py     # Script principal del pipeline
├── configs/                    # Configuraciones
│   ├── default.yaml           # Configuración por defecto
│   ├── gnn_config.yaml        # Específica para modelos GNN
│   └── table_config.yaml      # Específica para modelos tabulares
├── utils/                      # Módulos utilitarios
│   ├── data_loader.py         # Carga de datos
│   ├── clinical_processor.py  # Procesamiento de variables clínicas
│   ├── data_generators.py     # Generación de datos (tablas/grafos)
│   └── create_cv_folds.R      # Generación de folds CV (ya implementado)
└── test/                       # Tests unitarios y de integración
    ├── test_data_loader.py    # Tests para carga de datos
    ├── test_pipeline_integration.py  # Tests de integración
    └── run_tests.py           # Script para ejecutar tests
```

## 🚀 Uso Básico

### Generar datos con configuración por defecto
```bash
python prepare_data_unified.py \
    --config configs/default.yaml \
    --output results/experiment1
```

### Generar solo grafos para GNN
```bash
python prepare_data_unified.py \
    --config configs/gnn_config.yaml \
    --output results/gnn_experiment \
    --formats graph
```

### Generar solo tablas para modelos tradicionales
```bash
python prepare_data_unified.py \
    --config configs/table_config.yaml \
    --output results/table_experiment \
    --formats table
```

### Ver configuración sin ejecutar
```bash
python prepare_data_unified.py \
    --config configs/default.yaml \
    --output /tmp/test \
    --dry-run
```

## ⚙️ Configuración

### Estructura de Configuración

Las configuraciones se organizan en secciones:

- **`data_paths`**: Rutas a archivos de datos
- **`variable_processing`**: Procesamiento de variables clínicas y selección de genes
- **`data_generation`**: Configuración de formatos de salida (tablas/grafos)
- **`cross_validation`**: Configuración de CV (usa folds del archivo clínico)
- **`output`**: Configuración de archivos de salida
- **`logging`**: Configuración de logging

### Configuraciones Disponibles

1. **`default.yaml`**: Configuración balanceada que genera tanto tablas como grafos
2. **`table_config.yaml`**: Optimizada para modelos NN tradicionales con selección manual de genes
3. **`gnn_config.yaml`**: Optimizada para modelos GNN con embeddings y configuración permisiva

### Diferencias Principales entre Configuraciones

| Aspecto | `table_config.yaml` | `gnn_config.yaml` | `default.yaml` |
|---------|-------------------|------------------|----------------|
| **Input files** | Solo archivo clínico | Incluye VAF y embeddings | Incluye VAF y embeddings |
| **Output** | Solo tablas | Solo grafos | Tablas + grafos |
| **Selección genes** | Manual (10 genes específicos) | Automática (hasta 100) | Automática (hasta 50) |
| **VAF usage** | No se usa | Filtro (≥3%) | Filtro (≥5%) |
| **Gene embeddings** | No se usan | Disponibles | Disponibles |
| **Criterios inclusión** | Restrictivos (10% freq) | Permisivos (2% freq) | Balanceados (5% freq) |

## 📊 Formatos de Salida

### Tablas (para modelos tradicionales)
- **Archivos**: CSV con matriz de features por paciente
- **Features**: Variables clínicas procesadas + genes seleccionados (manual o automático)
- **Selección de genes**: Configurada en `variable_processing.gene_selection`
- **Cross-validation**: Folds incluidos en el archivo de salida
- **Supervivencia**: Variables OS/OS_STATUS como targets (no predictores)
- **Datos**: Completos (sin imputación)

### Grafos (para modelos GNN)
- **Archivos**: Pickle con objetos PyTorch Geometric
- **Nodos**: 
  - Pacientes: features clínicas ± mutaciones binarizadas
  - Genes: estadísticas de VAF y/o embeddings pre-entrenados
- **Edges**: 
  - Paciente-gen: edges binarios (1 si hay mutación filtrada por VAF, 0 si no)
  - **NO multi-edge**: Cada paciente tiene un grafo único y simple
- **Cross-validation**: Folds incluidos en el archivo principal
- **Supervivencia**: Variables OS/OS_STATUS como targets (no predictores)
- **Gene embeddings**: Opcional, desde archivo scGPT

### Personalización

Puedes crear tu propia configuración copiando una existente y modificándola:

```bash
cp configs/default.yaml configs/mi_experimento.yaml
# Editar mi_experimento.yaml según necesidades
python prepare_data_unified.py --config configs/mi_experimento.yaml --output results/mi_experimento
```

## 🧪 Tests

### Ejecutar todos los tests
```bash
cd test/
python run_tests.py
```

### Ejecutar tests específicos
```bash
python run_tests.py --test test_data_loader
```

### Ver tests disponibles
```bash
python run_tests.py --list
```

## 📊 Archivos de Input

El sistema trabaja con **dos archivos principales** que ya contienen el preprocesado necesario:

### Archivo 1: Variables Clínicas, Mutaciones Binarizadas y Folds CV
- **Ruta**: `data/processed/clinical_mutations_combined.tsv`
- **Contenido**: Variables clínicas + kariotipos + mutaciones binarizadas + folds de cross-validation
- **Uso**: 
  - Input para modelos NN tradicionales (output tabla)
  - Variables que definen nodos de pacientes en GNN
  - Información de folds para cross-validation

### Archivo 2: Mutaciones Priorizadas con VAF
- **Ruta**: `data/processed/mutations_prioritized.tsv`  
- **Contenido**: Mutaciones priorizadas por paciente con información de VAF
- **Uso**: 
  - **Para GNN**: VAF se usa como **filtro** para seleccionar qué mutaciones incluir
  - **Para tablas**: No se usa (las mutaciones ya están binarizadas en el archivo 1)

## 🔍 Uso de VAF (Variant Allele Frequency)

**Importante**: La VAF se usa únicamente como **criterio de filtrado**, NO como peso de edges:

- ✅ **Filtro de calidad**: Solo mutaciones con VAF ≥ umbral se incluyen en el grafo
- ✅ **Nodos de genes**: Estadísticas de VAF pueden incluirse como features de nodos
- ❌ **NO como peso de edges**: Los edges son binarios (0/1 si hay mutación)

```yaml
edge_criteria:
  use_vaf_filter: true          # Usar VAF para filtrar mutaciones
  min_vaf_threshold: 0.05       # Solo mutaciones con VAF ≥ 5%
  # Los edges resultantes son binarios, no pesados por VAF
```

### Ejemplo Comparativo:

**❌ Uso INCORRECTO (edges pesados por VAF):**
```
Paciente1 ----[weight=0.23]---- GenTP53
Paciente2 ----[weight=0.45]---- GenTP53
```

**✅ Uso CORRECTO (VAF como filtro, edges binarios):**
```
# Solo se incluyen mutaciones con VAF ≥ umbral
Paciente1 ----[1]---- GenTP53  # VAF=0.23 ≥ 0.05 → incluida
Paciente2 ----[1]---- GenTP53  # VAF=0.45 ≥ 0.05 → incluida
Paciente3      X      GenTP53  # VAF=0.02 < 0.05 → excluida
```

## 🧬 Selección de Genes

La selección de genes se configura en la sección `variable_processing.gene_selection`:

### Para Modelos Tabulares (Manual)
```yaml
variable_processing:
  gene_selection:
    # Lista específica de genes para modelos tabulares
    selected_genes: 
      - "TP53"
      - "ASXL1" 
      - "SF3B1"
      - "DNMT3A"
      - "SRSF2"
```

### Para Modelos GNN (Automática)
```yaml
variable_processing:
  gene_selection:
    # null = incluir todos los genes disponibles
    selected_genes: null
    
    # Criterios para selección automática
    auto_selection_criteria:
      min_mutation_frequency: 0.02  # 2% mínimo
      max_genes: 100  # Máximo 100 genes
```

## 🔍 Gene Embeddings (para GNN)

Los modelos GNN pueden usar embeddings pre-entrenados como features iniciales para los nodos de genes:

```yaml
data_paths:
  # Archivo con embeddings de genes (opcional)
  gene_embeddings_file: "data/scGPT/gene_embeddings.npy"

data_generation:
  graph_format:
    gene_node_features:
      feature_type: "both"  # vaf_stats + embeddings
      embedding_source: "scgpt"
      embedding_dimension: 512
```

### Tipos de Features para Nodos de Genes:
- **`vaf_stats`**: Solo estadísticas de VAF (media, std, freq)
- **`embeddings`**: Solo embeddings pre-entrenados
- **`both`**: Combinación de estadísticas + embeddings

```yaml
table_format:
  selected_genes: [
    "TP53",      # Gen supresor tumoral
    "ASXL1",     # Regulador epigenético  
    "SF3B1",     # Factor de splicing
    "DNMT3A",    # Metiltransferasa
    "TET2",      # Demetilasa
    "SRSF2",     # Factor de splicing
    "RUNX1",     # Factor transcripcional
    "BCOR",      # Corepressor
    "STAG2"      # Cohesinas
  ]
```

- **`null`**: Incluye todos los genes disponibles en el archivo
- **Lista específica**: Solo incluye los genes nombrados
- **Validación**: El sistema verifica que los genes existen como columnas en el archivo

## 🎯 Variables de Supervivencia

**Importante**: Las variables de supervivencia son **outputs del modelo**, NO predictores:

```yaml
survival_variables:
  time_variable: "OS"          # Tiempo hasta muerte/último seguimiento
  status_variable: "OS_STATUS" # 0=censored, 1=death
```

### Variables NO Usadas como Predictores:
- ❌ **Scores de riesgo**: IPSS, IPSS_R, IPSS_M (ya incorporan pronóstico)
- ❌ **Variables de supervivencia**: OS, OS_STATUS (son targets)
- ❌ **Otras supervivencias**: EFS, EFS_STATUS

## ⚙️ Configuración Simplificada

### Estructura Basada en Propiedades Estadísticas

La configuración se organiza por **tipo de variable** según sus propiedades estadísticas:

```yaml
variable_processing:
  # Variables continuas - requieren escalado
  continuous_variables:
    - "Age"
    - "HB"        # Hemoglobina
    - "ANC"       # Neutrófilos absolutos
    - "PLT"       # Plaquetas
    - "BM_BLAST"  # Blastos en médula ósea
    - "PB_BLAST"  # Blastos en sangre periférica
    # NOTA: Scores de riesgo (IPSS, IPSS_R, etc.) NO se usan como predictores
  
  # Variables categóricas - requieren encoding
  categorical_variables:
    - "Sex"
    - "Cytogenetics_IPSS"
    - "WHO_Classification"
  
  # Variables binarias - no requieren procesamiento especial
  binary_variables:
    - "Complex_Karyotype"
    # Las mutaciones binarizadas se detectan automáticamente
  
  # Variables de supervivencia (son OUTPUTS del modelo, NO predictores)
  survival_variables:
    time_variable: "OS"          # Tiempo hasta el evento
    status_variable: "OS_STATUS" # Estado del evento (0=censored, 1=event)
  
  # Opciones de procesamiento
  processing_options:
    scaling_method: "standard"     # Para continuas
    encoding_method: "onehot"     # Para categóricas  
    # NOTA: No hay imputación - se trabaja con datos completos
```

## 🔄 Estado del Desarrollo

### ✅ Completado (Paso 2A)
- [x] Estructura base del sistema
- [x] Archivo principal `prepare_data_unified.py`
- [x] Configuraciones YAML para diferentes casos de uso
- [x] Módulos placeholder en `utils/`
- [x] Framework de tests básico
- [x] Documentación inicial

### 🚧 Por Implementar
- [ ] **Paso 2B**: Migrar carga de datos de R a Python
- [ ] **Paso 2C**: Implementar sistema de configuración flexible
- [ ] **Paso 2D**: Migrar procesamiento de variables clínicas
- [ ] **Paso 2F**: Implementar generadores de datos (tablas/grafos)
- [ ] **Paso 2G**: Integrar todo en script principal
- [ ] **Paso 2H**: Validación y comparación end-to-end

## 🔗 Dependencias

### Python
- pandas
- numpy
- scikit-learn
- torch
- torch-geometric
- networkx
- pyyaml

### R (para folds pre-generados)
- Los folds ya están generados en el paso 1

## 📝 Notas Importantes

1. **Folds de CV**: El sistema usa folds de cross-validation pre-generados en R. No los regenera, solo los carga y usa.

2. **Compatibilidad**: Diseñado para ser compatible con el código R existente mientras migra gradualmente a Python.

3. **Flexibilidad**: Cada componente puede configurarse independientemente para máxima flexibilidad experimental.

4. **Validación**: Incluye validación automática de datos generados para asegurar calidad.

## 🆘 Soporte

Para problemas o preguntas sobre el pipeline unificado, consulta:
- Los logs generados en `prepare_data_unified.log`
- Los tests en la carpeta `test/`
- Las configuraciones de ejemplo en `configs/`
