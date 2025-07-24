Plan de Implementación: Calculadora de Riesgo MDS con Bokeh
Resumen General
Se propone desarrollar una aplicación web interactiva usando Bokeh para calcular el riesgo de pacientes con Síndrome Mielodisplásico (MDS), replicando la funcionalidad de https://mds-risk-model.com/. La aplicación ejecutará tres modelos: un modelo PyTorch Lightning GNN pre-entrenado y dos modelos adicionales que se implementarán posteriormente como scripts de Python.

Requisitos
Requisitos Funcionales
Interfaz de entrada de datos del paciente: Campos para datos clínicos, citogenéticos y mutacionales
Cálculo de riesgo en tiempo real: Procesamiento de datos y predicción con tres modelos
Visualización de resultados: Gráficos interactivos mostrando estratificación de riesgo, curvas de supervivencia y hazard ratios
Comparación de modelos: Mostrar resultados de los tres modelos de forma comparativa
Exportación de resultados: Posibilidad de descargar reportes en PDF
Requisitos Técnicos
Framework: Bokeh Server para la aplicación web
Modelos: Integración con PyTorch Lightning para el modelo GNN
Datos: Procesamiento en tiempo real de variables clínicas y mutacionales
Performance: Respuesta rápida (<5 segundos) para el cálculo de riesgo
Pasos de Implementación
Paso 1: Estructura del Proyecto y Configuración Base
Archivos a crear/modificar:

app/ (directorio principal de la aplicación)
main.py (aplicación principal de Bokeh)
config.yaml (configuración de modelos y variables)
models/ (directorio para modelos)
gnn_model.py (wrapper para modelo PyTorch Lightning)
legacy_models.py (placeholder para otros dos modelos)
utils/ (utilidades)
data_processor.py (procesamiento de datos de entrada)
visualization.py (funciones de visualización)
static/ (archivos estáticos CSS/JS si es necesario)
Configuración:

Configurar entorno Bokeh Server
Definir estructura de datos de entrada basada en el código existente
Configurar carga de modelos pre-entrenados

Paso 2: Implementación del Wrapper para Modelo GNN
Archivos principales:

app/models/gnn_model.py
Tareas:

Crear clase wrapper que cargue el modelo PyTorch Lightning desde checkpoint
Implementar método de predicción que tome datos del paciente y retorne score de riesgo
Integrar con el código existente de evalGNN_classes.py y las clases de evaluación
Manejar el preprocesamiento de datos (transformación a grafos heterogéneos)
Implementar categorización de riesgo basada en scores
Detalles de implementación:

Utilizar BaseEvalGNN o crear wrapper similar
Cargar checkpoint más reciente de model_checkpoints
Procesar variables clínicas usando ClinicalProcessor
Transformar datos de mutaciones a formato de grafo



Paso 3: Interfaz de Usuario con Bokeh
Archivos principales:

app/main.py
app/utils/visualization.py
Componentes de la interfaz:

Panel de entrada de datos:

Campos numéricos: Edad, % blastos en médula ósea, WBC, ANC, etc.
Selectores categóricos: Sexo, cariotipo, riesgo citogenético
Checklist de mutaciones: TP53, SF3B1, RUNX1, etc. (basado en genes IPSS-M)
VAF (Variant Allele Frequency) para mutaciones presentes
Panel de resultados:

Score de riesgo numérico
Categoría de riesgo (Very Low, Low, Moderate Low, etc.)
Hazard ratio comparado con paciente promedio
Gráfico de distribución de riesgo con posición del paciente
Panel de supervivencia:

Tabla con medianas de supervivencia por categoría
Gráficos de curvas de supervivencia (Overall Survival, Leukemia-Free Survival)
Probabilidades de transformación a AML por año


Paso 4: Procesamiento de Datos y Lógica de Negocio
Archivos principales:

app/utils/data_processor.py
Funcionalidades:

Validación de datos de entrada:

Rangos válidos para variables numéricas
Combinaciones válidas de mutaciones
Manejo de valores faltantes
Transformación de datos:

Normalización de variables continuas
Codificación de variables categóricas
Creación de estructura de grafo para modelo GNN
Cálculo de scores:

Interfaz unificada para los tres modelos
Manejo de errores y valores fuera de rango
Logging de predicciones

Paso 5: Visualizaciones Interactivas
Archivos principales:

app/utils/visualization.py
Gráficos a implementar:

Distribución de riesgo: Histograma con marcador de posición del paciente
Curvas de supervivencia: Gráficos interactivos con Bokeh
Tabla de estratificación: Tabla dinámica con highlighting de categoría del paciente
Comparación de modelos: Panel que muestre resultados de los tres modelos lado a lado

Paso 6: Configuración de Despliegue
Archivos principales:

Dockerfile
app/deploy.py
Tareas:

Crear Dockerfile específico para la aplicación Bokeh
Configurar Bokeh Server para producción
Implementar logging y monitoreo
Configurar variables de entorno para paths de modelos

Paso 7: Placeholder para Modelos Adicionales
Archivos principales:

app/models/legacy_models.py
Implementación inicial:

Crear interfaces mock que retornen valores placeholder
Definir estructura esperada para cuando se implementen los modelos reales
Documentar requirements para los scripts de Python futuros

Paso 8: Testing y Validación
Archivos principales:

app/tests/ (directorio de tests)
Tests a implementar:

Test de carga de modelos
Test de procesamiento de datos
Test de interfaz de usuario (casos límite)
Test de comparación con resultados conocidos
Consideraciones Técnicas
Dependencias del Modelo GNN
Basado en el análisis del código, el modelo GNN requiere:
Variables clínicas: SEX, BM_BLAST, WBC, ANC, karyotype_risk, etc.
Datos de mutaciones: boolean o VAF para genes específicos
Embeddings de genes (pre-computados con scGPT)
Performance
Los modelos PyTorch Lightning pueden requerir GPU para mejor performance
Implementar cache para embeddings de genes pre-computados
Considerar batch processing si se esperan múltiples usuarios simultáneos
Escalabilidad
Estructura modular permite fácil adición de nuevos modelos
Configuración mediante YAML permite ajustes sin cambios de código
Separación clara entre lógica de modelos y presentación
Archivos de Configuración
app/config.yaml
Este plan proporciona una base sólida para implementar la calculadora de riesgo MDS, manteniendo compatibilidad con el código existente y permitiendo futuras extensiones con los otros dos modelos.




Propuesta de Implementación - Paso 1: Estructura del Proyecto y Configuración Base
Sub-paso 1.1: Creación de la Estructura de Directorios
Objetivo: Establecer la arquitectura base del proyecto

Directorios a crear:

Justificación: Separación clara de responsabilidades siguiendo principios de arquitectura modular.

Sub-paso 1.2: Configuración Principal (config.yaml)
Objetivo: Definir configuración centralizada basada en el análisis del código existente

Estructura propuesta:

Sección 1: Modelos
Modelo GNN: Configuración para cargar checkpoints de PyTorch Lightning
Rutas a checkpoints (boolean y VAF)
Hiperparámetros (hidden_dim, learning_rate, etc.)
Configuración de GPU/CPU
Modelos Legacy: Placeholders para futuros modelos
Configuración mock inicial
Estructura esperada para implementación futura
Sección 2: Variables Clínicas
Basado en el análisis del código:

Variables continuas: AGE, BM_BLAST, WBC, ANC, HB, PLT, MONOCYTES
Variables categóricas: SEX, CYTO_IPSSR (Very Good, Good, Intermediate, Poor, Very Poor)
Variables binarias: plus8, del7, del20q, del7q, delY, del5q, complex
Rangos de validación: Min/max para cada variable numérica
Sección 3: Genes de Mutación
Lista de genes IPSS-M identificados en el código:

Genes principales: TP53, SF3B1, RUNX1, ASXL1, SRSF2, U2AF1, etc.
Configuración VAF: rangos válidos (0-100%)
Categorización por frecuencia y relevancia clínica
Sección 4: Categorización de Riesgo
Basado en IPSS-M estándar:

Very Low (≤ -1.5), Low (-1.5 to -0.5), Moderate Low (-0.5 to 0)
Moderate High (0 to 0.5), High (0.5 to 1.5), Very High (> 1.5)
Colores y etiquetas para visualización
Sección 5: Configuración de Bokeh Server
Puerto, host, configuración de desarrollo/producción
Configuración de logging y debugging

Sub-paso 1.3: Aplicación Principal (main.py)
Objetivo: Punto de entrada de la aplicación Bokeh

Componentes principales:

1. Inicialización
Carga de configuración desde config.yaml
Inicialización de modelos (con manejo de errores)
Setup de logging
2. Layout Principal
Panel izquierdo: Entrada de datos del paciente
Pestañas: Datos Clínicos, Mutaciones, Cariotipo
Panel derecho: Resultados y visualizaciones
Pestañas: Scores, Supervivencia, Comparación de Modelos
3. Callbacks
Función principal de cálculo de riesgo
Actualización de visualizaciones en tiempo real
Validación de entrada de datos
4. Manejo de Estado
Session state para datos del paciente
Cache para resultados de modelos
Manejo de errores y mensajes de usuario

Sub-paso 1.4: Wrapper del Modelo GNN (models/gnn_model.py)
Objetivo: Interfaz unificada para el modelo PyTorch Lightning

Componentes:

Clase MDSGNNModel
Inicialización: Carga de checkpoint y configuración
Método predict(): Toma datos del paciente, retorna score de riesgo
Método get_risk_category(): Convierte score a categoría de riesgo
Método validate_input(): Valida datos de entrada
Integración con Código Existente
Utilizar clases de evalGNN_classes.py
Integrar con ClinicalProcessor para preprocesamiento
Manejar transformación a grafos heterogéneos
Manejo de Datos
Transformación de datos clínicos a formato esperado por GNN
Creación de estructura de grafo paciente-genes
Normalización y escalado de variables


Sub-paso 1.5: Modelos Legacy (models/legacy_models.py)
Objetivo: Placeholders para futuros modelos con estructura definida

Componentes:

Clase Base: BaseModel
Interface común para todos los modelos
Métodos abstractos: predict(), validate_input(), get_risk_category()
Clase MockModel1 y MockModel2
Implementaciones mock que retornan valores aleatorios realistas
Documentación de inputs/outputs esperados
Estructura preparada para implementación futura
Documentación
Comentarios detallados sobre qué tipo de modelos se esperan
Ejemplos de cómo implementar nuevos modelos
Guías de integración
Sub-paso 1.6: Procesador de Datos (utils/data_processor.py)
Objetivo: Centralizar lógica de procesamiento y validación de datos

Componentes:

Clase DataProcessor
Validación: Rangos válidos, tipos de datos, valores requeridos
Transformación: Normalización, encoding de categóricas
Preparación: Conversión a formatos requeridos por modelos
Funciones de Utilidad
Conversión de datos de UI a formato de modelo
Manejo de valores faltantes
Cálculo de estadísticas descriptivas
Integración con Código Existente
Utilizar ClinicalProcessor del código existente
Adaptar funciones de preprocesamiento de GNN
Mantener compatibilidad con estructura de datos actual

Sub-paso 1.7: Visualizaciones (utils/visualization.py)
Objetivo: Funciones para generar gráficos interactivos con Bokeh

Componentes:

Gráficos de Resultados
risk_distribution_plot(): Histograma con posición del paciente
survival_curves_plot(): Curvas de Kaplan-Meier interactivas
risk_stratification_table(): Tabla dinámica con highlighting
Gráficos de Comparación
model_comparison_plot(): Comparación de scores de modelos
hazard_ratio_plot(): Visualización de hazard ratios
Utilidades de Bokeh
Paletas de colores consistentes
Estilos y temas personalizados
Funciones helper para layouts responsivos

Sub-paso 1.8: Archivos Estáticos (static/)
Objetivo: Recursos adicionales para la aplicación

Archivos:

custom.css: Estilos personalizados para Bokeh
logo.png: Logo de la aplicación (si aplica)
survival_data.json: Datos de supervivencia por categoría de riesgo (basados en literatura)


Sub-paso 1.9: Configuración de Testing (tests/)
Objetivo: Estructura base para testing automatizado

Archivos:

test_config.py: Tests de carga de configuración
test_models.py: Tests de modelos mock
test_data_processor.py: Tests de validación y transformación
conftest.py: Fixtures compartidas para pytest

Sub-paso 1.10: Configuración de Logging y Monitoreo
Objetivo: Sistema de logging para debugging y monitoreo

Componentes:

Configuración de logging en múltiples niveles
Rotación de logs automática
Métricas de uso y performance
Manejo de errores con stack traces
Consideraciones de Implementación
Compatibilidad con Código Existente
Reutilizar al máximo las utilidades existentes en 3.GeneticScores
Mantener compatibilidad con estructura de datos actual
Adaptar funciones de procesamiento sin duplicar código
Escalabilidad
Estructura modular permite fácil adición de nuevos modelos
Configuración externa evita hard-coding
Separación clara entre lógica de negocio y presentación
Performance
Cache de modelos pre-cargados
Validación de datos optimizada
Lazy loading de componentes pesados
Mantenibilidad
Documentación exhaustiva en código
Separación clara de responsabilidades
Tests automatizados desde el inicio
Esta propuesta establece una base sólida y extensible que mantiene compatibilidad con el código existente mientras proporciona la flexibilidad necesaria para el desarrollo futuro de la aplicación.