"""
Aplicación principal para MDS Risk Calculator con Bokeh.

Este módulo implementa la aplicación web principal usando Bokeh Server para
calcular el riesgo de pacientes con Síndrome Mielodisplásico (MDS).
"""

import logging
import time
import yaml
from pathlib import Path
from typing import Dict, Any, Optional
import pandas as pd

from bokeh.application import Application
from bokeh.application.handlers import FunctionHandler
from bokeh.io import curdoc
from bokeh.layouts import column, row
from bokeh.models import (
    Tabs, TabPanel, Div, Button, TextInput, Select, CheckboxGroup, 
    NumericInput, PreText, Spacer
)
from bokeh.plotting import figure
from bokeh.server.server import Server

from models.gnn_model import MDSGNNModel
from models.base_model import PatientData, create_empty_patient_data
from utils.graph_builder import GraphBuilder
from ui.widgets import WidgetFactory
from ui.layout_builder import create_layout
from utils.gene_importance import GeneImportanceEvaluator

# Configuración de logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class MDSCalculatorApp:
    """
    Aplicación principal para cálculo de riesgo MDS.
    
    Maneja el estado de la aplicación, componentes de UI y callbacks
    para crear una interfaz interactiva completa.
    """
    
    def __init__(self, config_path: str = "app/config.yaml"):
        """
        Inicializa la aplicación con configuración y modelos.
        
        Args:
            config_path: Ruta al archivo de configuración YAML
        """
        self.config_path = config_path
        self.config = self._load_config()
        
        self.cyto_ipssr_mapping = {
            'Very Poor': 'Very-Poor',
            'Poor': 'Poor',
            'Intermediate': 'Int',
            'Good': 'Good',
            'Very Good': 'Very-Good'
        }

        # Estado de la aplicación
        self.current_patient_data = create_empty_patient_data()
        self.last_prediction = None
        self.calculation_in_progress = False
        
        # Componentes de UI
        self.widgets = {}
        self.panels = {}
        self.layout = None
        
        # Modelo GNN
        self.gnn_model = None
        
        # Inicializar aplicación
        self._initialize_model()

        ## Define widgets
        self.widgets = WidgetFactory(self.config).create_all_widgets()
        self.layout = create_layout(self.widgets)
        #self._create_widgets()
        #self._create_layout()

        logger.info("MDSCalculatorApp inicializada correctamente")
    
    def _load_config(self) -> Dict[str, Any]:
        """Carga configuración desde archivo YAML."""
        try:
            config_path = Path(self.config_path)
            if not config_path.exists():
                logger.warning(f"Archivo de configuración no encontrado: {config_path}")
                
                
            with open(config_path, 'r', encoding='utf-8') as f:
                config = yaml.safe_load(f)
            
            logger.info(f"Configuración cargada desde {config_path}")
            return config
            
        except Exception as e:
            logger.error(f"Error cargando configuración: {e}") 
    
    def _initialize_model(self):
        """Inicializa el modelo GNN con manejo de errores."""
        try:
            logger.info("Inicializando modelo GNN...")
            self.gnn_model = MDSGNNModel(self.config_path)
            logger.info("Modelo GNN inicializado correctamente")
            
        except Exception as e:
            logger.error(f"Error inicializando modelo GNN: {e}")
            logger.warning("Continuando sin modelo GNN - solo validación disponible")
            self.gnn_model = None
    
    def setup_callbacks(self):
        """Configura todos los callbacks de la aplicación."""
        # Callback principal de cálculo
        self.widgets['calculate_button'].on_click(self._calculate_risk)
        
        # Callbacks de validación en tiempo real
        for widget_name in ['AGE', 'BM_BLAST', 'WBC', 'ANC', 'HB', 'PLT', 'MONOCYTES']:
            self.widgets[widget_name].on_change('value', self._validate_inputs)

        for widget_name in ['SEX', 'CYTO_IPSSR']:
            self.widgets[widget_name].on_change('value', self._validate_inputs)
            
        self.widgets['mutations'].on_change('active', self._validate_inputs)
        self.widgets['karyotype'].on_change('active', self._validate_inputs)
    
    def _validate_inputs(self, attr, old, new):
        """Valida inputs en tiempo real y actualiza estado."""
        try:
            # Verificar campos requeridos
            missing_fields = self._check_required_fields()
            
            if missing_fields:
                self.widgets['status_message'].text = f"Campos requeridos vacíos: {', '.join(missing_fields)}"
                self.widgets['calculate_button'].disabled = True
                return
            
            # Actualizar datos del paciente
            self._update_patient_data()
            
            # Validar con el modelo si está disponible
            if self.gnn_model:
                validation = self.gnn_model.validate_input(self.current_patient_data)
                
                if validation.errors:
                    self.widgets['status_message'].text = f"Errores: {', '.join(validation.errors)}"
                    self.widgets['calculate_button'].disabled = True
                elif validation.warnings:
                    self.widgets['status_message'].text = f"Advertencias: {', '.join(validation.warnings)}"
                    self.widgets['calculate_button'].disabled = False
                else:
                    self.widgets['status_message'].text = "Datos válidos - Listo para calcular"
                    self.widgets['calculate_button'].disabled = False
            else:
                self.widgets['status_message'].text = "Datos completos - Modelo no disponible (solo validación básica)"
                self.widgets['calculate_button'].disabled = False  # Permitir cálculo mock
                
        except Exception as e:
            logger.error(f"Error en validación: {e}")
            self.widgets['status_message'].text = f"Error en validación: {str(e)}"
            self.widgets['calculate_button'].disabled = True
    
    def _update_patient_data(self):
        """Actualiza los datos del paciente desde los widgets."""
        # Verificar que todos los campos numéricos tengan valores
        numeric_fields = ['AGE', 'BM_BLAST', 'WBC', 'ANC', 'HB', 'PLT', 'MONOCYTES']
        
        for field in numeric_fields:
            if self.widgets[field].value is None:
                # Si algún campo está vacío, no actualizar datos y deshabilitar cálculo
                self.widgets['status_message'].text = f"Campo requerido vacío: {field}"
                self.widgets['calculate_button'].disabled = True
                return
        
        # Obtener mutaciones seleccionadas
        supported_genes = self.config['data']['gene_embeddings']['supported_genes']
        selected_mutations = [supported_genes[i] for i in self.widgets['mutations'].active]
        
        # Actualizar PatientData solo si todos los campos están completos
        self.current_patient_data = PatientData(
            AGE=float(self.widgets['AGE'].value),
            BM_BLAST=float(self.widgets['BM_BLAST'].value),
            WBC=float(self.widgets['WBC'].value),
            ANC=float(self.widgets['ANC'].value),
            HB=float(self.widgets['HB'].value),
            PLT=float(self.widgets['PLT'].value),
            MONOCYTES=float(self.widgets['MONOCYTES'].value),
            SEX=self.widgets['SEX'].value,
            CYTO_IPSSR=self.cyto_ipssr_mapping[self.widgets['CYTO_IPSSR'].value],
            plus8=0 in self.widgets['karyotype'].active,
            del7=1 in self.widgets['karyotype'].active,
            del20q=2 in self.widgets['karyotype'].active,
            del7q=3 in self.widgets['karyotype'].active,
            delY=4 in self.widgets['karyotype'].active,
            del5q=5 in self.widgets['karyotype'].active,
            complex='complex' if 6 in self.widgets['karyotype'].active else 'non-complex',
            mutations=selected_mutations
        )
    
    def _calculate_risk(self):
        """Calcula el riesgo usando el modelo GNN."""
        if self.calculation_in_progress:
            return
        
        # Verificar que todos los campos estén completos antes de proceder
        missing_fields = self._check_required_fields()
        if missing_fields:
            self.widgets['status_message'].text = f"No se puede calcular: campos vacíos - {', '.join(missing_fields)}"
            return
            
        self.calculation_in_progress = True
        self.widgets['calculate_button'].disabled = True
        self.widgets['status_message'].text = "Calculando predicción de riesgo..."
        
        try:
            # Actualizar datos del paciente
            self._update_patient_data()
            
            if not self.gnn_model:
                self._show_error("Modelo GNN no disponible")
                return
            
            # Realizar predicción
            start_time = time.time()
            prediction = self.gnn_model.predict(self.current_patient_data)
            calculation_time = time.time() - start_time
            
            # Guardar resultado
            self.last_prediction = prediction
            
            # Actualizar UI con resultados
            self._update_results_display(prediction, calculation_time)
            
            self.widgets['status_message'].text = f"Cálculo completado en {calculation_time:.2f}s"
            
        except Exception as e:
            logger.error(f"Error en cálculo de riesgo: {e}")
            self._show_error(f"Error calculando riesgo: {str(e)}")
            
        finally:
            self.calculation_in_progress = False
            self.widgets['calculate_button'].disabled = False
    
    def _check_required_fields(self) -> list:
        """Verifica qué campos requeridos están vacíos."""
        missing_fields = []
        
        # Campos numéricos requeridos
        numeric_fields = {
            'AGE': 'Edad',
            'BM_BLAST': '% Blastos',
            'WBC': 'WBC',
            'ANC': 'ANC',
            'HB': 'Hemoglobina',
            'PLT': 'Plaquetas',
            'MONOCYTES': 'Monocitos'
        }
        
        for field_key, field_name in numeric_fields.items():
            if self.widgets[field_key].value is None:
                missing_fields.append(field_name)
        
        # Los campos categóricos siempre tienen valor por defecto, no necesitan verificación
        
        return missing_fields
    
    def _update_results_display(self, prediction, calculation_time):
        """Actualiza la visualización de resultados."""
        # Actualizar scores principales
        self.widgets['risk_score'].text = f"<h3>Score de Riesgo: {prediction.raw_score:.3f}</h3>"
        
        # Obtener color de categoría
        if self.gnn_model:
            category_label = self.gnn_model.get_risk_category_label(prediction.risk_category)
            category_color = self.gnn_model.get_risk_category_color(prediction.risk_category)
        else:
            category_label = prediction.risk_category
            category_color = "#666666"
            
        self.widgets['risk_category'].text = (
            f"<h3 style='color: {category_color}'>Categoría: {category_label}</h3>"
        )
                # Actualizar detalles
        details_text = f"""Detalles de la Predicción:

Score de Riesgo: {prediction.raw_score:.4f}
Categoría: {prediction.risk_category} ({category_label})

Modelo: {prediction.model_version}
Tiempo de Cálculo: {calculation_time:.3f} segundos

Datos del Paciente:
- Edad: {self.current_patient_data.AGE} años
- Blastos en MO: {self.current_patient_data.BM_BLAST}%
- WBC: {self.current_patient_data.WBC} ×10⁹/L
- Mutaciones: {', '.join(self.current_patient_data.mutations) if self.current_patient_data.mutations else 'Ninguna'}
"""
        if prediction.warnings:
            details_text += f"\nAdvertencias:\n" + "\n".join(f"- {w}" for w in prediction.warnings)
            
        self.widgets['hazard_ratio'].text = f"<h3>Hazard Ratio: {prediction.hazard_ratio:.3f}</h3>"

        # Actualizar gráfico simple
        self._update_risk_plot(prediction)
    
    def _update_risk_plot(self, prediction):
        """Actualiza el gráfico de distribución de riesgo."""
        # Limpiar gráfico anterior
        self.widgets['risk_plot'].renderers = []
        
        # Crear histograma simple de ejemplo
        import numpy as np
        from scipy.stats import gaussian_kde

        # Generar datos de ejemplo para distribución
        cohort_scores = self._load_cohort_scores()

        gene_evaluator = GeneImportanceEvaluator(self.gnn_model)
        shap_dict = gene_evaluator.evaluate_gene_importance(self.current_patient_data, self.gnn_model.gene_embeddings)
        shap_values = shap_dict.get('shap_values', {})
        base_score = shap_dict.get('base', {})
        if cohort_scores:
            # Crear histograma con los datos reales
            kde = gaussian_kde(cohort_scores)
            
            # Generar puntos suaves para la curva
            x_min, x_max = min(cohort_scores), max(cohort_scores)
            x_range = x_max - x_min
            x_vals = np.linspace(x_min - 0.1 * x_range, x_max + 0.1 * x_range, 200)
            y_vals = kde(x_vals)
            max_y = max(y_vals)
                
            # Dibujar histograma como línea
            self.widgets['risk_plot'].line(x_vals, y_vals, line_width=2, color='black', alpha=0.7)

            # Añadir líneas verticales para categorías de riesgo
            if 'risk_categories' in self.config and 'thresholds' in self.config['risk_categories']:
                thresholds = self.config['risk_categories']['thresholds']
                colors = self.config['risk_categories'].get('colors', {})
                labels = self.config['risk_categories'].get('labels_short', {})
                
                prev_threshold = x_min - 0.05 * x_range
                for category, threshold in thresholds.items():

                    if isinstance(threshold, (int, float)) and x_min <= threshold <= x_max:
                        # Color específico para cada categoría
                        segment_color = line_color = colors.get(category, '#666666')

                        ## Colorear segmento de la curva
                        mask = (x_vals >= prev_threshold) & (x_vals <= threshold)
                        x_segment = x_vals[mask]
                        y_segment = y_vals[mask]
                        x_fill = np.concatenate([x_segment, x_segment[::-1]])
                        y_fill = np.concatenate([y_segment, np.zeros_like(y_segment)])
                        
                        self.widgets['risk_plot'].patch(
                            x_fill, y_fill,
                            alpha=0.3, color=segment_color
                        )

                        # Añadir linea vertical
                        line_y = kde(threshold)[0]

                        if category != "very_high":
                            # Línea vertical
                            self.widgets['risk_plot'].line(
                                [threshold, threshold], [0, line_y],
                                line_width=2, color=line_color, alpha=0.8, line_dash='solid'
                            )
                        segment_center = np.mean([prev_threshold, threshold])
                        # Etiqueta de la categoría
                        category_label = labels.get(category, category)
                        self.widgets['risk_plot'].text(
                            x=[segment_center], y=[max_y], 
                            text=[category_label],
                            text_font_size="9pt", text_align="center",
                            text_color=line_color
                        )
                        prev_threshold = threshold                
            # Marcar posición del paciente
            patient_score = prediction.raw_score
            # Encontrar altura aproximada en la distribución
            patient_y = kde(patient_score)[0]  # Altura en la curva de densidad
                
            self.widgets['risk_plot'].circle(
                [patient_score], [patient_y], 
                size=15, color='red', alpha=0.8
            )

            if len(shap_values) > 0:
                ## Ordenar diccionario de SHAP
                shap_values_sorted = dict(sorted(shap_values.items(), key=lambda item: item[1], reverse=True))

                base_y = kde(base_score)[0]
                self.widgets['risk_plot'].circle(
                    [base_score], [base_y], 
                    size=15, color='grey', alpha=0.8
                )
                self.widgets['risk_plot'].text(
                        x=[base_score], y=[base_y + 0.02], 
                        text=["Base"],
                        text_font_size="9pt", text_align="center",
                        text_color='black'
                    )
                gene_value = base_score
                for gene, value in shap_values_sorted.items():
                    gene_value += value
                    gene_y = kde(gene_value)[0]
                    self.widgets['risk_plot'].circle(
                        [gene_value], [gene_y], 
                        size=15, color='pink', alpha=0.8
                    )
                    self.widgets['risk_plot'].text(
                        x=[gene_value], y=[gene_y + 0.02], 
                        text=[gene],
                        text_font_size="9pt", text_align="center",
                        text_color='black'
                    )


            # Calcular percentil del paciente
            percentile = (np.sum(np.array(cohort_scores) <= patient_score) / len(cohort_scores)) * 100
                
            # Configurar labels
            self.widgets['risk_plot'].xaxis.axis_label = "Score de Riesgo"
            self.widgets['risk_plot'].yaxis.axis_label = ""
            self.widgets['risk_plot'].title.text = f"Percentil {percentile:.1f}%"
            
            # Configurar rangos del gráfico
            self.widgets['risk_plot'].x_range.start = x_min - 0.05 * x_range
            self.widgets['risk_plot'].x_range.end = x_max + 0.05 * x_range
            self.widgets['risk_plot'].y_range.start = 0
            self.widgets['risk_plot'].y_range.end = max_y * 1.2  # Espacio para etiquetas

        else:
            # Fallback si no hay datos
            self.widgets['risk_plot'].text(
                x=0, y=0.5, text=["No hay datos de cohorte disponibles"], 
                text_font_size="12pt", text_align="center"
            )
    
    def _load_cohort_scores(self) -> list:
        """Carga los scores reales de la cohorte desde archivo."""
        try:
            scores_file = Path(self.config["risk_categories"]["risk_distribution"]["file_path"])
            if scores_file.exists():
                with open(scores_file, 'r') as f:
                    data = pd.read_csv(scores_file, sep='\t')
                    return data['Score'].tolist()
            else:
                logger.warning("Archivo de scores de cohorte no encontrado")
                return []
                
        except Exception as e:
            logger.error(f"Error cargando scores de cohorte: {e}")
            return []


    def _show_error(self, message: str):
        """Muestra mensaje de error en la interfaz."""
        self.widgets['status_message'].text = f"ERROR: {message}"
        self.widgets['risk_score'].text = "<h3 style='color: red'>Score de Riesgo: ERROR</h3>"
        self.widgets['risk_category'].text = "<h3 style='color: red'>Categoría: ERROR</h3>"
        self.widgets['hazard_ratio'].text = "<h3 style='color: red'>Hazard Ratio: ERROR</h3>"
        self.widgets['prediction_details'].text = f"Error en el cálculo: {message}"


def create_application(config_path: str = "app/config.yaml") -> Application:
    """
    Factory function para crear la aplicación Bokeh.
    
    Args:
        config_path: Ruta al archivo de configuración
        
    Returns:
        Application: Aplicación Bokeh configurada
    """
    def create_doc(doc):
        """Función para crear el documento de Bokeh."""
        try:
            # Crear instancia de la aplicación
            app = MDSCalculatorApp(config_path)
            
            # Configurar callbacks
            app.setup_callbacks()
            
            # Agregar layout al documento
            doc.add_root(app.layout)
            
            # Configurar título
            doc.title = app.config['app'].get('title', 'MDS Risk Calculator')
            
            logger.info("Aplicación Bokeh creada correctamente")
            
        except Exception as e:
            logger.error(f"Error creando aplicación: {e}")
            # Crear layout de error
            error_layout = column(
                Div(text="<h1 style='color: red'>Error de Inicialización</h1>"),
                Div(text=f"<p>Error: {str(e)}</p>"),
                Div(text="<p>Revise la configuración y los logs para más detalles.</p>")
            )
            doc.add_root(error_layout)
            doc.title = "MDS Calculator - Error"
    
    return Application(FunctionHandler(create_doc))


def run_server(config_path: str = "app/config.yaml", 
               port: Optional[int] = None,
               allow_websocket_origin: Optional[list] = None):
    """
    Ejecuta el servidor Bokeh.
    
    Args:
        config_path: Ruta al archivo de configuración
        port: Puerto del servidor (opcional)
        allow_websocket_origin: Orígenes permitidos para WebSocket
    """
    try:
        # Cargar configuración
        with open(config_path, 'r', encoding='utf-8') as f:
            config = yaml.safe_load(f)
        
        server_config = config.get('app', {}).get('server', {})
        
        # Configuración del servidor
        server_port = port or server_config.get('port', 5006)
        allowed_origins = allow_websocket_origin or server_config.get('allow_websocket_origin', [f"localhost:{server_port}"])
        
        # Crear aplicación
        app = create_application(config_path)
        
        # Configurar servidor
        server = Server({'/': app}, 
                       port=server_port,
                       allow_websocket_origin=allowed_origins,
                       num_procs=server_config.get('num_procs', 1))
        
        logger.info(f"Iniciando servidor en puerto {server_port}")
        logger.info(f"Aplicación disponible en: http://localhost:{server_port}")
        
        # Iniciar servidor
        server.start()
        
        # Mantener servidor corriendo
        if __name__ == '__main__':
            server.io_loop.start()
        
    except Exception as e:
        logger.error(f"Error iniciando servidor: {e}")
        raise


# Punto de entrada principal
if __name__ == '__main__':
    import sys
    
    # Configurar logging para producción
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler('logs/mds_calculator.log'),
            logging.StreamHandler(sys.stdout)
        ]
    )
    
    # Ejecutar servidor
    try:
        run_server()
    except KeyboardInterrupt:
        logger.info("Servidor detenido por el usuario")
    except Exception as e:
        logger.error(f"Error fatal: {e}")
        sys.exit(1)
