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
        
        #self.widgets['hazard_ratio'].text = f"<h3>Hazard Ratio: {prediction.hazard_ratio:.2f}</h3>"
        
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
            
        self.widgets['prediction_details'].text = details_text
        
        # Actualizar gráfico simple
        self._update_risk_plot(prediction)
    
    def _update_risk_plot(self, prediction):
        """Actualiza el gráfico de distribución de riesgo."""
        # Limpiar gráfico anterior
        self.widgets['risk_plot'].renderers = []
        
        # Crear histograma simple de ejemplo
        import numpy as np
        
        # Generar datos de ejemplo para distribución
        x_vals = np.linspace(-3, 3, 100)
        y_vals = np.exp(-0.5 * x_vals**2) / np.sqrt(2 * np.pi)  # Distribución normal
        
        # Agregar línea de distribución
        self.widgets['risk_plot'].line(x_vals, y_vals, line_width=2, color='blue', alpha=0.7)
        
        # Marcar posición del paciente
        patient_y = np.exp(-0.5 * prediction.raw_score**2) / np.sqrt(2 * np.pi)
        self.widgets['risk_plot'].circle(
            [prediction.raw_score], [patient_y], 
            size=15, color='red', alpha=0.8
        )
        
        # Configurar labels
        self.widgets['risk_plot'].xaxis.axis_label = "Score de Riesgo"
        self.widgets['risk_plot'].yaxis.axis_label = "Densidad"
        self.widgets['risk_plot'].title.text = f"Posición del Paciente (Score: {prediction.raw_score:.3f})"
    
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
