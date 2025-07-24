# -*- coding: utf-8 -*-
"""
Definición de widgets para la aplicación MDS Calculator.

Este módulo contiene las funciones para crear todos los widgets de la interfaz,
incluyendo campos de entrada clínicos, mutaciones, cariotipo y widgets de resultados.
"""

from bokeh.models import (
    Div, Button, NumericInput, Select, CheckboxGroup, 
    PreText
)
from bokeh.plotting import figure


class WidgetFactory:
    """
    Factory class para crear widgets de la aplicación MDS Calculator.
    """
    
    def __init__(self, config: dict):
        """
        Inicializa el factory con la configuración.
        
        Args:
            config: Diccionario de configuración de la aplicación
        """
        self.config = config
        self.widgets = {}
        
        # Mapeo para valores de cariotipo
        self.cyto_ipssr_mapping = {
            'Very Poor': 'Very-Poor',
            'Poor': 'Poor',
            'Intermediate': 'Int',
            'Good': 'Good',
            'Very Good': 'Very-Good'
        }
    
    def create_all_widgets(self) -> dict:
        """
        Crea todos los widgets de la interfaz.
        
        Returns:
            dict: Diccionario con todos los widgets creados
        """
        # Título principal
        self.widgets['title'] = Div(
            text=f"<h1>{self.config['app'].get('title', 'MDS Risk Calculator')}</h1>",
            width=800,
            height=50
        )
        
        # Crear widgets de entrada
        self._create_clinical_widgets()
        self._create_mutation_widgets()
        self._create_karyotype_widgets()
        
        # Crear widgets de resultados
        self._create_result_widgets()
        
        # Botón de cálculo
        self.widgets['calculate_button'] = Button(
            label="Calcular Riesgo",
            button_type="primary",
            width=200,
            height=40
        )
        
        # Mensaje de estado
        self.widgets['status_message'] = Div(
            text="Ingrese los datos del paciente y presione 'Calcular Riesgo'",
            width=600,
            height=30
        )
        
        return self.widgets
    
    def _create_clinical_widgets(self):
        """Crea widgets para datos clínicos."""
        clinical_ranges = self.config.get('validation', {}).get('clinical_ranges', {})
        categorical_values = self.config.get('validation', {}).get('categorical_values', {})
        
        # Valores por defecto si no hay configuración
        default_ranges = {
            'AGE': {'min': 0, 'max': 120},
            'BM_BLAST': {'min': 0, 'max': 100},
            'WBC': {'min': 0, 'max': 100},
            'ANC': {'min': 0, 'max': 50},
            'HB': {'min': 0, 'max': 25},
            'PLT': {'min': 0, 'max': 2000},
            'MONOCYTES': {'min': 0, 'max': 10}
        }
        
        default_categorical = {
            'SEX': ['M', 'F']
        }
        
        # Usar valores por defecto si no están en config
        for field, range_vals in default_ranges.items():
            if field not in clinical_ranges:
                clinical_ranges[field] = range_vals
        
        for field, values in default_categorical.items():
            if field not in categorical_values:
                categorical_values[field] = values
        
        # Variables numéricas
        self.widgets['AGE'] = NumericInput(
            title="Edad (años):", 
            value=65.0,
            low=clinical_ranges['AGE']['min'],
            high=clinical_ranges['AGE']['max'],
            width=150
        )

        self.widgets['BM_BLAST'] = NumericInput(
            title="% Blastos en médula ósea:", 
            value=5.0,
            low=clinical_ranges['BM_BLAST']['min'],
            high=clinical_ranges['BM_BLAST']['max'],
            mode='float',
            width=150
        )

        self.widgets['WBC'] = NumericInput(
            title="WBC (×10⁹/L):", 
            value=3.0,
            low=clinical_ranges['WBC']['min'],
            high=clinical_ranges['WBC']['max'],
            mode='float',
            width=150
        )

        self.widgets['ANC'] = NumericInput(
            title="ANC (×10⁹/L):", 
            value=1.5,
            low=clinical_ranges['ANC']['min'],
            high=clinical_ranges['ANC']['max'],
            mode='float',
            width=150
        )

        self.widgets['HB'] = NumericInput(
            title="Hemoglobina (g/dL):", 
            value=10.0,
            low=clinical_ranges['HB']['min'],
            high=clinical_ranges['HB']['max'],
            mode='float',
            width=150
        )
        
        self.widgets['PLT'] = NumericInput(
            title="Plaquetas (×10⁹/L):", 
            value=150.0,
            low=clinical_ranges['PLT']['min'],
            high=clinical_ranges['PLT']['max'],
            mode='float',
            width=150
        )

        self.widgets['MONOCYTES'] = NumericInput(
            title="Monocitos (×10⁹/L):", 
            value=0.3,
            low=clinical_ranges['MONOCYTES']['min'],
            high=clinical_ranges['MONOCYTES']['max'],
            mode='float',
            width=150
        )
        
        # Variables categóricas
        self.widgets['SEX'] = Select(
            title="Sexo:",
            value='M',
            options=categorical_values['SEX'],
            width=150
        )

        self.widgets['CYTO_IPSSR'] = Select(
            title="Riesgo Citogenético (IPSS-R):",
            value='Intermediate',
            options=list(self.cyto_ipssr_mapping.keys()),
            width=200
        )
    
    def _create_karyotype_widgets(self):
        """Crea widgets para aberraciones cromosómicas."""
        karyotype_options = [
            "Trisomía 8 (+8)",
            "Deleción cromosoma 7 (del7)",
            "Deleción 20q (del20q)",
            "Deleción 7q (del7q)",
            "Deleción cromosoma Y (delY)",
            "Deleción 5q (del5q)",
            "Cariotipo complejo"
        ]
        
        self.widgets['karyotype'] = CheckboxGroup(
            labels=karyotype_options,
            active=[],
            width=300
        )
    
    def _create_mutation_widgets(self):
        """Crea widgets para mutaciones genéticas."""
        supported_genes = self.config.get('data', {}).get('gene_embeddings', {}).get('supported_genes', [
            'ASXL1', 'CBL', 'DNMT3A', 'ETV6', 'EZH2', 'IDH2', 
            'KRAS', 'NRAS', 'RUNX1', 'SF3B1', 'SRSF2', 'STAG2', 
            'TET2', 'TP53', 'U2AF1', 'ZRSR2'
        ])
        
        self.widgets['mutations'] = CheckboxGroup(
            labels=supported_genes,
            active=[],
            width=300,
            height=400
        )
    
    def _create_result_widgets(self):
        """Crea widgets para mostrar resultados."""
        # Resultado principal
        self.widgets['risk_score'] = Div(
            text="<h3>Score de Riesgo: --</h3>",
            width=300,
            height=50
        )
        
        self.widgets['risk_category'] = Div(
            text="<h3>Categoría: --</h3>",
            width=300,
            height=50
        )
        
        self.widgets['hazard_ratio'] = Div(
            text="<h3>Hazard Ratio: --</h3>",
            width=300,
            height=50
        )
        
        # Detalles de la predicción
        self.widgets['prediction_details'] = PreText(
            text="Presione 'Calcular Riesgo' para ver los resultados detallados.",
            width=500,
            height=200
        )
        
        # Placeholder para gráficos
        self.widgets['risk_plot'] = figure(
            title="Distribución de Riesgo",
            width=500,
            height=300,
            tools="pan,zoom_in,zoom_out,reset"
        )
        
        # Plot inicial vacío
        self.widgets['risk_plot'].circle([0], [0], size=0, alpha=0)
    
    def get_widget(self, name: str):
        """
        Obtiene un widget específico por nombre.
        
        Args:
            name: Nombre del widget
            
        Returns:
            Widget de Bokeh o None si no existe
        """
        return self.widgets.get(name)
    
    def get_widgets_by_category(self, category: str) -> dict:
        """
        Obtiene widgets por categoría.
        
        Args:
            category: Categoría ('clinical', 'mutations', 'karyotype', 'results')
            
        Returns:
            dict: Diccionario con widgets de la categoría especificada
        """
        if category == 'clinical':
            return {
                name: widget for name, widget in self.widgets.items()
                if name in ['AGE', 'BM_BLAST', 'WBC', 'ANC', 'HB', 'PLT', 'MONOCYTES', 'SEX', 'CYTO_IPSSR']
            }
        elif category == 'mutations':
            return {
                name: widget for name, widget in self.widgets.items()
                if name == 'mutations'
            }
        elif category == 'karyotype':
            return {
                name: widget for name, widget in self.widgets.items()
                if name == 'karyotype'
            }
        elif category == 'results':
            return {
                name: widget for name, widget in self.widgets.items()
                if name in ['risk_score', 'risk_category', 'hazard_ratio', 'prediction_details', 'risk_plot']
            }
        elif category == 'controls':
            return {
                name: widget for name, widget in self.widgets.items()
                if name in ['calculate_button', 'status_message']
            }
        else:
            return {}


def create_widgets(config: dict) -> dict:
    """
    Función de conveniencia para crear todos los widgets.
    
    Args:
        config: Diccionario de configuración de la aplicación
        
    Returns:
        dict: Diccionario con todos los widgets creados
    """
    factory = WidgetFactory(config)
    return factory.create_all_widgets()


def create_clinical_widgets(config: dict) -> dict:
    """
    Crea únicamente los widgets para datos clínicos.
    
    Args:
        config: Diccionario de configuración de la aplicación
        
    Returns:
        dict: Diccionario con widgets clínicos
    """
    factory = WidgetFactory(config)
    factory._create_clinical_widgets()
    return factory.get_widgets_by_category('clinical')


def create_mutation_widgets(config: dict) -> dict:
    """
    Crea únicamente los widgets para mutaciones.
    
    Args:
        config: Diccionario de configuración de la aplicación
        
    Returns:
        dict: Diccionario con widgets de mutaciones
    """
    factory = WidgetFactory(config)
    factory._create_mutation_widgets()
    return factory.get_widgets_by_category('mutations')


def create_karyotype_widgets(config: dict) -> dict:
    """
    Crea únicamente los widgets para cariotipo.
    
    Args:
        config: Diccionario de configuración de la aplicación
        
    Returns:
        dict: Diccionario con widgets de cariotipo
    """
    factory = WidgetFactory(config)
    factory._create_karyotype_widgets()
    return factory.get_widgets_by_category('karyotype')


def create_result_widgets(config: dict) -> dict:
    """
    Crea únicamente los widgets para resultados.
    
    Args:
        config: Diccionario de configuración de la aplicación
        
    Returns:
        dict: Diccionario con widgets de resultados
    """
    factory = WidgetFactory(config)
    factory._create_result_widgets()
    return factory.get_widgets_by_category('results')
