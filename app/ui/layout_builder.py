# -*- coding: utf-8 -*-
"""
Definición de layout para la aplicación MDS Calculator.

Este módulo contiene la función para crear el layout principal de la aplicación,
basado exactamente en el código original sin funcionalidad adicional.
"""

from bokeh.layouts import column, row
from bokeh.models import Tabs, TabPanel, Div, Spacer


def create_layout(widgets):
    """
    Crea el layout principal de la aplicación.
    
    Args:
        widgets: Diccionario con todos los widgets de la aplicación
        
    Returns:
        Layout principal de Bokeh
    """
    # Panel de entrada - Datos Clínicos
    clinical_panel = TabPanel(
        child=column(
            row(widgets['AGE'], widgets['BM_BLAST']),
            row(widgets['WBC'], widgets['ANC']),
            row(widgets['HB'], widgets['PLT']),
            row(widgets['MONOCYTES'], widgets['SEX']),
            row(widgets['CYTO_IPSSR'])
        ),
        title="Datos Clínicos"
    )
    
    # Panel de entrada - Mutaciones
    mutation_panel = TabPanel(
        child=column(
            Div(text="<h4>Seleccione las mutaciones presentes:</h4>"),
            widgets['mutations']
        ),
        title="Mutaciones"
    )
    
    # Panel de entrada - Cariotipo
    karyotype_panel = TabPanel(
        child=column(
            Div(text="<h4>Seleccione las aberraciones cromosómicas:</h4>"),
            widgets['karyotype']
        ),
        title="Cariotipo"
    )
    
    # Pestañas de entrada
    input_tabs = Tabs(tabs=[clinical_panel, mutation_panel, karyotype_panel])
    
    # Panel de resultados - Scores
    scores_panel = TabPanel(
        child=column(
            widgets['risk_score'],
            widgets['risk_category'],
          #  widgets['hazard_ratio'],
            widgets['risk_plot']
        ),
        title="Scores"
    )
    
    # Panel de resultados - Detalles
    details_panel = TabPanel(
        child=column(
            widgets['prediction_details']
        ),
        title="Detalles"
    )
    
    # Panel de resultados - Visualización
    viz_panel = TabPanel(
        child=column(
            widgets['risk_plot']
        ),
        title="Visualización"
    )
    
    # Pestañas de resultados
   # result_tabs = Tabs(tabs=[scores_panel, details_panel, viz_panel])
    result_tabs = Tabs(tabs=[scores_panel])

    # Layout principal
    left_panel = column(
        Div(text="<h2>Datos del Paciente</h2>"),
        input_tabs,
        Spacer(height=20),
        widgets['calculate_button'],
        widgets['status_message'],
        width=450
    )
    
    right_panel = column(
        Div(text="<h2>Resultados</h2>"),
        result_tabs,
        width=550
    )
    
    layout = column(
        widgets['title'],
        row(left_panel, Spacer(width=50), right_panel)
    )
    
    return layout
