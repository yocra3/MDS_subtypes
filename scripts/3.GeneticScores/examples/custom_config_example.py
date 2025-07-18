#!/usr/bin/env python3
"""
Ejemplo de configuración personalizada para el pipeline MDS.

Este ejemplo demuestra:
1. Creación de configuraciones personalizadas desde código
2. Validación de configuraciones antes de ejecución
3. Configuraciones específicas para diferentes casos de uso
"""

import yaml
import json
from pathlib import Path
import logging
from typing import Dict, Any, List, Optional

# Configurar logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class ConfigGenerator:
    """Generador de configuraciones personalizadas para el pipeline MDS."""
    
    def __init__(self, base_data_path: str = "data"):
        """Inicializa el generador de configuraciones."""
        self.base_data_path = Path(base_data_path)
        logger.info(f"🔧 Inicializando generador de configuraciones")
        logger.info(f"📁 Ruta base de datos: {self.base_data_path}")
    
    def create_minimal_config(self) -> Dict[str, Any]:
        """Crea una configuración mínima funcional."""
        
        config = {
            'data_sources': {
                'clinical_data': str(self.base_data_path / 'IPSSMol/df_clinical.tsv')
            },
            
            'variable_processing': {
                'gene_selection': ['TP53', 'ASXL1', 'RUNX1'],  # Solo 3 genes principales
                'continuous_variables': ['AGE'],
                'categorical_variables': [],
                'binary_variables': ['AML_TRANSF'],
                'survival_variables': {
                    'time_variable': 'LFS_YEARS',
                    'status_variable': 'LFS_STATUS'
                }
            },
            
            'cross_validation': {
                'fold_column': 'fold'
            },
            
            'data_generation': {
                'output_formats': ['table'],  # Solo tabla
                'table_format': {
                    'output_filename': 'minimal_table'
                }
            },
            
            'output': {
                'file_prefix': 'minimal'
            }
        }
        
        return config
    
    def create_research_config(self, 
                              genes: List[str],
                              clinical_vars: Dict[str, List[str]],
                              use_embeddings: bool = True,
                              include_mutations: bool = True) -> Dict[str, Any]:
        """
        Crea una configuración para investigación personalizada.
        
        Args:
            genes: Lista de genes a incluir
            clinical_vars: Diccionario con variables clínicas por tipo
            use_embeddings: Si usar embeddings de genes
            include_mutations: Si incluir datos de mutaciones VAF
        """
        
        config = {
            'data_sources': {
                'clinical_data': str(self.base_data_path / 'IPSSMol/df_clinical.tsv')
            },
            
            'variable_processing': {
                'gene_selection': genes,
                'continuous_variables': clinical_vars.get('continuous', []),
                'categorical_variables': clinical_vars.get('categorical', []),
                'binary_variables': clinical_vars.get('binary', []),
                'survival_variables': {
                    'time_variable': 'LFS_YEARS',
                    'status_variable': 'LFS_STATUS'
                }
            },
            
            'cross_validation': {
                'fold_column': 'fold'
            },
            
            'data_generation': {
                'output_formats': ['table', 'graph'],
                'table_format': {
                    'output_filename': 'research_table'
                },
                'graph_format': {
                    'output_filename': 'research_graph'
                }
            },
            
            'output': {
                'file_prefix': 'research'
            }
        }
        
        # Añadir embeddings si se requiere
        if use_embeddings:
            config['data_sources']['gene_embeddings'] = str(self.base_data_path / 'scGPT/gene_embeddings.csv')
        
        # Añadir mutaciones VAF si se requiere
        if include_mutations:
            config['data_sources']['mutation_vaf_data'] = str(self.base_data_path / 'IPSSMol/maf.tsv')
        
        return config
    
    def create_clinical_only_config(self, include_all_vars: bool = False) -> Dict[str, Any]:
        """Crea configuración para análisis solo con variables clínicas."""
        
        if include_all_vars:
            # Incluir todas las variables clínicas disponibles
            continuous_vars = [
                'AGE', 'HB', 'PLT', 'ANC', 'BM_BLAST', 'LDH', 
                'FERRITIN', 'B2_MICRO', 'WBC', 'MONO'
            ]
            categorical_vars = [
                'SEX', 'WHO_SUBTYPE', 'CYTOGENETICS', 'IPSS_R_RISK'
            ]
            binary_vars = [
                'AML_TRANSF', 'PRIOR_MPN', 'PRIOR_CHEMO', 'PRIOR_RT'
            ]
        else:
            # Variables clínicas básicas
            continuous_vars = ['AGE', 'HB', 'PLT', 'ANC', 'BM_BLAST']
            categorical_vars = ['SEX']
            binary_vars = ['AML_TRANSF']
        
        config = {
            'data_sources': {
                'clinical_data': str(self.base_data_path / 'IPSSMol/df_clinical.tsv')
            },
            
            'variable_processing': {
                'gene_selection': [],  # Sin genes
                'continuous_variables': continuous_vars,
                'categorical_variables': categorical_vars,
                'binary_variables': binary_vars,
                'survival_variables': {
                    'time_variable': 'LFS_YEARS',
                    'status_variable': 'LFS_STATUS'
                }
            },
            
            'cross_validation': {
                'fold_column': 'fold'
            },
            
            'data_generation': {
                'output_formats': ['table'],  # Solo tabla para análisis clínico
                'table_format': {
                    'output_filename': 'clinical_table'
                }
            },
            
            'output': {
                'file_prefix': 'clinical_only'
            }
        }
        
        return config
    
    def create_gnn_optimized_config(self, 
                                   gene_panel: str = "comprehensive",
                                   embedding_type: str = "scgpt") -> Dict[str, Any]:
        """
        Crea configuración optimizada para modelos GNN.
        
        Args:
            gene_panel: "minimal", "standard", "comprehensive"
            embedding_type: "scgpt", "onehot"
        """
        
        # Definir paneles de genes
        gene_panels = {
            'minimal': ['TP53', 'ASXL1', 'RUNX1', 'SRSF2', 'SF3B1'],
            'standard': [
                'TP53', 'ASXL1', 'RUNX1', 'SRSF2', 'SF3B1',
                'TET2', 'DNMT3A', 'IDH2', 'CBL', 'EZH2',
                'U2AF1', 'ZRSR2', 'STAG2', 'RAD21', 'SMC3'
            ],
            'comprehensive': [
                'TP53', 'ASXL1', 'RUNX1', 'SRSF2', 'SF3B1',
                'TET2', 'DNMT3A', 'IDH2', 'CBL', 'EZH2',
                'U2AF1', 'ZRSR2', 'STAG2', 'RAD21', 'SMC3',
                'BCOR', 'BCORL1', 'PHF6', 'WT1', 'NPM1',
                'FLT3', 'NRAS', 'KRAS', 'PTPN11', 'JAK2'
            ]
        }
        
        config = {
            'data_sources': {
                'clinical_data': str(self.base_data_path / 'IPSSMol/df_clinical.tsv'),
                'mutation_vaf_data': str(self.base_data_path / 'IPSSMol/maf.tsv')
            },
            
            'variable_processing': {
                'gene_selection': gene_panels[gene_panel],
                'continuous_variables': ['AGE', 'HB', 'PLT', 'ANC', 'BM_BLAST'],
                'categorical_variables': ['SEX'],
                'binary_variables': ['AML_TRANSF'],
                'survival_variables': {
                    'time_variable': 'LFS_YEARS',
                    'status_variable': 'LFS_STATUS'
                }
            },
            
            'cross_validation': {
                'fold_column': 'fold'
            },
            
            'data_generation': {
                'output_formats': ['graph'],  # Solo grafo para GNN
                'graph_format': {
                    'output_filename': f'gnn_{gene_panel}_graph'
                }
            },
            
            'output': {
                'file_prefix': f'gnn_{gene_panel}'
            }
        }
        
        # Añadir embeddings según el tipo
        if embedding_type == "scgpt":
            config['data_sources']['gene_embeddings'] = str(self.base_data_path / 'scGPT/gene_embeddings.csv')
        # Para onehot no se añade source de embeddings
        
        return config
    
    def validate_config(self, config: Dict[str, Any]) -> Dict[str, List[str]]:
        """Valida una configuración antes de su uso."""
        
        logger.info("🔍 Validando configuración...")
        
        issues = {
            'errors': [],
            'warnings': []
        }
        
        # Validar estructura básica
        required_sections = ['data_sources', 'variable_processing', 'data_generation']
        for section in required_sections:
            if section not in config:
                issues['errors'].append(f"Sección requerida faltante: {section}")
        
        # Validar data_sources
        if 'data_sources' in config:
            if 'clinical_data' not in config['data_sources']:
                issues['errors'].append("data_sources.clinical_data es requerido")
            else:
                # Verificar que el archivo existe
                clinical_file = Path(config['data_sources']['clinical_data'])
                if not clinical_file.exists():
                    issues['warnings'].append(f"Archivo clínico no encontrado: {clinical_file}")
        
        # Validar variable_processing
        if 'variable_processing' in config:
            var_proc = config['variable_processing']
            
            # Verificar survival_variables
            if 'survival_variables' not in var_proc:
                issues['errors'].append("variable_processing.survival_variables es requerido")
            else:
                surv_vars = var_proc['survival_variables']
                if 'time_variable' not in surv_vars:
                    issues['errors'].append("survival_variables.time_variable es requerido")
                if 'status_variable' not in surv_vars:
                    issues['errors'].append("survival_variables.status_variable es requerido")
            
            # Verificar selección de genes
            genes = var_proc.get('gene_selection', [])
            if len(genes) == 0:
                issues['warnings'].append("No se seleccionaron genes (análisis solo clínico)")
            elif len(genes) > 50:
                issues['warnings'].append(f"Muchos genes seleccionados ({len(genes)}), podría afectar performance")
        
        # Validar data_generation
        if 'data_generation' in config:
            data_gen = config['data_generation']
            
            if 'output_formats' not in data_gen:
                issues['errors'].append("data_generation.output_formats es requerido")
            else:
                formats = data_gen['output_formats']
                valid_formats = ['table', 'graph']
                invalid_formats = [f for f in formats if f not in valid_formats]
                if invalid_formats:
                    issues['errors'].append(f"Formatos de salida inválidos: {invalid_formats}")
        
        # Resumen de validación
        if issues['errors']:
            logger.error(f"❌ {len(issues['errors'])} errores encontrados")
            for error in issues['errors']:
                logger.error(f"  - {error}")
        
        if issues['warnings']:
            logger.warning(f"⚠️ {len(issues['warnings'])} advertencias encontradas")
            for warning in issues['warnings']:
                logger.warning(f"  - {warning}")
        
        if not issues['errors'] and not issues['warnings']:
            logger.info("✅ Configuración válida")
        
        return issues
    
    def save_config(self, config: Dict[str, Any], output_path: str, 
                   format_type: str = "yaml") -> Path:
        """Guarda una configuración en archivo."""
        
        output_file = Path(output_path)
        output_file.parent.mkdir(parents=True, exist_ok=True)
        
        if format_type.lower() == "yaml":
            if not output_file.suffix:
                output_file = output_file.with_suffix('.yaml')
            
            with open(output_file, 'w', encoding='utf-8') as f:
                yaml.dump(config, f, default_flow_style=False, allow_unicode=True)
        
        elif format_type.lower() == "json":
            if not output_file.suffix:
                output_file = output_file.with_suffix('.json')
            
            with open(output_file, 'w', encoding='utf-8') as f:
                json.dump(config, f, indent=2, ensure_ascii=False)
        
        else:
            raise ValueError(f"Formato no soportado: {format_type}")
        
        logger.info(f"💾 Configuración guardada en: {output_file}")
        return output_file

def demo_config_generation():
    """Demuestra la generación de diferentes configuraciones."""
    
    logger.info("🎯 Iniciando demo de generación de configuraciones...")
    
    # Inicializar generador
    generator = ConfigGenerator()
    
    configs_dir = Path('configs/custom_examples')
    configs_dir.mkdir(parents=True, exist_ok=True)
    
    # 1. Configuración mínima
    logger.info("📝 Generando configuración mínima...")
    minimal_config = generator.create_minimal_config()
    issues = generator.validate_config(minimal_config)
    generator.save_config(minimal_config, configs_dir / 'minimal_example.yaml')
    
    # 2. Configuración de investigación personalizada
    logger.info("📝 Generando configuración de investigación...")
    research_genes = ['TP53', 'ASXL1', 'RUNX1', 'SRSF2', 'SF3B1', 'TET2', 'DNMT3A']
    research_clinical = {
        'continuous': ['AGE', 'HB', 'PLT', 'ANC', 'BM_BLAST', 'LDH'],
        'categorical': ['SEX', 'WHO_SUBTYPE'],
        'binary': ['AML_TRANSF', 'PRIOR_MPN']
    }
    
    research_config = generator.create_research_config(
        genes=research_genes,
        clinical_vars=research_clinical,
        use_embeddings=True,
        include_mutations=True
    )
    generator.validate_config(research_config)
    generator.save_config(research_config, configs_dir / 'research_example.yaml')
    
    # 3. Configuración solo clínica
    logger.info("📝 Generando configuración clínica...")
    clinical_config = generator.create_clinical_only_config(include_all_vars=True)
    generator.validate_config(clinical_config)
    generator.save_config(clinical_config, configs_dir / 'clinical_only_example.yaml')
    
    # 4. Configuración optimizada para GNN
    logger.info("📝 Generando configuración GNN...")
    gnn_config = generator.create_gnn_optimized_config(
        gene_panel="standard",
        embedding_type="scgpt"
    )
    generator.validate_config(gnn_config)
    generator.save_config(gnn_config, configs_dir / 'gnn_optimized_example.yaml')
    
    logger.info("✅ Demo completado exitosamente!")
    logger.info(f"📁 Configuraciones generadas en: {configs_dir}")
    
    return configs_dir

def main():
    """Función principal del ejemplo de configuraciones."""
    
    try:
        # Cambiar al directorio correcto
        script_dir = Path(__file__).parent.parent
        original_dir = Path.cwd()
        
        import os
        os.chdir(script_dir)
        
        logger.info(f"📁 Directorio de trabajo: {Path.cwd()}")
        
        # Ejecutar demo
        configs_dir = demo_config_generation()
        
        logger.info("\n🎉 ¡Generación de configuraciones completada!")
        logger.info(f"📋 Revisa las configuraciones en: {configs_dir}")
        logger.info("💡 Puedes usar estas configuraciones como base para tus análisis")
        
        # Mostrar ejemplo de uso
        print("\n" + "="*60)
        print("EJEMPLO DE USO CON CONFIGURACIONES GENERADAS:")
        print("="*60)
        print("# Usar configuración mínima:")
        print("python prepare_data_unified.py \\")
        print("    --config configs/custom_examples/minimal_example.yaml \\")
        print("    --output results/minimal_test")
        print("")
        print("# Usar configuración de investigación:")
        print("python prepare_data_unified.py \\")
        print("    --config configs/custom_examples/research_example.yaml \\")
        print("    --output results/research_analysis")
        print("="*60)
        
        return 0
        
    except Exception as e:
        logger.error(f"❌ Error en generación de configuraciones: {e}")
        return 1
    
    finally:
        os.chdir(original_dir)

if __name__ == "__main__":
    exit(main())
