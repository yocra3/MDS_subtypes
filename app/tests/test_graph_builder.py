"""
Tests para el módulo GraphBuilder.

Valida la construcción de grafos heterogéneos para modelos GNN de predicción de riesgo MDS.
"""

import pytest
import torch
import pandas as pd
import numpy as np
from unittest.mock import Mock, patch, MagicMock
from torch_geometric.data import HeteroData
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import StandardScaler, OneHotEncoder

from app.utils.graph_builder import GraphBuilder, create_dummy_gene_embeddings, load_gene_embeddings_from_file
from app.models.base_model import PatientData
import tempfile


class TestGraphBuilder:
    """Test suite para la clase GraphBuilder."""
    
    @pytest.fixture
    def mock_config(self):
        """Configuración mock para testing."""
        return {
            'genes': ['TP53', 'ASXL1', 'SF3B1', 'RUNX1', 'TET2'],
            'embedding_dim': 16
        }
    
    @pytest.fixture
    def mock_transformer(self):
        """Mock ColumnTransformer para testing."""
        transformer = Mock(spec=ColumnTransformer)
        # Simular transformación que devuelve array 2D
        transformer.transform.return_value = np.array([[0.1, 0.2, 0.3, 0.4, 0.5]])
        return transformer
    
    @pytest.fixture
    def sample_patient_data(self):
        """Datos de paciente de ejemplo para testing."""
        return PatientData(
            age=65,
            sex="M",
            bm_blast=8.5,
            wbc=3.2,
            anc=1.5,
            hb=9.8,
            plt=85000,
            monocytes=0.8,
            plus8=False,
            del7=True,
            del20q=False,
            del7q=False,
            delY=False,
            del5q=False,
            complex=False,
            cyto_ipssr="Intermediate",
            mutations=['TP53', 'ASXL1']
        )
    
    @pytest.fixture
    def sample_gene_embeddings(self, mock_config):
        """Embeddings de genes de ejemplo."""
        embedding_dim = mock_config['embedding_dim']
        embeddings_dict = {}
        for gene in mock_config['genes']:
            embeddings_dict[gene] = np.random.randn(embedding_dim)
        embeddings_dict['Gene_0'] = np.zeros(embedding_dim)
        return embeddings_dict
    
    @pytest.fixture
    def graph_builder(self, mock_config, mock_transformer):
        """Instancia de GraphBuilder para testing."""
        builder = GraphBuilder(config=mock_config, transformer=mock_transformer)
        return builder
    
    def test_init_default(self):
        """Test inicialización con parámetros por defecto."""
        builder = GraphBuilder()
        assert builder.config == {}
        assert builder.transformer is None
    
    def test_init_with_config(self, mock_config, mock_transformer):
        """Test inicialización con configuración y transformer."""
        builder = GraphBuilder(config=mock_config, transformer=mock_transformer)
        assert builder.config == mock_config
        assert builder.transformer == mock_transformer
    
    def test_build_patient_graph_success(self, graph_builder, sample_gene_embeddings):
        """Test construcción exitosa de grafo de paciente."""
        # Datos de entrada
        clinical_features = torch.tensor([[0.1, 0.2, 0.3, 0.4, 0.5]])
        mutated_genes = ['TP53', 'ASXL1']
        
        # Ejecutar
        graph = graph_builder.build_patient_graph(
            clinical_features, mutated_genes, sample_gene_embeddings
        )

        # Verificaciones
        assert isinstance(graph, HeteroData)
        assert 'patient' in graph.node_types
        assert 'gene' in graph.node_types
        assert graph['patient'].x.shape[0] == 1
        assert graph['gene'].x.shape[0] == len(mutated_genes)
        assert graph['gene'].gene_names == mutated_genes
    
    def test_build_patient_graph_no_mutations(self, graph_builder, sample_gene_embeddings):
        """Test construcción de grafo sin mutaciones."""
        clinical_features = torch.tensor([[0.1, 0.2, 0.3, 0.4, 0.5]])
        mutated_genes = []
        
        graph = graph_builder.build_patient_graph(
            clinical_features, mutated_genes, sample_gene_embeddings
        )
        
        assert graph['gene'].x.shape[0] == 1
        assert graph['gene'].gene_names == ['Gene_0']  
    
    def test_build_patient_graph_from_data(self, graph_builder, sample_patient_data, sample_gene_embeddings):
        """Test construcción de grafo desde PatientData."""
        with patch.object(graph_builder, 'build_patient_graph') as mock_build:
            mock_graph = Mock(spec=HeteroData)
            mock_build.return_value = mock_graph
            
            result = graph_builder.build_patient_graph_from_data(
                sample_patient_data, sample_gene_embeddings
            )
            
            # Verificar que se llamó build_patient_graph con argumentos correctos
            mock_build.assert_called_once()
            args = mock_build.call_args[0]
            assert args[1] == sample_patient_data.mutations  # mutated_genes
            assert args[2] == sample_gene_embeddings  # gene_embeddings (ahora dict)
            assert result == mock_graph
    
    def test_get_gene_embeddings_valid_genes(self, graph_builder, sample_gene_embeddings):
        """Test obtención de embeddings con genes válidos."""
        mutated_genes = ['TP53', 'ASXL1']
        
        gene_embeddings = graph_builder._get_gene_embeddings(mutated_genes, sample_gene_embeddings)
        
        # Verificar que se obtuvieron embeddings para los genes solicitados
        assert gene_embeddings.shape[0] == len(mutated_genes)
        assert gene_embeddings.shape[1] == 16  # embedding_dim
        assert isinstance(gene_embeddings, torch.Tensor)
    
    def test_get_gene_embeddings_invalid_genes(self, graph_builder, sample_gene_embeddings, caplog):
        """Test obtención de embeddings con genes inválidos."""
        mutated_genes = ['TP53', 'INVALID_GENE', 'ASXL1']
        
        gene_embeddings = graph_builder._get_gene_embeddings(mutated_genes, sample_gene_embeddings)
        
        # Debe crear embeddings para todos los genes (válidos e inválidos)
        assert gene_embeddings.shape[0] == len(mutated_genes) - 1

        # Verificar warning log
        assert "Gen INVALID_GENE no encontrado en embeddings" in caplog.text
    
    def test_create_edges_with_genes(self, graph_builder):
        """Test creación de edges con genes mutados."""
        num_genes = 3
        
        edge_index = graph_builder._create_edges(num_genes)
        
        # Verificar estructura de edges bidireccionales
        assert edge_index.shape == (2, num_genes)  
        assert edge_index.dtype == torch.long
        
        # Verificar edges
        expected_edges = torch.tensor([[0, 1, 2], [0, 0, 0]])
        assert torch.equal(edge_index, expected_edges)
    
    def test_create_edges_no_genes(self, graph_builder):
        """Test creación de edges sin genes mutados."""
        edge_index = graph_builder._create_edges(0)
        
        assert edge_index.shape == (2, 0)
        assert edge_index.dtype == torch.long
    
    def test_create_hetero_data(self, graph_builder):
        """Test creación de objeto HeteroData."""
        patient_x = torch.tensor([[0.1, 0.2, 0.3]])
        gene_x = torch.tensor([[1.0, 1.1], [2.0, 2.1]])
        edge_index = torch.tensor([[0, 0], [0, 1]])
        mutated_genes = ['TP53', 'ASXL1']
        
        graph = graph_builder._create_hetero_data(
            patient_x, gene_x, edge_index, mutated_genes
        )
        
        # Verificar estructura del grafo
        assert isinstance(graph, HeteroData)
        assert torch.equal(graph['patient'].x, patient_x)
        assert torch.equal(graph['gene'].x, gene_x)
        assert graph['gene'].gene_names == mutated_genes
        assert torch.equal(graph['gene', 'patient'].edge_index, edge_index)
    
    def test_patient_data_to_tensor(self, graph_builder, sample_patient_data):
        """Test conversión de PatientData a tensor."""
        result = graph_builder._patient_data_to_tensor(sample_patient_data)
        
        # Verificar que se llamó al transformer
        graph_builder.transformer.transform.assert_called_once()
        
        # Verificar estructura del DataFrame pasado al transformer
        call_args = graph_builder.transformer.transform.call_args[0][0]
        assert isinstance(call_args, pd.DataFrame)
        assert call_args.shape[0] == 1
        assert 'age' in call_args.columns
        assert 'bm_blast' in call_args.columns
        assert 'sex' in call_args.columns
        
        # Verificar tensor resultante
        assert isinstance(result, torch.Tensor)
        assert result.dtype == torch.float32
        assert result.shape == (1, 5)  # Según mock_transformer
    
    def test_validate_graph_valid(self, graph_builder):
        """Test validación de grafo válido."""
        # Crear grafo válido
        graph = HeteroData()
        graph['patient'].x = torch.tensor([[1.0, 2.0]])
        graph['gene'].x = torch.tensor([[3.0, 4.0], [5.0, 6.0]])
        graph['gene', 'patient'].edge_index = torch.tensor([[0, 1], [0, 0]])

        result = graph_builder.validate_graph(graph)
        assert result is True
    
    def test_validate_graph_missing_patient_nodes(self, graph_builder):
        """Test validación de grafo sin nodos paciente."""
        graph = HeteroData()
        graph['gene'].x = torch.tensor([[3.0, 4.0]])
        
        result = graph_builder.validate_graph(graph)
        assert result is False
    
    def test_validate_graph_missing_gene_nodes(self, graph_builder):
        """Test validación de grafo sin nodos gene."""
        graph = HeteroData()
        graph['patient'].x = torch.tensor([[1.0, 2.0]])
        
        result = graph_builder.validate_graph(graph)
        assert result is False
    
    def test_validate_graph_empty_genes(self, graph_builder):
        """Test validación de grafo con nodos gene vacíos."""
        graph = HeteroData()
        graph['patient'].x = torch.tensor([[1.0, 2.0]])
        graph['gene'].x = torch.empty((0, 2))
        
        result = graph_builder.validate_graph(graph)
        assert result is False
    
    def test_validate_graph_missing_edges(self, graph_builder):
        """Test validación de grafo sin edges."""
        graph = HeteroData()
        graph['patient'].x = torch.tensor([[1.0, 2.0]])
        graph['gene'].x = torch.tensor([[3.0, 4.0]])
        
        result = graph_builder.validate_graph(graph)
        assert result is False


class TestUtilityFunctions:
    """Tests para funciones de utilidad."""
    
    def test_create_dummy_gene_embeddings(self):
        """Test creación de embeddings dummy."""
        gene_list = ['TP53', 'ASXL1', 'SF3B1']
        embedding_dim = 8
        
        embeddings = create_dummy_gene_embeddings(gene_list, embedding_dim)
        
        # Verificar que es un diccionario
        assert isinstance(embeddings, dict)
        
        # Verificar que contiene todos los genes + Gene_0
        expected_genes = gene_list + ['Gene_0']
        assert set(embeddings.keys()) == set(expected_genes)
        
        # Verificar dimensiones de embeddings
        for gene, embedding in embeddings.items():
            assert embedding.shape == (embedding_dim,)
            assert isinstance(embedding, np.ndarray)
        
        # Verificar que el gene dummy son zeros
        assert np.allclose(embeddings['Gene_0'], np.zeros(embedding_dim))
        
        # Verificar que el resto no son zeros
        for gene in gene_list:
            assert not np.allclose(embeddings[gene], np.zeros(embedding_dim))
    
    def test_create_dummy_gene_embeddings_default_dim(self):
        """Test creación de embeddings dummy con dimensión por defecto."""
        gene_list = ['TP53', 'ASXL1']
        
        embeddings = create_dummy_gene_embeddings(gene_list)
        
        # Verificar que es un diccionario con los genes correctos
        expected_genes = gene_list + ['Gene_0']
        assert set(embeddings.keys()) == set(expected_genes)
        
        # Verificar dimensión por defecto
        for embedding in embeddings.values():
            assert embedding.shape == (16,)  # dim=16 por defecto
    
    @patch('torch.load')
    def test_load_gene_embeddings_from_file_success(self, mock_torch_load):
        """Test carga exitosa de embeddings desde archivo."""
        # Mock data como diccionario
        # Crear un DataFrame de embeddings con los nombres de genes como índice
        mock_embeddings_df = pd.DataFrame(
            np.vstack([
            torch.randn(16).numpy(),
            torch.randn(16).numpy(),
            np.zeros(16)
            ]),
            index=['TP53', 'ASXL1', 'Gene_0']
        )
        # Simular que torch.load devuelve un DataFrame, pero la función espera un archivo TSV.
        # Así que creamos un archivo temporal TSV para testear correctamente.
        file_path = None
        with tempfile.NamedTemporaryFile(mode='w+', suffix='.tsv', delete=False) as tmpfile:
            file_path = tmpfile.name
            # Guardar el DataFrame como TSV
            mock_embeddings_df.to_csv(tmpfile, sep='\t', header=True)
        
        result = load_gene_embeddings_from_file(file_path)
               
        # Verificar que es un diccionario con numpy arrays
        assert isinstance(result, dict)
        for gene, embedding in result.items():
            assert isinstance(embedding, np.ndarray)
            assert embedding.shape == (16,)
    
    @patch('torch.load')
    def test_load_gene_embeddings_from_file_error(self, mock_torch_load):
        """Test error al cargar embeddings desde archivo."""
        mock_torch_load.side_effect = FileNotFoundError("File not found")
        
        file_path = "/path/to/nonexistent.pt"
        
        with pytest.raises(FileNotFoundError):
            load_gene_embeddings_from_file(file_path)


class TestIntegration:
    """Tests de integración para flujo completo."""
    
    @pytest.fixture
    def real_transformer(self):
        """ColumnTransformer real para testing de integración."""
        # Crear transformer simple para testing
        numeric_features = ['age', 'bm_blast', 'wbc', 'anc', 'hb', 'plt', 'monocytes']
        categorical_features = ['sex', 'cyto_ipssr']
        binary_features = ['plus8', 'del7', 'del20q', 'del7q', 'delY', 'del5q', 'complex']
        
        transformer = ColumnTransformer([
            ('num', StandardScaler(), numeric_features),
            ('cat', OneHotEncoder(drop='first', sparse_output=False), categorical_features),
            ('bin', 'passthrough', binary_features)
        ])
        
        # Simular que ya está entrenado
        transformer.fit(pd.DataFrame({
            'age': [65, 70, 55],
            'bm_blast': [5, 15, 8],
            'wbc': [2.5, 4.0, 3.2],
            'anc': [1.2, 2.0, 1.5],
            'hb': [9.0, 11.0, 10.0],
            'plt': [80000, 120000, 100000],
            'monocytes': [0.5, 1.2, 0.8],
            'sex': ['M', 'F', 'M'],
            'cyto_ipssr': ['Good', 'Poor', 'Intermediate'],
            'plus8': [False, True, False],
            'del7': [True, False, True],
            'del20q': [False, False, False],
            'del7q': [False, True, False],
            'delY': [False, False, False],
            'del5q': [False, False, True],
            'complex': [False, True, False]
        }))
        
        return transformer
    
    def test_full_workflow(self, real_transformer):
        """Test flujo completo desde PatientData hasta grafo validado."""
        # Configurar GraphBuilder
        config = {'genes': ['TP53', 'ASXL1', 'SF3B1', 'RUNX1']}
        builder = GraphBuilder(config=config, transformer=real_transformer)
        
        # Datos de paciente
        patient_data = PatientData(
            age=68,
            sex="F",
            bm_blast=12.5,
            wbc=2.8,
            anc=1.3,
            hb=8.5,
            plt=75000,
            monocytes=1.1,
            plus8=False,
            del7=True,
            del20q=False,
            del7q=False,
            delY=False,
            del5q=True,
            complex=False,
            cyto_ipssr="Poor",
            mutations=['TP53', 'SF3B1']
        )
        
        # Embeddings de genes como diccionario
        gene_embeddings = create_dummy_gene_embeddings(config['genes'], 16)
        
        # Ejecutar flujo completo
        graph = builder.build_patient_graph_from_data(patient_data, gene_embeddings)
        
        # Validar resultado
        assert builder.validate_graph(graph)
        assert graph['gene'].x.shape[0] == 3  # TP53, SF3B1 y Gene_0
        assert graph['gene'].gene_names == ['TP53', 'SF3B1', 'Gene_0']
        assert graph['patient'].x.shape[0] == 1
        
        # Verificar que las features clínicas fueron transformadas
        assert graph['patient'].x.shape[1] > 0
