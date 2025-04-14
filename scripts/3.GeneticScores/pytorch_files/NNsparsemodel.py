"""
Pytorch template for creating a sparse NN model
"""
import torch
from torch import nn
from sparseLayer import SparseGO

class sparseSurvivalNN(nn.Module):
    """Create a sparse NN network with a sparse and a full NN layers.

    Args:
    input_shape: An integer indicating number of input features (e.g. genes).
    go_number: An integer indicating number of GO features
    relation_dict: sparsity codification
    hidden_units: An integer indicating the number of hidden units in the second layer.
    output_shape: An integer indicating number of output units.
    """
    def __init__(self, input_shape: int, go_number:int, hidden_units: int, relation_dict: dict, device = None) -> None:
        super().__init__()
        self.sparse = SparseGO(input_shape, go_number, go_number, relation_dict, device = device)
        self.NN = nn.Linear(in_features = go_number + 1, out_features = hidden_units)
        self.out = nn.Linear(in_features = hidden_units, out_features = 1)
        self.decoder = nn.Linear(in_features = go_number, out_features = input_shape)
        self.relu = nn.ReLU()
        self.sigmoid = nn.Sigmoid()

  
    def forward(self, x1: torch.Tensor, x2: torch.Tensor):
        encoded = self.sparse(x1)
        decoded = self.sigmoid(self.decoder(encoded))
        phi = self.out(self.relu(self.NN(torch.cat((x2, self.relu(encoded)), 1))))
        return phi, decoded

    def predict(self, x1: torch.Tensor, x2: torch.Tensor):
        return self.out(self.relu(self.NN(torch.cat((x2, self.relu(self.sparse(x1))), 1))))
    


class sparseSurvivalNN_cluster(nn.Module):
    """Create a sparse NN network with a sparse and a full NN layers.

    Args:
    input_shape: An integer indicating number of input features (e.g. genes).
    go_number: An integer indicating number of GO features
    relation_dict: sparsity codification
    hidden_units: An integer indicating the number of hidden units in the second layer.
    output_shape: An integer indicating number of output units.
    """
    def __init__(self, input_shape: int, go_number:int, hidden_units: int, relation_dict: dict, device = None) -> None:
        super().__init__()
        self.sparse = SparseGO(input_shape, go_number, go_number, relation_dict, device = device)
        self.NN = nn.Linear(in_features = go_number + 17, out_features = hidden_units)
        self.out = nn.Linear(in_features = hidden_units, out_features = 1)
        self.decoder = nn.Linear(in_features = go_number, out_features = input_shape)
        self.relu = nn.ReLU()
        self.sigmoid = nn.Sigmoid()

  
    def forward(self, x1: torch.Tensor, x2: torch.Tensor):
        encoded = self.sparse(x1)
        decoded = self.sigmoid(self.decoder(encoded))
        phi = self.out(self.relu(self.NN(torch.cat((x2, self.relu(encoded)), 1))))
        return phi, decoded

    def predict(self, x1: torch.Tensor, x2: torch.Tensor):
        return self.out(self.relu(self.NN(torch.cat((x2, self.relu(self.sparse(x1))), 1))))