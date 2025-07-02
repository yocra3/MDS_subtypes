"""
Pytorch template for creating a single NN model
"""
import torch
from torch import nn 

class survivalNN(nn.Module):
    """Create a full NN network with two hidden layers. The first layer will combine
    data from the mutations. The second will concatenate the mutation data with the 
    clinical data

    Args:
    input_shape: An integer indicating number of input features (e.g. genes).
    hidden_units: A tuple indicating number of hidden units in the first layer.
    output_shape: An integer indicating number of output units.
    """
    def __init__(self, input_shape: int, hidden_units: tuple) -> None:
        super().__init__()
        self.NN = nn.Linear(in_features = input_shape,
                    out_features=hidden_units[0])
        self.NN2 = nn.Linear(in_features = hidden_units[0] + 1,
            out_features = hidden_units[1])
        self.out = nn.Linear(in_features = hidden_units[1],
            out_features = 1)
        self.relu = nn.ReLU()
       
    def forward(self, x1: torch.Tensor, x2: torch.Tensor):
        return self.out(self.relu(self.NN2(torch.cat((x2, self.relu(self.NN(x1))), 1))))


