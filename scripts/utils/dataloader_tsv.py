"""
Contains functionality for creating PyTorch DataLoaders for 
loading tabular data.
"""

import torch
from torch.utils.data import Dataset, DataLoader
import pandas as pd

class TabularDataset(Dataset):
    def __init__(self, file_path):
        self.data = pd.read_csv(file_path, sep = "\t")
        
        # Assuming the last column is the output/target
        self.clin = self.data.iloc[:, 0].values
        self.mutations = self.data.iloc[:, 1:-2].values
        self.time = self.data.iloc[:, -2].values
        self.event = self.data.iloc[:, -1].values

    def __len__(self):
        return len(self.data)
    
    def __getitem__(self, idx):
        clin = torch.tensor(self.clin[idx], dtype=torch.float32)
        clin = torch.unsqueeze(clin, dim = -1)

        time = torch.tensor(self.time[idx], dtype=torch.float32)
        time = torch.unsqueeze(time, dim = -1)

        event = torch.tensor(self.event[idx], dtype=torch.float32)
        event = torch.unsqueeze(event, dim = -1)

        return torch.tensor(self.mutations[idx], dtype=torch.float32), clin, time, event
    
    def getFeaturesCols(self):
        return self.mutations.shape[1]


def data_loader_tsv(train_path, test_path, batch_size=32):
    train_dataset = TabularDataset(train_path)
    test_dataset = TabularDataset(test_path)

    train_loader = DataLoader(train_dataset, batch_size = batch_size, shuffle = True)
    test_loader = DataLoader(test_dataset, batch_size = batch_size, shuffle = True)
    
    return train_loader, test_loader, train_dataset.getFeaturesCols()

