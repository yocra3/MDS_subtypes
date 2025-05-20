#Usefull functions for GNN: train, test, Earlystopping, plot
import os

import torch
from torch import nn
import numpy as np
import pandas as pd

import copy

import matplotlib as mpl
import matplotlib.pyplot as plt

from torch_geometric.loader import DataLoader

from sklearn.metrics import roc_auc_score, f1_score, average_precision_score, balanced_accuracy_score

#load the custom splits
def obtain_custom_split(folder_path, Kfold_or_Rfold='Kfold'):

    if Kfold_or_Rfold=='Kfold':
        folder_path=os.path.join(folder_path,'Kfold')
    elif Kfold_or_Rfold=='Rfold':
        folder_path=os.path.join(folder_path,'Rfold')
    else:
        raise NotImplementedError

    # Subfolders names
    subfolders = [name for name in os.listdir(folder_path) if os.path.isdir(os.path.join(folder_path, name))]
    subfolders = sorted(subfolders, key=lambda x: int(x[4:])) #Order the seeds

    custom_cv_splits = []
    test_split=[]

    for folder in subfolders:
        folder_path_subfolder = os.path.join(folder_path, folder)

        folder_path_subfolder_aux = os.path.join(folder_path_subfolder, 'train.csv')
        train_idx = np.loadtxt(folder_path_subfolder_aux, dtype=str)

        folder_path_subfolder_aux = os.path.join(folder_path_subfolder, 'val.csv')
        val_idx = np.loadtxt(folder_path_subfolder_aux, dtype=str)

        custom_cv_splits.append((train_idx, val_idx))

        folder_path_subfolder_aux = os.path.join(folder_path_subfolder, 'test.csv')
        test_idx = np.loadtxt(folder_path_subfolder_aux, dtype=str)

        test_split.append(test_idx) #All equal

    return custom_cv_splits, test_split

#Obtain the data of a split list
def generate_loader_from_split(graphs, train_idx, val_idx, batch_size, device):

    #Train split
    graphs_train = []
    for barcode in train_idx:
        graphs_train.append(next((graph for graph in graphs if graph.barcode == barcode), None))

    #Val split
    graphs_val = []
    for barcode in val_idx:
        graphs_val.append(next((graph for graph in graphs if graph.barcode == barcode), None))

    dataloader_train = DataLoader([graph.to(device) for graph in graphs_train], batch_size=batch_size, shuffle=False)
    dataloader_val = DataLoader([graph.to(device) for graph in  graphs_val], batch_size=batch_size, shuffle=False)

    return dataloader_train, dataloader_val

def train(loader, model, criterion, optimizer):

    model.train()

    loss_mean = 0

    for batch in loader:  # Iterate in batches over the training dataset.
        out = model(batch.x, batch.edge_index, batch.batch)  # Perform a single forward pass.
        y=batch.y

        if str(criterion) == 'BCEWithLogitsLoss()':
            out=out.view(-1)
            y=y.float()

        loss = criterion(out, y)  # Compute the loss.
        loss_mean += loss.cpu().detach().numpy()

        loss.backward()  # Derive gradients.
        optimizer.step()  # Update parameters based on gradients.
        optimizer.zero_grad()  # Clear gradients.

    return loss_mean/len(loader.dataset)

@torch.no_grad()
def test(loader, model, criterion):
    model.eval()

    correct = 0
    auroc = 0
    aupr = 0
    f1 = 0
    balanced_accuracy = 0

    loss_mean = 0

    all_preds = []  # List to store all predictions
    all_labels = []  # List to store all true labels

    for batch in loader:  # Iterate in batches over the training/test dataset.
        out = model(batch.x, batch.edge_index, batch.batch)
        y=batch.y

        if str(criterion)== 'BCEWithLogitsLoss()':
            out=out.view(-1)
            y=y.float()

            loss = criterion(out, y)  # Compute the loss.
            loss_mean += loss.cpu().detach().numpy()
        
            out = out.sigmoid()  # Sigmoid of the output
            pred = (out > 0.5).float()

            correct += int((pred == y).sum())  # Check against ground-truth labels.
            out = out.cpu().detach().numpy()
        
        else:
            loss = criterion(out, y)  # Compute the loss.
            loss_mean += loss.cpu().detach().numpy()

            pred = out.argmax(dim=1)  # Use the class with highest probability.
            correct += int((pred == y).sum())  # Check against ground-truth labels.
            out = out.cpu().detach().numpy()[:,1]
        
        all_labels.extend(y.cpu().detach().numpy())  # True labels
        all_preds.extend(pred.cpu().detach().numpy())  # Predictions

        y = y.cpu().detach().numpy()
        pred_np=pred.cpu().detach().numpy()

        if len(y)!=1 and sum(y)!=0 and sum(y) != len(y):
            # AUROC score
            auroc += roc_auc_score(y, out)

            # AUCPR score
            aupr+=average_precision_score(y,out)

            # F1 score
            f1 += f1_score(y, pred_np) 

            # Balanced accuracy score
            balanced_accuracy += balanced_accuracy_score(y, pred_np)      

    num_data = len(loader.dataset)
    num_batches = len(loader)

    confusion_matrix = pd.crosstab(pd.Series(all_labels), pd.Series(all_preds), rownames=['Real'], colnames=['Predicted'])

    return loss_mean/num_batches, correct / num_data, pred, auroc/num_batches,aupr/num_batches, f1/num_batches, balanced_accuracy/num_batches, confusion_matrix


# class EarlyStopping:
#     # how much (and for how long) the validation loss diverges from the training loss. 
#     # This will break when the validation loss is indeed decreasing but is generally not close enough to the training loss. 
#     def __init__(self, tolerance=5, min_delta=0):

#         self.tolerance = tolerance
#         self.min_delta = min_delta
#         self.counter = 0
#         self.best_model = 0
#         self.val_acc_max = 0
#         self.early_stop = False

#     def __call__(self, train_loss, validation_loss, val_acc, epoch, model):
#         if abs(validation_loss - train_loss) > self.min_delta:
#             self.counter +=1
#             if self.counter >= self.tolerance:  
#                 self.early_stop = True
        
#         if val_acc >= self.val_acc_max:
#             self.best_model = copy.deepcopy(model.state_dict())
#             self.val_acc_max=val_acc
#             self.best_model_epoch=epoch


# 2nd option            
class EarlyStopping:
    #  watch for the trend in validation loss alone, 
    # i.e., if the training is not resulting in lowering of the validation loss then terminate it. 
    def __init__(self, tolerance=1, min_delta=0):
        self.patience = tolerance
        self.min_delta = min_delta
        self.counter = 0
        self.counter_constant = 0
        self.min_validation_loss = float('inf')
        self.best_model = 0
        self.val_acc_max = 0
        self.best_model_epoch=0
        self.early_stop = False
        self.min_cont_delta=0.00001
        self.last_val_loss = 0

    def __call__(self, train_loss, validation_loss, val_acc, epoch, model):
        
        if abs(validation_loss - self.last_val_loss) <= self.min_cont_delta: #Only if loss remain equal for a long period 
            self.counter_constant +=1
            if self.counter_constant >= self.patience:
                self.early_stop = True

        elif validation_loss < self.min_validation_loss:
            self.min_validation_loss = validation_loss 
            self.last_val_loss = validation_loss
            self.counter = 0
            self.counter_constant = 0

        elif validation_loss > (self.min_validation_loss + self.min_delta):
            self.counter += 1
            self.counter_constant = 0
            self.last_val_loss = validation_loss
            if self.counter >= self.patience:
                self.early_stop=True

        else:
            self.counter_constant=0
            self.last_val_loss = validation_loss

        if val_acc >= self.val_acc_max:
            self.best_model = copy.deepcopy(model.state_dict())
            self.val_acc_max=val_acc
            self.best_model_epoch=epoch

# Usage
# early_stopper = EarlyStopper(patience=3, min_delta=10)
#    if early_stopper.early_stop(validation_loss):             
#         break

####### PLOTTING W\O WANDB
mpl.rcParams['figure.max_open_warning'] = 50

def plot(losses_train, losses_val, train_accs, test_accs, output_folder, output_name):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 10))
    # Plotting loss
    ax1.plot(losses_train, c='b', label='Train loss')
    ax1.plot(losses_val, c='g', label='Test loss')
    ax1.set_ylabel('Loss')
    ax1.set_title('Loss train-b, val-g)')
    ax1.legend()

    # Plotting accuracies
    ax2.plot(train_accs, c='b', label='Train Accuracy')
    ax2.plot(test_accs, c='r', label='Test Accuracy')
    ax2.set_title('Accuracies (train-b, test-r)')
    ax2.set_ylabel('Accuracy')
    ax2.set_xlabel('Epoch')
    ax2.legend()
    ax2.axhline(y=0.5, color='gray', linestyle='--', linewidth=0.8)

    plt.savefig(output_folder+'/'+output_name, dpi=330 ,bbox_inches='tight',  pad_inches = 0.25)
