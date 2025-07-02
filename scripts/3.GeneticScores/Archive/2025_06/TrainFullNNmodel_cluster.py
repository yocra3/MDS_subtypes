"""
Train model
docker run --gpus all -v $PWD:$PWD -w $PWD -it mds_subtypes_python:1.1 python

"""

import torch
from torch import nn
import pandas as pd
import sys
import numpy as np
from torchinfo import summary
from pycox.models import CoxPH
from pycox.evaluation import EvalSurv
import matplotlib.pyplot as plt
import torchtuples as tt
from lifelines.utils import concordance_index

## Import custom modules
# sys.path.append('./scripts/utils')
# from dataloader_tsv import data_loader_tsv
# from engine import train
# from tensorboard_writer import create_writer

sys.path.append('./scripts/3.GeneticScores')
from NNmodel_cluster import survivalNN_cluster

# Setup target device
device = "cuda" if torch.cuda.is_available() else "cpu"
## Remove warnings
import torch._dynamo
torch._dynamo.config.suppress_errors = True


## Prepare data
df_train = pd.read_csv("results/mutations/boolean_cluster_mutations_train.tsv", sep = "\t")
df_test = pd.read_csv("results/mutations/boolean_cluster_mutations_test.tsv", sep = "\t")

create_input = lambda df: (torch.tensor(df.iloc[:, 18:-4].to_numpy(), dtype=torch.float32), 
                           torch.tensor(df.iloc[:, 1:18].to_numpy(), dtype=torch.float32))
get_target = lambda df: (torch.tensor(df['OS_YEARS'].values, dtype=torch.float32), 
                          torch.tensor(df['OS_STATUS'].values, dtype=torch.float32))
get_target2 = lambda df: (df['OS_YEARS'].values, df['OS_STATUS'].values)

x_train = create_input(df_train)
x_test = create_input(df_test)

durations_test, events_test = get_target2(df_test)
y_train = get_target(df_train)
y_test = get_target(df_test)

## Define architecture
n_feats = x_train[0].shape[1]
hidden_units = (59, 10) ## Equals to number of gene sets
net = survivalNN_cluster(n_feats, hidden_units).to(device)

summary(model = net,
        input_size=((32, 115), (32, 17)), # make sure this is "input_size", not "input_shape"
        col_names=["input_size", "output_size", "num_params", "trainable"],
)

model = CoxPH(net, tt.optim.Adam)
batch_size = 32

## Train model
### Define learning rate
lrfinder = model.lr_finder(x_train, y_train, batch_size, tolerance=10)
lr_plot = lrfinder.plot()
plt.savefig('figures/NN_models/lr_finder_fullNN.png', dpi=300, format='png')

lrfinder.get_best_lr()
model.optimizer.set_lr(0.0001)
epochs = 30

log = model.fit(x_train, y_train, batch_size, epochs, verbose = True, 
                val_data = (x_test, y_test), val_batch_size = batch_size)


_ = log.plot()
plt.savefig('figures/NN_models/boolean_cluster_NNfull_training.png', dpi=300, format='png')

## Evaluate model
_ = model.compute_baseline_hazards()
surv = model.predict_surv_df(x_test)
pd.Series.is_monotonic = pd.Series.is_monotonic_increasing
ev = EvalSurv(surv, durations_test, events_test, censor_surv='km')
ev.concordance_td()

risk = model.predict(x_test)

concordance_index(durations_test, -risk.cpu(), events_test)
# 0.6911225624230041

df_test_comp = df_test.dropna(subset=["IPSSM_SCORE"])
durations_test_comp, events_test_comp = get_target2(df_test_comp)

concordance_index(durations_test_comp, -df_test_comp['IPSSM_SCORE'], events_test_comp)
# 0.739408163967124

model.save_net("results/mutations/fullNN_boolean_cluster_model.pt")
model.save_model_weights("results/mutations/fullNN_boolean_cluster_weights.pt")


## Use vaf data
df_vaf_train = pd.read_csv("results/mutations/vaf_cluster_mutations_train.tsv", sep = "\t")
df_vaf_test = pd.read_csv("results/mutations/vaf_cluster_mutations_test.tsv", sep = "\t")


x_vaf_train = create_input(df_vaf_train)
x_vaf_test = create_input(df_vaf_test)

durations_vaf_test, events_vaf_test = get_target2(df_vaf_test)
y_vaf_train = get_target(df_vaf_train)
y_vaf_test = get_target(df_vaf_test)

## Define model
net2 = survivalNN_cluster(n_feats, hidden_units).to(device)
optimizer = torch.optim.Adam(params=net2.parameters(), lr = 1e-4)

model_vaf = CoxPH(net2, optimizer)
batch_size = 32

## Train model
epochs = 60

log_vaf = model_vaf.fit(x_vaf_train, y_vaf_train, batch_size, epochs, verbose = True, 
                val_data = (x_vaf_test, y_vaf_test), val_batch_size = batch_size)


_ = log_vaf.plot()
plt.savefig('figures/NN_models/vaf_cluster_NNfull_training.png', dpi=300, format='png')


_ = model_vaf.compute_baseline_hazards()
surv_vaf = model_vaf.predict_surv_df(x_vaf_test)
ev_vaf = EvalSurv(surv_vaf, durations_vaf_test, events_vaf_test, censor_surv='km')
ev_vaf.concordance_td()

risk_vaf = model_vaf.predict(x_vaf_test)

concordance_index(durations_vaf_test, -risk_vaf.cpu(), events_vaf_test)
# 0.698685188103654

df_vaf_test_comp = df_vaf_test.dropna(subset=["IPSSM_SCORE"])
x_vaf_test_comp = create_input(df_vaf_test_comp)

risk_vaf_comp = model_vaf.predict(x_vaf_test_comp)
durations_vaf_test_comp, events_vaf_test_comp = get_target2(df_vaf_test_comp)

concordance_index(durations_vaf_test_comp, -risk_vaf_comp.cpu(), events_vaf_test_comp)
# 0.6988032600845971
concordance_index(durations_vaf_test_comp, -df_test_comp['IPSSM_SCORE'], events_vaf_test_comp)
# 0.739408163967124

model_vaf.save_net("results/mutations/fullNN_vaf_model.pt")
model_vaf.save_model_weights("results/mutations/fullNN_vaf_weights.pt")







# ## Define loss function
# def cox_ph_loss(predictions, durations, events):
#     risk = predictions.squeeze()
#     log_risk = torch.logcumsumexp(risk, dim=0)
#     uncensored_likelihood = risk - log_risk
#     censored_likelihood = uncensored_likelihood * events
#     return -torch.mean(censored_likelihood)

hidden_units = 103
## Train models
for experiment, exp_list in experiments.items():
    print(f"Training {experiment} model...")
    train_data, test_data, n_feats = exp_list
    ## Define model
    model = survivalNN(n_feats, hidden_units).to(device)
    model = torch.compile(model)
    ## Define function loss and optimizer
    optimizer = torch.optim.Adam(params=model.parameters(), lr = 1e-4)
    # Start training with help from engine.py
    train_res = train(model = model,
             train_dataloader=train_data,
             test_dataloader=test_data,
             loss_fn=cox_ph_loss,
             optimizer=optimizer,
             epochs = epochs,
             writer = "",
             device=device)

## Recover best model (full dataset 1%)
train_data, test_data, n_feats = experiments["vaf"]


test_pd =  pd.read_csv( "results/mutations/vaf_mutations_filtered_test.tsv", sep = "\t")
test_clin =  torch.tensor(test_pd.iloc[:, 0].values,dtype=torch.float32)
test_clin = torch.unsqueeze(test_clin, dim = -1)

test_muts =  torch.tensor(test_pd.iloc[:, 1:-2].values,dtype=torch.float32)
## Define model
model0 = survivalNN(n_feats, hidden_units).to(device)
model0 = torch.compile(model0)
epochs = 200

scores_in = model0(test_muts.to(device), test_clin.to(device))
scores_in = scores_in.cpu().detach().numpy()

## Define function loss and optimizer
optimizer = torch.optim.Adam(params=model0.parameters(), lr = 1e-4)
# Start training with help from engine.py
train_res = train(model = model0,
             train_dataloader=train_data,
             test_dataloader=test_data,
             loss_fn=cox_ph_loss,
             optimizer=optimizer,
             epochs = epochs,
             writer = None,
             device=device)

scores_end = model0(test_muts.to(device), test_clin.to(device))
scores_end = scores_end.cpu().detach().numpy()

concordance_index(test_pd.iloc[:, -2], test_pd.iloc[:, 0], test_pd.iloc[:, -1])
concordance_index(test_pd.iloc[:, -2], scores_in, test_pd.iloc[:, -1])


test_pd['Pred_probs'] = pred_probs.cpu().detach().numpy()
test_pd.to_csv("results/NNmodels/collapse_feats_1_percen_test_prediction.tsv", sep = "\t",
               index=False)