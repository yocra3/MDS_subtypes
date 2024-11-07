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
from collections import defaultdict
from pycox.models.loss import CoxPHLoss

sys.path.append('./scripts/3.GeneticScores/pytorch_files/')
from NNsparsemodel import sparseSurvivalNN_cluster

# Setup target device
device = "cuda" if torch.cuda.is_available() else "cpu"
## Remove warnings
import torch._dynamo
torch._dynamo.config.suppress_errors = True

## Load gene - GO mapping
gene_map =  pd.read_csv("results/mutations/go_gene_map_filtered.tsv", sep = "\t")
n_gos = len(gene_map["PathwayID"].unique())

## Prepare data
df_train = pd.read_csv("results/mutations/boolean_cluster_mutations_train.tsv", sep = "\t")
df_test = pd.read_csv("results/mutations/boolean_cluster_mutations_test.tsv", sep = "\t")

create_input = lambda df: (torch.tensor(df.iloc[:, 18:-4].to_numpy(), dtype=torch.float32), 
                           torch.tensor(df.iloc[:, 1:18].to_numpy(), dtype=torch.float32))

get_mutations = lambda df: torch.tensor(df.iloc[:, 18:-4].to_numpy(), dtype=torch.float32)


get_target = lambda df: (torch.tensor(df['OS_YEARS'].values, dtype=torch.float32), 
                          torch.tensor(df['OS_STATUS'].values, dtype=torch.float32))
get_target2 = lambda df: (df['OS_YEARS'].values, df['OS_STATUS'].values)



x_train = create_input(df_train)
x_test = create_input(df_test)

durations_test, events_test = get_target2(df_test)
y_train = (get_target(df_train), get_mutations(df_train))
y_test = (get_target(df_test), get_mutations(df_test))

## Create go-genes dictionary
genes = df_train.columns
genes = genes[18:-4]

gos = gene_map['PathwayID'].unique()
go_dict, gen_dict = dict(zip(gos, range(len(gos)))), dict(zip(genes, range(len(genes))))
rel_dict = defaultdict(list)
for go, gen in zip(gene_map['PathwayID'], gene_map['gene']):
    rel_dict[gen_dict[gen]].append(go_dict[go])



## Define loss
class LossAECoxPH(nn.Module):
    def __init__(self, alpha):
        super().__init__()
        assert (alpha >= 0) and (alpha <= 1), 'Need `alpha` in [0, 1].'
        self.alpha = alpha
        self.loss_surv = CoxPHLoss()
        self.loss_ae = nn.MSELoss()
        
    def forward(self, phi, decoded, target_loghaz, target_ae):
        idx_durations, events = target_loghaz
        loss_surv = self.loss_surv(phi, idx_durations, events)
        loss_ae = self.loss_ae(decoded, target_ae)
        return self.alpha * loss_surv + (1 - self.alpha) * loss_ae

loss = LossAECoxPH(0.6) ## Proportion of each loss

## Define architecture
n_genes = len(genes)
hidden_units = 10 ## Equals to number of gene sets
bol_net = sparseSurvivalNN_cluster(n_genes, n_gos, hidden_units, relation_dict = rel_dict, device = device)

summary(model = bol_net,
        input_size=((32, 115), (32, 17)), # make sure this is "input_size", not "input_shape"
        col_names=["input_size", "output_size", "num_params", "trainable"],
)
bol_model = CoxPH(bol_net, tt.optim.Adam, loss = loss)
batch_size = 32


## Train model
### Define learning rate
lrfinder = bol_model.lr_finder(x_train, y_train, batch_size, tolerance=10)
lr_plot = lrfinder.plot()
plt.savefig('figures/NN_models/lr_finder_sparseNN.png', dpi=300, format='png')

lrfinder.get_best_lr()
bol_model.optimizer.set_lr(0.0001)
epochs = 70

metrics = dict(
    loss_surv = LossAECoxPH(1),
    loss_ae   = LossAECoxPH(0)
)
bol_log = bol_model.fit(x_train, y_train, batch_size, epochs, verbose = True, 
                val_data = (x_test, y_test), val_batch_size = batch_size, 
                metrics = metrics)

res_bol = bol_model.log.to_pandas()
_ = res_bol[['train_loss', 'val_loss']].plot()
plt.savefig('figures/NN_models/boolean_cluster_sparseNN_training.png', dpi=300, format='png')

_ = res_bol[['train_loss_surv', 'val_loss_surv']].plot()
plt.savefig('figures/NN_models/boolean_cluster_sparseNN_survloss.png', dpi=300, format='png')

_ = res_bol[['train_loss_ae', 'val_loss_ae']].plot()
plt.savefig('figures/NN_models/boolean_cluster_sparseNN_aeloss.png', dpi=300, format='png')

## Evaluate model
_ = bol_model.compute_baseline_hazards()
surv = bol_model.predict_surv_df(x_test)
pd.Series.is_monotonic = pd.Series.is_monotonic_increasing
ev = EvalSurv(surv, durations_test, events_test, censor_surv='km')
ev.concordance_td()

risk = bol_model.predict(x_test)

concordance_index(durations_test, -risk.cpu(), events_test)
# 0.689545085784682

df_test_comp = df_test.dropna(subset=["IPSSM_SCORE"])
durations_test_comp, events_test_comp = get_target2(df_test_comp)

concordance_index(durations_test_comp, -df_test_comp['IPSSM_SCORE'], events_test_comp)
# 0.739408163967124

bol_model.save_net("results/mutations/sparseNN_cluster_boolean_model.pt")
bol_model.save_model_weights("results/mutations/sparseNN_cluster_boolean_weights.pt")

## Use vaf data
df_vaf_train = pd.read_csv("results/mutations/vaf_cluster_mutations_train.tsv", sep = "\t")
df_vaf_test = pd.read_csv("results/mutations/vaf_cluster_mutations_test.tsv", sep = "\t")


x_vaf_train = create_input(df_vaf_train)
x_vaf_test = create_input(df_vaf_test)

durations_vaf_test, events_vaf_test = get_target2(df_vaf_test)

y_vaf_train = (get_target(df_vaf_train), get_mutations(df_vaf_train))
y_vaf_test = (get_target(df_vaf_test), get_mutations(df_vaf_test))

## Define model
vaf_net = sparseSurvivalNN_cluster(n_genes, n_gos, hidden_units, relation_dict = rel_dict, device = device)
optimizer = torch.optim.Adam(params=vaf_net.parameters(), lr = 1e-4)

model_vaf = CoxPH(vaf_net, tt.optim.Adam, loss = loss)
batch_size = 32

## Train model
epochs = 20

log_vaf = model_vaf.fit(x_vaf_train, y_vaf_train, batch_size, epochs, verbose = True, 
                val_data = (x_vaf_test, y_vaf_test), val_batch_size = batch_size, 
                metrics = metrics)


res_vaf = model_vaf.log.to_pandas()
_ = res_vaf[['train_loss', 'val_loss']].plot()
plt.savefig('figures/NN_models/vaf_cluster_sparseNN_training.png', dpi=300, format='png')

_ = res_vaf[['train_loss_surv', 'val_loss_surv']].plot()
plt.savefig('figures/NN_models/vaf_cluster_sparseNN_survloss.png', dpi=300, format='png')

_ = res_vaf[['train_loss_ae', 'val_loss_ae']].plot()
plt.savefig('figures/NN_models/vaf_cluster_sparseNN_aeloss.png', dpi=300, format='png')


_ = model_vaf.compute_baseline_hazards()
surv_vaf = model_vaf.predict_surv_df(x_vaf_test)
ev_vaf = EvalSurv(surv_vaf, durations_vaf_test, events_vaf_test, censor_surv='km')
ev_vaf.concordance_td()

risk_vaf = model_vaf.predict(x_vaf_test)

concordance_index(durations_vaf_test, -risk_vaf.cpu(), events_vaf_test)
# 0.6919638832967759

model_vaf.save_net("results/mutations/sparseNN_cluster_vaf_model.pt")
model_vaf.save_model_weights("results/mutations/sparseNN_cluster_vaf_weights.pt")
