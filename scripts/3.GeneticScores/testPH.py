"""
Module Summary:
Evaluate the proportional hazards assumption using Cox Proportional Hazards model.
This script runs a CoxPH model on the provided data and checks the proportional hazards assumption.
docker run --gpus all -v $PWD:$PWD -w $PWD -it mds_subtypes_python:1.5 python

"""

import os
import sys
import random
import numpy as np
import pandas as pd

from lifelines import CoxPHFitter
from lifelines.statistics import proportional_hazard_test

patient_vals = pd.read_csv("results/gnn/preprocess/patient_variables.tsv", sep = "\t")
patient_vals['PLT'] = np.log(patient_vals['PLT'] + 1) 
all_vars = patient_vals.drop(['OS_YEARS', 'OS_STATUS', 'AMLt_YEARS', 'AMLt_STATUS', 'ID', 'IPSSM_SCORE', 'train'], axis=1).copy()


cph = CoxPHFitter()
cph.fit(all_vars, duration_col='LFS_YEARS', event_col='LFS_STATUS')

results = proportional_hazard_test(cph, all_vars, time_transform='rank')
print(results.summary)

## Se cumple el PH
