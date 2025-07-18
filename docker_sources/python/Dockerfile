# Docker inheritance
FROM  graphcore/pytorch-geometric:3.3.0

# Install python packages
RUN pip3 install pandas lifelines

# Install pytorch version compatible with cuda 
RUN pip3 install --upgrade torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118
RUN pip3 install torchinfo tensorboard


# Install pycox to run survival analysis (1.1)
RUN pip3 install pycox

# Install pycox to run survival analysis (1.2)
RUN pip3 install torch_geometric

# Install scgpt to get gene embeddings (1.3)
RUN pip3 install scgpt

# Install lightning to speed implementation (1.4)
RUN pip3 install lightning

# Install optuna to define hyperparameters (1.5)
RUN pip3 install optuna optuna-integration[pytorch_lightning] plotly

# Install pytest to test Claude
RUN pip3 install pytest
