# Docker inheritance
FROM  graphcore/pytorch-geometric:3.3.0

# Install python packages
RUN pip3 install pandas lifelines

# Install pytorch version compatible with cuda 
RUN pip3 install --upgrade torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118
RUN pip3 install torchinfo tensorboard


# Install pycox to run survival analysis (1.1)
RUN pip3 install pycox


