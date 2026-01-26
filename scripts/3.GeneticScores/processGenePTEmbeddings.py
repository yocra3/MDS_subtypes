"""
Get gene embeddings using the pre-trained model from GenePT
docker run  -v $PWD:$PWD -w $PWD -it mds_subtypes_python:1.4 python

"""

# Load libraries
import pickle
import pandas as pd

# Load data
genept = pickle.load(open("data/GenePT_embedding_v2/GenePT_gene_protein_embedding_model_3_text.pickle", "rb"))


genept_df = pd.DataFrame.from_dict(genept, orient='index')
genept_df.to_csv('results/preprocess/gene_embedding_genept.tsv', sep='\t', index=True)

