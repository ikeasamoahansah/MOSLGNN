import pandas as pd
import numpy as np
from src import config


def load_ppi(file_path):
    file_path = config.DATA_DIR / file_path
    # add separation for txt file
    df = pd.read_csv(file_path, sep=' ')
    
    return df


def load_pi(file_path):
    file_path = config.DATA_DIR / file_path
    # add separation for txt file
    df = pd.read_csv(file_path, sep=' ')

    protein_info = {}

    for protein, name in zip(df['#string_protein_id'], df['preferred_name']):
        protein_info[protein] = name

    return df, protein_info


