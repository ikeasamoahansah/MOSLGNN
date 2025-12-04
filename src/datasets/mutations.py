import pandas as pd
import numpy as np
from src import config


def load_mut(file_path):
    file_path = config.DATA_DIR / file_path
    df = pd.read_csv(file_path)

    # select relevant columns
    df = df[['Hugo_Symbol', 'Entrez_Gene_Id', 'Variant_Type']]

    return df
