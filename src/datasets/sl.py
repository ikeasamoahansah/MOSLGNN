import pandas as pd
import numpy as np
from src import config

def load_sl_data(file_path):
    file_path = config.DATA_DIR / file_path
    df = pd.read_csv(file_path)

    # separate computated sl pairs from others
    real, comp = separate_sl_pairs(df)

    return real, comp


def load_non_sl_data(file_path):
    file_path = config.DATA_DIR / file_path
    df = pd.read_csv(file_path)

    return df


def separate_sl_pairs(df):
    comp_sl = pd.DataFrame(
        df[df['rel_source'] == "Computational Prediction"]
    )


    rows_drop = df[df['rel_surce'] == "Computational Prediction"].index

    df.drop(rows_drop, inplace=True)
    
    return df, comp_sl
