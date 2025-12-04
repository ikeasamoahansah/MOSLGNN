import pandas as pd
import numpy as np
import re

from src import config


def load_cell_info(file_path):
    file_path = config.DATA_DIR / file_path
    df = pd.read_csv(file_path)
    
    # remove trailing (number) from column names
    df = remove_trailing_number(df)

    # rename column
    df = rename_unnamed(df)

    # convert all remnaining columns to numeric
    for col in df.columns:
        df[col] = pd.to_numeric(df[col], errors='coerce')
   
    print("Processed cell line data from depmap......100% complete")
    return df


def remove_trailing_number(df):
    """
    Remove trailing " (number)" from column names, e.g. "A1BG (1)" -> "A1BG"
    """
    return df.rename(
        columns=lambda c: re.sub(r"\s*\(\d+\)$", "", c) if isinstance(c, str) else c
    )


def rename_unnamed(df):
    """
    Rename Unnamed:0 column to ModelID and set index
    """
    if 'Unnamed: 0' in df.columns:
        df = df.set_index('Unnamed: 0')
        df.index.name = 'ModelID'

    return df

