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


def load_string_data(string_df, score_threshold=550):
    """
    Load STRING protein-protein interaction data.

    STRING database format:
    protein1 | protein2 | combined_score | experimental | database |
    coexpression | neighborhood | fusion | cooccurrence | textmining

    Download from: https://string-db.org/cgi/download
    File: 9606.protein.links.detailed.v12.0.txt (for human)

    Parameters:
    -----------
    string_path : str
        Path to STRING interactions file
    score_threshold : int
        Minimum combined score (0-1000). Default 0 includes all.
        Recommended: 400 (medium confidence), 700 (high confidence)
    """

    # Filter by score threshold
    if score_threshold > 0:
        string_df = string_df[string_df['combined_score'] >= score_threshold]
        print(f"Filtered to {len(string_df)} interactions with score >= {score_threshold}")

    print(f"Loaded {len(string_df)} protein interactions")

    # Store as dictionary for fast lookup
    # Key: (gene_a, gene_b), Value: interaction scores
    string_data = {}

    for _, row in string_df.iterrows():
        # Extract gene symbols from protein IDs
        # Format: "9606.ENSP00000000233" -> need mapping to gene symbol
        protein1 = row['protein1']
        protein2 = row['protein2']

        # For now, use protein IDs directly
        # In practice, you'd map ENSP IDs to gene symbols
        pair = tuple(sorted([protein1, protein2]))

        string_data[pair] = {
            'combined_score': row['combined_score'],
            'experimental': row.get('experimental', 0),
            'database': row.get('database', 0),
            'coexpression': row.get('coexpression', 0),
            'neighborhood': row.get('neighborhood', 0),
            'fusion': row.get('fusion', 0),
            'cooccurrence': row.get('cooccurrence', 0),
            'textmining': row.get('textmining', 0),
        }

    print(f"STRING data indexed: {len(string_data)} unique protein pairs")
    return string_data
