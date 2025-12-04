import pandas as pd
import numpy as np


def compute_string_features(string_data, gene_a, gene_b):
    """
    Extract protein-protein interaction features from STRING database.

    STRING provides multiple evidence channels for interactions:
    - experimental: lab-verified interactions
    - database: curated interaction databases
    - coexpression: genes with correlated expression patterns
    - neighborhood: genes close on chromosome
    - fusion: gene fusion evidence
    - cooccurrence: phylogenetic co-occurrence
    - textmining: text mining from scientific literature

    Parameters:
    -----------
    gene_a, gene_b : str
        Gene symbols

    Returns:
    --------
    features : dict
    """
    features = {}

    if string_data is None:
        print("Warning: STRING data not loaded. Use load_string_data() first.")
        return _empty_string_features()

    # Check both orderings since interactions are bidirectional
    pair = tuple(sorted([gene_a, gene_b]))

    if pair not in string_data:
        return _empty_string_features()

    interaction = string_data[pair]

    # Normalize scores to 0-1 range (STRING scores are 0-1000)
    features['string_combined_score'] = interaction['combined_score'] / 1000.0
    features['string_experimental'] = interaction.get('experimental', 0) / 1000.0
    features['string_database'] = interaction.get('database', 0) / 1000.0
    features['string_coexpression'] = interaction.get('coexpression', 0) / 1000.0
    features['string_neighborhood'] = interaction.get('neighborhood', 0) / 1000.0
    features['string_fusion'] = interaction.get('fusion', 0) / 1000.0
    features['string_cooccurrence'] = interaction.get('cooccurrence', 0) / 1000.0
    features['string_textmining'] = interaction.get('textmining', 0) / 1000.0

    # Binary flags for different confidence levels
    features['string_has_interaction'] = 1 if interaction['combined_score'] > 0 else 0
    features['string_medium_confidence'] = 1 if interaction['combined_score'] >= 400 else 0
    features['string_high_confidence'] = 1 if interaction['combined_score'] >= 700 else 0

    # Evidence diversity (how many channels provide evidence)
    evidence_channels = [
        interaction.get('experimental', 0),
        interaction.get('database', 0),
        interaction.get('coexpression', 0),
        interaction.get('neighborhood', 0),
        interaction.get('fusion', 0),
        interaction.get('cooccurrence', 0),
    ]
    features['string_evidence_count'] = sum(1 for score in evidence_channels if score > 0)

    # Physical vs functional interaction
    # Physical: experimental + database
    # Functional: coexpression + cooccurrence
    physical_score = (interaction.get('experimental', 0) +
                        interaction.get('database', 0)) / 2000.0
    functional_score = (interaction.get('coexpression', 0) +
                        interaction.get('cooccurrence', 0)) / 2000.0

    features['string_physical_interaction'] = physical_score
    features['string_functional_association'] = functional_score
    features['string_physical_vs_functional'] = physical_score - functional_score

    return features


def empty_features():
    """Return zero-filled STRING features when data is missing."""
    return {
        'string_combined_score': 0,
        'string_experimental': 0,
        'string_database': 0,
        'string_coexpression': 0,
        'string_neighborhood': 0,
        'string_fusion': 0,
        'string_cooccurrence': 0,
        'string_textmining': 0,
        'string_has_interaction': 0,
        'string_medium_confidence': 0,
        'string_high_confidence': 0,
        'string_evidence_count': 0,
        'string_physical_interaction': 0,
        'string_functional_association': 0,
        'string_physical_vs_functional': 0,
    }
