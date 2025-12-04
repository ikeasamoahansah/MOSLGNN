import pandas as pd
import numpy as np

from cell_line_features import compute_codependency_features
from mutation_features import compute_mutation_context_features
from ppi_features import compute_string_features 
from pathway_features import compute_kegg_features


def extract_features_for_pair(gene_a, gene_b, genesdf, cell_line_mutations, string_data, kegg_pathways):
    """
    Extract all features for a gene pair.

    Parameters:
    -----------
    gene_a, gene_b : str
        Gene symbols

    Returns:
    --------
    features : dict
    """
    features = {}
    features.update(compute_codependency_features(gene_a, gene_b, genesdf))
    features.update(compute_mutation_context_features(gene_a, gene_b, cell_line_mutations))
    features.update(compute_string_features(string_data, gene_a, gene_b))
    features.update(compute_kegg_features(gene_a, gene_b, kegg_pathways))

    if 'string_combined_score' in features and 'depmap_pearson_correlation' in features:
        features['combined_string_depmap_interaction'] = (
            features['string_combined_score'] * abs(features['depmap_pearson_correlation'])
        )

        # Physical interaction with complementary essentiality
        if 'depmap_complementary' in features:
            features['combined_physical_complementary'] = (
                features.get('string_physical_interaction', 0) *
                features['depmap_complementary']
            )
    return features
