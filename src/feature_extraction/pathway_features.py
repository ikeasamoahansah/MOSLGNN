import pandas as pd
import numpy as np


def compute_kegg_features(gene_a, gene_b, kegg_pathways):
    """
    Extract KEGG pathway features for a gene pair.

    Parameters:
    -----------
    gene_a, gene_b : str
        Gene symbols

    Returns:
    --------
    features : dict
    """
    features = {}

    if kegg_pathways is None:
        return {}

    # Get pathways for each gene
    pathways_a = set(kegg_pathways.get(gene_a, []))
    pathways_b = set(kegg_pathways.get(gene_b, []))

    if len(pathways_a) == 0 and len(pathways_b) == 0:
        return _empty_kegg_features()

    # 1. SHARED PATHWAYS
    # SL genes often function in parallel pathways (e.g., DNA repair)
    shared_pathways = pathways_a & pathways_b
    features['kegg_shared_pathways'] = len(shared_pathways)
    features['kegg_in_same_pathway'] = 1 if len(shared_pathways) > 0 else 0

    # 2. JACCARD SIMILARITY
    # Measure pathway overlap
    union_pathways = pathways_a | pathways_b
    if len(union_pathways) > 0:
        features['kegg_pathway_jaccard'] = len(shared_pathways) / len(union_pathways)
    else:
        features['kegg_pathway_jaccard'] = 0

    # 3. PATHWAY COUNTS
    features['kegg_pathways_a'] = len(pathways_a)
    features['kegg_pathways_b'] = len(pathways_b)
    features['kegg_total_pathways'] = len(union_pathways)

    # 4. PATHWAY CATEGORY OVERLAP
    # KEGG pathways have hierarchical structure: hsa##### where first 2 digits = category
    # Categories: 00=Metabolism, 01=Genetic Info, 02=Environmental Info,
    #             03=Cellular Processes, 04=Organismal Systems, 05=Diseases
    categories_a = set([p[:5] for p in pathways_a if len(p) >= 5])
    categories_b = set([p[:5] for p in pathways_b if len(p) >= 5])

    shared_categories = categories_a & categories_b
    features['kegg_shared_categories'] = len(shared_categories)

    # 5. SPECIFIC PATHWAY TYPES (relevant for SL)
    dna_repair_pathways = {'hsa03410', 'hsa03420', 'hsa03430', 'hsa03440', 'hsa03450', 'hsa03460'}
    cell_cycle_pathways = {'hsa04110', 'hsa04111', 'hsa04112', 'hsa04113', 'hsa04114', 'hsa04115'}

    features['kegg_both_in_dna_repair'] = 1 if (pathways_a & dna_repair_pathways) and (pathways_b & dna_repair_pathways) else 0
    features['kegg_both_in_cell_cycle'] = 1 if (pathways_a & cell_cycle_pathways) and (pathways_b & cell_cycle_pathways) else 0

    # 6. COMPLEMENTARY PATHWAYS
    # Measure if genes are in different but related pathways
    features['kegg_complementary_pathways'] = len(pathways_a) + len(pathways_b) - 2 * len(shared_pathways)

    return features


def empty_kegg_features():
    """Return zero-filled KEGG features when data is missing."""
    return {
        'kegg_shared_pathways': 0,
        'kegg_in_same_pathway': 0,
        'kegg_pathway_jaccard': 0,
        'kegg_pathways_a': 0,
        'kegg_pathways_b': 0,
        'kegg_total_pathways': 0,
        'kegg_shared_categories': 0,
        'kegg_both_in_dna_repair': 0,
        'kegg_both_in_cell_cycle': 0,
        'kegg_complementary_pathways': 0,
    }
