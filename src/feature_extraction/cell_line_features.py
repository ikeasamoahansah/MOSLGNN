import pandas as pd
import numpy as np
from scipy.stats import pearsonr, spearmanr


def compute_codependency_features(gene_a, gene_b, genesdf):
    """
    Compute co-dependency features between two genes.

    This is the KEY for SL prediction: genes that are synthetic lethal
    show correlated dependency patterns across cell lines.

    Parameters:
    -----------
    gene_a, gene_b : str
        Gene symbols
    genesdf: dataframe
        Gene effect data

    Returns:
    --------
    features : dict
    """
    features = {}

    # Check if genes exist in data
    if gene_a not in genesdf.columns or gene_b not in genesdf.columns:
        return empty_features()

    effects_a = genesdf[gene_a].values
    effects_b = genesdf[gene_b].values

    # Remove NaN values before computing correlations and other statistics
    mask = ~(np.isnan(effects_a) | np.isnan(effects_b))
    effects_a = effects_a[mask]
    effects_b = effects_b[mask]

    if len(effects_a) < 10:  # Not enough data after removing NaNs
        print("Not enough data after removing NaNs")
        return empty_features()

    # 1. PEARSON CORRELATION (linear relationship)
    # Positive correlation = both essential/non-essential together
    # Negative correlation = compensatory relationship
    pearson_corr, _ = pearsonr(effects_a, effects_b)
    features['depmap_pearson_correlation'] = pearson_corr

    # 2. SPEARMAN CORRELATION (monotonic relationship)
    spearman_corr, _ = spearmanr(effects_a, effects_b)
    features['depmap_spearman_correlation'] = spearman_corr

    # 3. CONDITIONAL DEPENDENCY
    # When gene A is essential (negative score), is gene B also essential?
    threshold_a = np.percentile(effects_a, 25)  # Bottom 25% = most essential
    essential_a_mask = effects_a < threshold_a

    if np.sum(essential_a_mask) > 0:
        # Average effect of gene B when gene A is essential
        conditional_effect_b = np.mean(effects_b[essential_a_mask])
        features['depmap_conditional_dependency'] = conditional_effect_b
    else:
        features['depmap_conditional_dependency'] = 0

    # 4. MUTUAL ESSENTIALITY
    # How often are both genes essential in the same cell lines?
    threshold_b = np.percentile(effects_b, 25)
    essential_b_mask = effects_b < threshold_b

    mutual_essentiality = np.sum(essential_a_mask & essential_b_mask) / len(effects_a)
    features['depmap_mutual_essentiality'] = mutual_essentiality

    # 5. DIFFERENTIAL ESSENTIALITY
    # Standard deviation of the difference
    diff = effects_a - effects_b
    features['depmap_essentiality_diff_std'] = np.std(diff)
    features['depmap_essentiality_diff_mean'] = np.mean(np.abs(diff))

    # 6. INDIVIDUAL GENE STATISTICS
    features['depmap_mean_effect_a'] = np.mean(effects_a)
    features['depmap_mean_effect_b'] = np.mean(effects_b)
    features['depmap_std_effect_a'] = np.std(effects_a)
    features['depmap_std_effect_b'] = np.std(effects_b)

    # 7. ESSENTIALITY SCORES
    # Negative mean = essential gene
    features['depmap_is_essential_a'] = 1 if np.mean(effects_a) < -0.5 else 0
    features['depmap_is_essential_b'] = 1 if np.mean(effects_b) < -0.5 else 0

    # 8. COMPLEMENTARY ESSENTIALITY (key for SL)
    # One gene essential, other not
    features['depmap_complementary'] = abs(
        features['depmap_is_essential_a'] - features['depmap_is_essential_b']
    )

    return features


def empty_depmap_features():
    """Return zero-filled features when data is missing."""
    return {
        'depmap_pearson_correlation': 0,
        'depmap_spearman_correlation': 0,
        'depmap_conditional_dependency': 0,
        'depmap_mutual_essentiality': 0,
        'depmap_essentiality_diff_std': 0,
        'depmap_essentiality_diff_mean': 0,
        'depmap_mean_effect_a': 0,
        'depmap_mean_effect_b': 0,
        'depmap_std_effect_a': 0,
        'depmap_std_effect_b': 0,
        'depmap_is_essential_a': 0,
        'depmap_is_essential_b': 0,
        'depmap_complementary': 0,
    }
