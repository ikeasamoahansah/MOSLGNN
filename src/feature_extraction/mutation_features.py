import pandas as pd
import numpy as np


def process_detailed_mutations(mutations_df):
    """
    Process detailed mutation data into binary matrix.

    Parameters:
    -----------
    mutations_df : DataFrame
        Detailed mutations with columns: ModelID, HugoSymbol, VariantInfo, etc.

    Returns:
    --------
    mutation_matrix : DataFrame
        Binary matrix (cell_lines x genes)
    """
    print("  Processing detailed mutation data...")

    # Filter for damaging mutations (optional)
    # You can filter by: isCOSMIChotspot, isDeleterious, etc.
    if 'isDeleterious' in mutations_df.columns:
        mutations_df = mutations_df[mutations_df['isDeleterious'] == True]
        print(f"  Filtered to {len(mutations_df)} damaging mutations")

    # Create binary matrix
    # Pivot: rows=ModelID (cell lines), columns=HugoSymbol (genes), values=1 (mutated)
    mutation_matrix = mutations_df.pivot_table(
        index='ModelID',
        columns='HugoSymbol',
        values='VariantInfo',  # Any column works since we just need presence
        aggfunc='count',
        fill_value=0
    )

    # Convert counts to binary (mutated=1, not mutated=0)
    mutation_matrix = (mutation_matrix > 0).astype(int)

    return mutation_matrix


mutation_freq_cache = None
mutation_mask_cache = None

def precompute_mutation_stats():
    """
    Precompute mutation statistics for all genes to speed up feature extraction.
    Call this once after loading data, before extracting features for pairs.
    """
    if cell_line_mutations is None:
        print("No mutation data loaded, skipping precomputation.")
        return

    print("Precomputing mutation statistics...")

    # Store mutation frequencies for all genes
    mutation_freq_cache = {}
    for gene in cell_line_mutations.columns:
        mut_values = cell_line_mutations[gene].values
        mutation_freq_cache[gene] = np.mean(mut_values)

    # Precompute masks for faster indexing
    # Store which cell lines have mutations for each gene
    mutation_mask_cache = {}
    for gene in cell_line_mutations.columns:
        mutation_mask_cache[gene] = cell_line_mutations[gene].values > 0

    print(f"Precomputed stats for {len(mutation_freq_cache)} genes")


def compute_mutation_context_features(gene_a, gene_b, cell_line_mutations):
    """
    Compute context-specific features based on mutations (OPTIMIZED).

    SL is often context-dependent: BRCA1-mutant cells depend on PARP1,
    but BRCA1-wildtype cells don't.

    This is KEY for understanding synthetic lethality mechanisms!

    Parameters:
    -----------
    gene_a, gene_b : str
        Gene symbols

    Returns:
    --------
    features : dict
    """
    features = {}

    if cell_line_mutations is None:
        return {}

    # Quick check if genes exist
    if (gene_a not in genesdf.columns or
        gene_b not in genesdf.columns):
        return {}

    # Use numpy arrays directly (faster than pandas indexing)
    effects_a = genesdf[gene_a].values
    effects_b = genesdf[gene_b].values

    # Precompute valid mask once
    valid_mask = ~(np.isnan(effects_a) | np.isnan(effects_b))
    effects_a_clean = effects_a[valid_mask]
    effects_b_clean = effects_b[valid_mask]

    # Get mutation data
    has_mut_a = gene_a in cell_line_mutations.columns
    has_mut_b = gene_b in cell_line_mutations.columns

    # Use cached masks if available, otherwise compute
    if mutation_mask_cache:
        mut_a_mask = mutation_mask_cache.get(gene_a, None) if has_mut_a else None
        mut_b_mask = mutation_mask_cache.get(gene_b, None) if has_mut_b else None
    else:
        mut_a_mask = (cell_line_mutations[gene_a].values > 0) if has_mut_a else None
        mut_b_mask = (cell_line_mutations[gene_b].values > 0) if has_mut_b else None

    # Apply valid mask to mutation masks
    if mut_a_mask is not None:
        mut_a_mask = mut_a_mask[valid_mask]
    if mut_b_mask is not None:
        mut_b_mask = mut_b_mask[valid_mask]

    # CONTEXT-DEPENDENT DEPENDENCY (vectorized operations)
    if mut_a_mask is not None:
        n_mutated_a = np.sum(mut_a_mask)
        if n_mutated_a > 5:  # Minimum threshold
            # Vectorized mean calculation
            effect_b_when_a_mutated = np.mean(effects_b_clean[mut_a_mask])
            effect_b_when_a_wt = np.mean(effects_b_clean[~mut_a_mask]) if n_mutated_a < len(mut_a_mask) else 0

            features['mutation_context_dependency_a_to_b'] = effect_b_when_a_wt - effect_b_when_a_mutated
            features['mutation_effect_b_in_mutant_a'] = effect_b_when_a_mutated

    if mut_b_mask is not None:
        n_mutated_b = np.sum(mut_b_mask)
        if n_mutated_b > 5:
            effect_a_when_b_mutated = np.mean(effects_a_clean[mut_b_mask])
            effect_a_when_b_wt = np.mean(effects_a_clean[~mut_b_mask]) if n_mutated_b < len(mut_b_mask) else 0

            features['mutation_context_dependency_b_to_a'] = effect_a_when_b_wt - effect_a_when_b_mutated
            features['mutation_effect_a_in_mutant_b'] = effect_a_when_b_mutated

    # MUTATION FREQUENCY (use cache if available)
    if mutation_freq_cache:
        if gene_a in mutation_freq_cache:
            features['mutation_frequency_a'] = mutation_freq_cache[gene_a]
        if gene_b in mutation_freq_cache:
            features['mutation_frequency_b'] = mutation_freq_cache[gene_b]
    else:
        if mut_a_mask is not None:
            features['mutation_frequency_a'] = np.mean(mut_a_mask)
        if mut_b_mask is not None:
            features['mutation_frequency_b'] = np.mean(mut_b_mask)

        # MUTUAL EXCLUSIVITY (vectorized)
        if mut_a_mask is not None and mut_b_mask is not None:
            both_mutated = np.sum(mut_a_mask & mut_b_mask)
            either_mutated = np.sum(mut_a_mask | mut_b_mask)

            features['mutation_co_occurrence_ratio'] = both_mutated / either_mutated if either_mutated > 0 else 0
            features['mutation_both_mutated_count'] = both_mutated
            features['mutation_either_mutated_count'] = either_mutated

        return features
