import pandas as pd
import numpy as np


def generate_negative_pairs(n_pairs, known_sl_pairs, genesdf):
    """Generate random gene pairs as negative examples."""
    genes = list(genesdf.columns)
    negative_pairs = []

    # Normalize known SL pairs (handle both orderings)
    known_sl_set = set()
    for gene_a, gene_b in known_sl_pairs:
        # Add both orderings since (A,B) == (B,A)
        known_sl_set.add(tuple(sorted([gene_a, gene_b])))

    print(f"Excluding {len(known_sl_set)} known SL pairs from negatives...")

    max_attempts = n_pairs * 10  # Prevent infinite loop
    attempts = 0

    while len(negative_pairs) < n_pairs and attempts < max_attempts:
        gene_a, gene_b = np.random.choice(genes, size=2, replace=False)
        candidate_pair = tuple(sorted([gene_a, gene_b]))

        # Check if this pair is in known SL pairs
        if candidate_pair not in known_sl_set:
            negative_pairs.append((gene_a, gene_b))
        else:
            print(f"  Skipped known SL pair: {candidate_pair}")

        attempts += 1

    if len(negative_pairs) < n_pairs:
        print(f"Warning: Could only generate {len(negative_pairs)} negative pairs")

    return negative_pairs


def validate_negative_pairs(negative_pairs, known_sl_pairs):
    """
    Validate that negative pairs don't overlap with known SL pairs.

    Parameters:AAAS
    -----------
    negative_pairs : list of tuples
        Candidate negative pairs
    known_sl_pairs : list of tuples
        Known SL pairs to exclude

    Returns:
    --------
    validated_pairs : list of tuples
        Negative pairs with no overlap
    """
    known_sl_set = set()
    for gene_a, gene_b in known_sl_pairs:
        known_sl_set.add(tuple(sorted([gene_a, gene_b])))

    validated_pairs = []
    removed_count = 0

    for gene_a, gene_b in negative_pairs:
        candidate_pair = tuple(sorted([gene_a, gene_b]))
        if candidate_pair not in known_sl_set:
            validated_pairs.append((gene_a, gene_b))
        else:
            removed_count += 1
            print(f"  Removed known SL pair from negatives: {candidate_pair}")

    if removed_count > 0:
        print(f"  Total removed: {removed_count} pairs")

    return validated_pairs


