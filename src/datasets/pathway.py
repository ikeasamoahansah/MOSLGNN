import pandas as pd
import numpy as np
from src import config
import re


def load_kegg_from_files(gene_list_file, gene_pathway_link_file):
    """
    Load KEGG pathway data from downloaded KEGG files.

    1. Gene list file (hsa_gene.list):
        Format: hsa:38\tCDS 11:108116705..108147603 ACAT1, ACAT, MAT, T2, THIL; acetyl-CoA acetyltransferase 1

    2. Gene-pathway link file (hsa_gene_pathway.list):
        Format: hsa:38\tpath:hsa00071

    Parameters:
    -----------
    gene_list_file : str
        name of file
    gene_pathway_link_file : str
        name of file
    """
    print("Loading KEGG data from files...")

    # Step 1: Load gene ID to gene symbol mapping
    print("  1. Loading gene list and extracting gene symbols...")
    gene_list_file = config.DATA_DIR / gene_list_file
    gene_pathway_link_file = config.DATA_DIR / gene_pathway_link_file

    gene_id_to_symbol = {}

    with open(gene_list_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                gene_id = parts[0]  # e.g., "hsa:38"
                description = parts[1]  # e.g., "CDS 11:... ACAT1, ACAT, MAT; description"

                # Extract gene symbol (first name before semicolon or comma)
                # Format examples:
                # "CDS 11:... ACAT1, ACAT, MAT, T2, THIL; acetyl-CoA acetyltransferase 1"
                # "miRNA 11:... MIR4491, mir-4491; microRNA 4491"

                # Find the part after the coordinates and before semicolon
                if ';' in description:
                    gene_part = description.split(';')[0]
                else:
                    gene_part = description

                # Extract gene symbols (they come after the location info)
                # Look for patterns like "ACAT1, ACAT, MAT"
                # Match gene symbols (uppercase letters/numbers)
                symbol_match = re.search(r'[A-Z][A-Z0-9\-]+(?:,\s*[A-Z][A-Z0-9\-]+)*', gene_part)

                if symbol_match:
                    # Get first symbol (preferred name)
                    symbols = symbol_match.group().split(',')
                    gene_symbol = symbols[0].strip()
                    gene_id_to_symbol[gene_id] = gene_symbol

    print(f"Found {len(gene_id_to_symbol)} genes with symbols")

    # Step 2: Load gene-pathway links
    print("2. Loading gene-pathway links...")
    gene_to_pathways = {}

    with open(gene_pathway_link_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                gene_id = parts[0]  # e.g., "hsa:38"
                pathway_id = parts[1].replace('path:', '')  # e.g., "hsa00071"

                if gene_id not in gene_to_pathways:
                    gene_to_pathways[gene_id] = []
                gene_to_pathways[gene_id].append(pathway_id)

    print(f"Found {len(gene_to_pathways)} genes with pathway annotations")

    # Step 3: Create final mapping: gene_symbol -> pathways
    print("3. Creating gene symbol to pathway mapping...")
    kegg_pathways = {}

    for gene_id, pathway_list in gene_to_pathways.items():
        gene_symbol = gene_id_to_symbol.get(gene_id, None)
        if gene_symbol:
            kegg_pathways[gene_symbol] = pathway_list

    print(f"✓ KEGG data loaded: {len(kegg_pathways)} genes with pathway annotations")
    print(f"Example: {list(kegg_pathways.items())[:3]}")
    return kegg_pathways, gene_to_pathways, gene_id_to_symbol


def load_mapping(file_path):
    gene_id_to_gene = {}

    file_path = config.DATA_DIR / file_path
    # add separataion for tabbed txt file
    df = pd.read_csv(file_path, sep='\t')

    for _, row in df.iterrows():
        gene_id_to_gene[row['gene_id']] = row['preferred_name']

    print(f"Loaded {len(gene_id_to_gene)} gene mappings")

    return gene_id_to_gene
