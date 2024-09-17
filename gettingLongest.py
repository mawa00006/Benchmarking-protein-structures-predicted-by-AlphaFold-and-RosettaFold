"""
Script for extracting the longest sequence out of .cif files in the specified folder and saving them in a new file
"""

import os
from Bio.PDB.MMCIFParser import MMCIFParser

# Dictionary for converting 3-letter amino acid codes to 1-letter codes
three_to_one = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y',
    'SEC': 'A', 'PYL': 'A',  # Uncommon amino acids substituted by alanine for alphafold
    'ASX': 'A', 'GLX': 'A', 'XAA': 'A', 'UNK': 'A'  # Ambiguous/unknown substituted by alanine for alphafold
}


def convert_3letter_to_1letter(residue_name):
    """Convert 3-letter residue name to 1-letter amino acid code."""
    return three_to_one.get(residue_name, '')  # blank for unknown residues


def get_longest_sequence_in_cif(cif_file_path):
    """Extract the longest sequence from the given mmCIF file."""
    parser = MMCIFParser(QUIET=True)

    # Parse the mmCIF file
    structure = parser.get_structure('', cif_file_path)

    longest_sequence = ''

    # Loop through all the chains in the structure
    for model in structure:
        for chain in model:
            # Extract the sequence of the current chain
            print(''.join([res.resname for res in chain.get_residues()]))
            chain_sequence = ''.join([convert_3letter_to_1letter(res.resname) for res in chain.get_residues()])

            # Find the longest sequence
            if len(chain_sequence) > len(longest_sequence):
                longest_sequence = chain_sequence

    print(longest_sequence)
    return longest_sequence


def process_cif_files(folder_path, output_file):
    """Process all mmCIF files in the folder to get the longest sequence for each file and save them to a text file."""
    with open(output_file, 'w') as f_out:
        # Loop through all files in the folder
        for file_name in os.listdir(folder_path):
            if file_name.endswith(".cif"):
                file_path = os.path.join(folder_path, file_name)
                print(f"Processing {file_name}...")

                # Get the longest sequence from the mmCIF file
                longest_sequence = get_longest_sequence_in_cif(file_path)

                # Write the results to the output file
                f_out.write(f">Longest sequence in {file_name}\n")
                f_out.write(f"Sequence: {longest_sequence}\n")
                f_out.write(f"Length: {len(longest_sequence)}\n")
                f_out.write("-" * 40 + "\n")
                print(f"Saved sequence from {file_name} to output file.")


# Specify the folder containing the mmCIF files and the output file path
cif_folder = 'mmcif_files'  # Change this to your folder path
output_file = 'longest_sequences.txt'  # The output text file

# Step 1: Process the mmCIF files and save the longest sequence for each to a text file
process_cif_files(cif_folder, output_file)

print(f"Longest sequences saved to {output_file}.")