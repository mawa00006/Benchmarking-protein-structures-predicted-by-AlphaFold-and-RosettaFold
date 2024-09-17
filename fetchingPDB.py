"""
This script is used for fetching the given pdb entries and deposit them as .cif files in a new folder
"""

import os
from Bio.PDB import PDBList


def download_mmCIF_files(pdb_ids, output_folder):
    """Download mmCIF files for the given PDB IDs and save them in the specified folder."""
    pdb_list = PDBList()

    # Make sure the output directory exists
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    for pdb_id in pdb_ids:
        print(f"Downloading mmCIF file for {pdb_id}...")

        # Download the mmCIF file
        mmcif_file_path = pdb_list.retrieve_pdb_file(pdb_id, pdir=output_folder, file_format='mmCif')

        # Rename the file to use .cif extension
        new_file_name = os.path.join(output_folder, f"{pdb_id.upper()}.cif")
        os.rename(mmcif_file_path, new_file_name)
        print(f"Saved as {new_file_name}")


def read_pdb_ids(file_path):
    """Read the list of PDB IDs from a text file."""
    with open(file_path, 'r') as f:
        pdb_ids = [line.strip().upper() for line in f.readlines()]
    return pdb_ids


# Input file containing PDB IDs and output folder for storing mmCIF files
pdb_file_path = 'group1_pdbs.txt'  # Path to your file with PDB IDs
output_folder = 'mmcif_files'  # Directory where mmCIF files will be stored

# Step 1: Read PDB IDs from the text file
pdb_ids = read_pdb_ids(pdb_file_path)

# Step 2: Download all mmCIF files into the specified folder
download_mmCIF_files(pdb_ids, output_folder)

print(f"Downloaded {len(pdb_ids)} mmCIF files into '{output_folder}' folder.")
