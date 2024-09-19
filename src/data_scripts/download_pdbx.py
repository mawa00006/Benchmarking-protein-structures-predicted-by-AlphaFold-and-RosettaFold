"""This script is used for fetching the given pdb entries and deposit them as .pdbx files in a new folder"""
import os
import biotite.database.rcsb as rcsb

# Path to the text file containing PDB IDs
pdb_id_file_path = '../../data/group1_pdbs.txt'  # Replace with your file path

# Directory to save the downloaded pdbx files
pdbx_dir = '../../data/pdbx_files'
if not os.path.exists(pdbx_dir):
    os.makedirs(pdbx_dir)

# Read PDB IDs from the file
with open(pdb_id_file_path, 'r') as file:
    pdb_ids = [line.strip() for line in file.readlines()]

# Download each pdbx file
for pdb_id in pdb_ids:
    try:
        # Fetch the pdbx file
        pdbx_file_path = rcsb.fetch(pdb_id, 'pdbx', pdbx_dir)
        print(f"Downloaded: {pdbx_file_path}")
    except Exception as e:
        print(f"Error downloading {pdb_id}: {e}")
