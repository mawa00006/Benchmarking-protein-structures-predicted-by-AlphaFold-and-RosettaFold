import os
#from tempfile import gettempdir
import biotite.database.rcsb as rcsb

# Path to the text file containing PDB IDs
pdb_id_file = 'group1_pdbs.txt'  # Replace with your file path

# Directory to save the downloaded mmCIF files
cif_dir = 'cif_files'
if not os.path.exists(cif_dir):
    os.makedirs(cif_dir)

# Read PDB IDs from the file
with open(pdb_id_file, 'r') as file:
    pdb_ids = [line.strip() for line in file.readlines()]

# Download each mmCIF file
for pdb_id in pdb_ids:
    try:
        # Fetch the mmCIF file
        cif_file_path = rcsb.fetch(pdb_id, 'cif', cif_dir)
        print(f"Downloaded: {cif_file_path}")
    except Exception as e:
        print(f"Error downloading {pdb_id}: {e}")