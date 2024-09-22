# DSSP installed from https://github.com/PDB-REDO/dssp

import os
import pickle

from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
from Bio.PDB import is_aa
import numpy as np

from tqdm import tqdm

def get_mask(string1, string2):
    # Create a mask of False for the whole main string
    mask = np.full(len(string2), False)

    # Find the positions of the substring in the main string
    index = string2.find(string1)

    # For some reason the longest structure did not fully match in some cases, by just taking the first n (length of
    # longest chain) amino acids this could be resolved
    if index == -1:
        index = 0

    # Set True for the matching positions
    mask[index:index + len(string1)] = True
    return mask

def get_longest_chain(id):

    sequence = []
    # Open the FASTA file
    with open(f'../../data/fasta_files/{id}.fasta', 'r') as file:
        # Loop through each line in the file
        for line in file:
            # Check if the line starts with '>' (Header -> skip)
            if line.startswith('>'):
                continue
            else:
                sequence.append(line)

    longest_chain = "".join(sequence)
    return longest_chain


def main():

    parser = MMCIFParser(QUIET=True)
    mmcif_path = "../../data/mmcif_files_longest_chain/"
    files = os.listdir(mmcif_path)

    for file in tqdm(files):
        if file.endswith(".cif"):
            #file = "8QES.cif"
            print(f"Running DSSP on {file}")
            structure = parser.get_structure('protein', os.path.join(mmcif_path, file))
            file_mmcif = file.replace("cif", "mmcif")
            structure_mmcif = parser.get_structure("protein", f"../../data/mmcif_files/{file_mmcif}")
            model = structure[0]
            chain_id = model.child_dict.keys()
            model_mmcif = structure_mmcif[0]#["A"]

            try:
                dssp = DSSP(model_mmcif, f"../../data/mmcif_files/{file_mmcif}")
            except Exception as e:
                print(e)
                continue

            secondary_structure = []
            amino_acids = []
            residues = dssp.property_list
            for residue in residues:
                secondary_structure.append(residue[2])
                amino_acids.append(residue[1])

            secondary_structure_string = "".join(secondary_structure)
            amino_acid_string = "".join(amino_acids)

            id = file.split(".")[0]
            longest_chain = get_longest_chain(id).replace("\n", "")
            x = len(longest_chain)
            mask = get_mask(longest_chain, amino_acid_string)

            filtered_string = ''.join([char for char, m in zip(secondary_structure_string, mask) if m])
            mask = [1 if c in ["T", "S", "-"] else 0 for c in filtered_string]
            if len(mask) == 0:
                print(0)

            with open(f"../../data/loop_annotations/{id}.pickle", "wb") as fp:
                pickle.dump(mask, fp)


if __name__ == "__main__":

    main()
