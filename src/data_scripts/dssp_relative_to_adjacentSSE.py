import os
import pickle
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.DSSP import DSSP
import numpy as np
from tqdm import tqdm

def get_mask(string1, string2):
    """
    Create a mask indicating where string1 (longest chain) matches within string2 (DSSP-annotated amino acid sequence).
    """
    mask = np.full(len(string2), False)
    index = string2.find(string1)

    if index == -1:
        index = 0

    mask[index:index + len(string1)] = True
    return mask

def get_longest_chain(id):
    """
    Retrieve the longest chain sequence from a FASTA file.
    """
    sequence = []
    with open(f'../../data/fasta_files/{id}.fasta', 'r') as file:
        for line in file:
            if line.startswith('>'):
                continue
            else:
                sequence.append(line.strip())

    longest_chain = "".join(sequence)
    return longest_chain

def get_residue_before_after_loop(loop_mask, secondary_structure):
    """
    Get secondary structure annotations for residues before and after each loop region.
    """
    loop_info = []
    in_loop = False
    start = end = None

    for i, mask_val in enumerate(loop_mask):
        if mask_val == 1 and not in_loop:
            in_loop = True
            start = i
        elif mask_val == 0 and in_loop:
            in_loop = False
            end = i - 1
            before = secondary_structure[start - 1] if start > 0 else "end"
            after = secondary_structure[end + 1] if end + 1 < len(secondary_structure) else "end"
            loop_info.append((start, end, before, after))

    if in_loop:
        before = secondary_structure[start - 1] if start > 0 else "end"
        after = "end"
        loop_info.append((start, len(loop_mask) - 1, before, after))

    return loop_info

def main():
    parser = MMCIFParser(QUIET=True)
    mmcif_path = "../../data/mmcif_files_longest_chain/"
    files = os.listdir(mmcif_path)

    for file in tqdm(files):
        if file.endswith(".cif"):
            print(f"Running DSSP on {file}")

            file_mmcif = file.replace("cif", "mmcif")
            structure_mmcif = parser.get_structure("protein", f"../../data/mmcif_files/{file_mmcif}")
            model_mmcif = structure_mmcif[0]

            try:
                dssp = DSSP(model_mmcif, f"../../data/mmcif_files/{file_mmcif}")
            except Exception as e:
                print(f"Error running DSSP for {file}: {e}")
                continue

            secondary_structure = []
            amino_acids = []
            residues = dssp.property_list

            for residue in residues:
                secondary_structure.append(residue[2])  # DSSP annotation
                amino_acids.append(residue[1])           # Amino acid

            secondary_structure_string = "".join(secondary_structure)
            amino_acid_string = "".join(amino_acids)

            id = file.split(".")[0]
            longest_chain = get_longest_chain(id).replace("\n", "")
            mask = get_mask(longest_chain, amino_acid_string)

            filtered_string = ''.join([char for char, m in zip(secondary_structure_string, mask) if m])

            loop_mask = [1 if c in ["T", "S", "-"] else 0 for c in filtered_string]

            if len(loop_mask) == 0:
                print(f"No loops found in {id}")
                continue

            # Get the secondary structure annotations for residues before and after each loop
            loop_info = get_residue_before_after_loop(loop_mask, filtered_string)

            # Prepare data for pickling: loop mask and DSSP annotations
            data_to_pickle = {
                'loop_mask': loop_mask,
                'loop_info': loop_info
            }

            # Save the data into a pickle file
            with open(f"../../data/loop_annotations_and_masks/{id}.pickle", "wb") as fp:
                pickle.dump(data_to_pickle, fp)

if __name__ == "__main__":
    main()
