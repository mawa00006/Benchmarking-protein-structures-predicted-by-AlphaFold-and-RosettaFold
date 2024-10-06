""" This script gets the secondary structure annotations for residues in a protein structure. """
# DSSP installed from https://github.com/PDB-REDO/dssp
import os
import pickle
import numpy as np

from tqdm import tqdm
from Bio.PDB.DSSP import DSSP
from Bio.PDB.MMCIFParser import MMCIFParser


def get_mask(longest_chain: str, all_chains: str) -> np.ndarray:
    """
    Create a mask indicating where the longest chain matches within the whole amino acid sequence.

    :param longest_chain: The longest chain sequence.
    :param all_chains: The whole structure amino acid sequence.
    :returns: A boolean mask array where True indicates matching positions.
    """
    mask = np.full(len(all_chains), False)

    # Find the positions of the substring in the main string
    index = all_chains.find(longest_chain)

    # For some reason the longest structure did not fully match in some cases, by just taking the first n (length of
    # longest chain) amino acids this could be resolved
    if index == -1:
        index = 0

    # Set True for the matching positions
    mask[index:index + len(longest_chain)] = True
    return mask


def get_longest_chain(id: str) -> str:
    """
    Retrieve the longest chain sequence from a FASTA file.

    :param id: The identifier of the protein, used to locate the corresponding FASTA file.
    :returns: The longest chain sequence as a string.
    """
    sequence = []
    with open(f'../../data/fasta_files/{id}.fasta', 'r') as file:
        for line in file:
            # Skip first line (header)
            if line.startswith('>'):
                continue
            else:
                sequence.append(line.strip())
    # Join the sequence lines into a single string
    longest_chain = "".join(sequence)
    return longest_chain


def get_residue_before_after_loop(loop_mask: np.ndarray, secondary_structure: str) -> list[tuple[int, int, str, str]]:
    """
    Get secondary structure annotations for residues before and after each loop region.

    :param loop_mask: A mask array indicating loop regions.
    :param secondary_structure: A string of secondary structure annotations.
    :returns: A list of tuples containing start and end indices of loops, and the secondary structure annotations before and after each loop.
    """
    loop_info = []
    in_loop = False
    start = None

    # Loop through the mask array
    for i, mask_val in enumerate(loop_mask):
        # If the mask value is 1 and not in a loop, start a new loop
        if mask_val == 1 and not in_loop:
            in_loop = True
            start = i
        # If the mask value is 0 and in a loop, end the loop
        elif mask_val == 0 and in_loop:
            in_loop = False
            end = i - 1
            # Get the secondary structure annotations before and after the loop
            before = secondary_structure[start - 1] if start > 0 else "end"
            after = secondary_structure[end + 1] if end + 1 < len(secondary_structure) else "end"
            # Store the loop information
            loop_info.append((start, end, before, after))
    # If the loop ends at the end of the sequence, store the loop information
    if in_loop:
        before = secondary_structure[start - 1] if start > 0 else "end"
        after = "end"
        loop_info.append((start, len(loop_mask) - 1, before, after))

    return loop_info


def main():
    # Initialize the MMCIFParser
    parser = MMCIFParser(QUIET=True)
    # Get all CIF files in the directory (for identified longest chains)
    cif_path = "../../data/cif_files_longest_chain/"
    files = os.listdir(cif_path)

    # Loop through all CIF files
    for file in tqdm(files):
        if file.endswith(".cif"):
            print(f"Running DSSP on {file}")
            # Load the whole structure from the mmCIF file
            file_mmcif = file.replace("cif", "mmcif")
            structure_mmcif = parser.get_structure("protein", f"../../data/mmcif_files/{file_mmcif}")
            model_mmcif = structure_mmcif[0]
            # Run DSSP on the structure
            try:
                dssp = DSSP(model_mmcif, f"../../data/mmcif_files/{file_mmcif}")
            except Exception as e:
                print(f"Error running DSSP for {file}: {e}")
                continue
            # Initialize lists to store secondary structure annotations and amino acids
            secondary_structure = []
            amino_acids = []
            # Get the secondary structure annotations and amino acids for each residue
            residues = dssp.property_list

            for residue in residues:
                secondary_structure.append(residue[2])  # DSSP annotation
                amino_acids.append(residue[1])  # Amino acid

            # Join the secondary structure and amino acid lists into strings
            secondary_structure_string = "".join(secondary_structure)
            amino_acid_string = "".join(amino_acids)

            # Get the longest chain sequence
            id = file.split(".")[0]
            longest_chain = get_longest_chain(id).replace("\n", "")
            # Get a mask indicating where the longest chain matches within the whole amino acid sequence
            mask = get_mask(longest_chain, amino_acid_string)

            # Filter the secondary structure annotations based on the mask
            filtered_string = ''.join([char for char, m in zip(secondary_structure_string, mask) if m])
            # Create a mask for the loop regions
            loop_mask = [1 if c in ["T", "S", "-"] else 0 for c in filtered_string]
            # If no loops are found, print a message and continue with the next file
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
