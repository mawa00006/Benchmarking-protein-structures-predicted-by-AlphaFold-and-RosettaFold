import os

import biotite.structure as struc
import biotite.structure.io as strucio
import numpy as np
import pickle
from Bio.SVDSuperimposer import SVDSuperimposer
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBParser import PDBParser
from tqdm import tqdm
from biotite.structure import AtomArray

def calculate_rmsd(s1, s2):
    all_dist = 0
    for l in range(len(s1)):
        c1 = s1[l]#.coord
        c2 = s2[l]#.coord
        # Calculate the RMSD between two sets of coordinates
        distance = np.sqrt((c1[0] - c2[0])** 2 + (c1[1] - c2[1])** 2 + (c1[2] - c1[2])** 2)
        all_dist += distance

    rmsd = np.sqrt((all_dist/len(s1)))
    return rmsd

def get_loop_annotation(filepath):
    with open(filepath, "rb") as fp:
        annotation = pickle.load(fp)
    return annotation


def filter_amino_acids_by_mask(atom_array: AtomArray, mask: np.ndarray) -> list[AtomArray]:
    """
    Filters the amino acids from an AtomArray based on a binary mask.
    Groups the atoms of the amino acids in contiguous blocks of `1`s from the mask.

    :param atom_array: The AtomArray containing the atoms of the protein,
                       with atoms grouped by amino acids in order.
    :param mask: A binary mask (numpy array) where `1` indicates keeping the corresponding amino acid.
    :returns: A list of AtomArray objects, where each AtomArray contains atoms
              of the contiguous amino acids marked by `1` in the mask.
    :raises ValueError: If the length of the mask does not match the number of unique residues in the AtomArray.
    """

    residue_ids = np.unique(atom_array.res_id)  # Get unique residue IDs in sequence
    if len(mask) != len(residue_ids):
        raise ValueError(f"Mask length must match the number of unique residues in the AtomArray. {len(mask)}, {len(residue_ids)}")

    # Step 1: Filter by mask
    filtered_residues = residue_ids[mask == 1]  # Get residue IDs where mask is 1

    # Step 2: Group by contiguous '1's
    blocks = group_residues_by_mask(mask, residue_ids)

    # Step 3: Extract atoms for each contiguous block
    batches = [extract_atoms_for_block(atom_array, block) for block in blocks]

    return batches


def group_residues_by_mask(mask: np.ndarray, residue_ids: np.ndarray) -> list[np.ndarray]:
    """
    Groups contiguous amino acids (residues) where the mask is `1`.

    :param mask: A binary mask indicating which residues to keep.
    :param residue_ids: The residue IDs corresponding to the protein's amino acids.
    :returns: A list of arrays, each containing the residue IDs for one contiguous block of `1`s.
    """
    blocks = []
    current_block = []

    for i, val in enumerate(mask):
        if val == 1:
            current_block.append(residue_ids[i])
        else:
            if current_block:
                blocks.append(np.array(current_block))
                current_block = []

    if current_block:  # Append the last block if not empty
        blocks.append(np.array(current_block))

    return blocks


def extract_atoms_for_block(atom_array: AtomArray, block: np.ndarray) -> AtomArray:
    """
    Extracts all atoms corresponding to a given block of residue IDs.

    :param atom_array: The AtomArray containing the atoms of the protein.
    :param block: A list or array of residue IDs corresponding to a contiguous block.
    :returns: An AtomArray containing the atoms for the given block of residues.
    """
    return atom_array[np.isin(atom_array.res_id, block)]

def main():

    loop_path = "../../data/loop_annotations"
    loop_annotation_files = os.listdir(loop_path)

    for loop_annotation_file in loop_annotation_files: #tqdm(loop_annotation_files):

        compound_id = loop_annotation_file.split(".")[0]
        #loop_annotation = get_loop_annotation(os.path.join(loop_path, loop_annotation_file))

        try:
            experimental_structure = strucio.load_structure(f"../../data/mmcif_files_longest_chain/{compound_id}.mmcif")
            alphafold_structure = strucio.load_structure(f"../../data/alphafold_extracted/{compound_id}/FOLD_{compound_id}_MODEL_0.cif")
            rosetta_structure = strucio.load_structure(f"../../data/rosetta/robetta_models_{compound_id}.pdb")[0]
            rosetta_structure = rosetta_structure[rosetta_structure.element != "H"]
            loop_annotation = get_loop_annotation(os.path.join(loop_path, loop_annotation_file))
        except Exception as e:
            #print(e)
            continue

        #experimental_loops = filter_amino_acids_by_mask(experimental_structure, loop_annotation)
        try:
            alphafold_loops = filter_amino_acids_by_mask(alphafold_structure, loop_annotation)
            rosetta_loops = filter_amino_acids_by_mask(rosetta_structure, loop_annotation)
        except Exception as e:
            print("Filtering:",compound_id)
            continue

        q = alphafold_loops
        b = rosetta_loops
        x_experimental = experimental_structure.coord
        y = structure1_chainA.coord[structure1_chainA.element != "H"]
        sup = SVDSuperimposer()
        sup.set(x, y)
        sup.run()
        y_on_x = sup.get_transformed()

        rmsd = calculate_rmsd(x, y_on_x)
        print(f"RMSD: {rmsd:.3f} Å")







    file1 = "../../data/rmsd_test/robetta_models_633493.pdb"
    #file1 = "testrmsd/fold_8ton_longest_chain_model_1.cif"
    file2 = "../../data/rmsd_test/fold_8ton_longest_chain_model_1.cif"


    # Parse the PDB files and get structure
    structure1 = strucio.load_structure(file1)[0]
    structure2 = strucio.load_structure(file2)


    # only compare chain A of the two input files
    chain_id = 'A'
    structure1_chainA = structure1[structure1.chain_id == chain_id]
    structure2_chainA = structure2[structure2.chain_id == chain_id]


    x = structure2_chainA.coord
    y = structure1_chainA.coord[structure1_chainA.element != "H"]
    sup = SVDSuperimposer()
    sup.set(x,y)
    sup.run()
    y_on_x = sup.get_transformed()

    rmsd = calculate_rmsd(x, y_on_x)
    print(f"RMSD: {rmsd:.3f} Å")

    structure_aligned = struc.superimpose(structure1_chainA, structure2_chainA)
    structure2_chainA = structure_aligned[0]

    if len(structure1_chainA) < len(structure2_chainA):
        rmsd = calculate_rmsd(structure1_chainA, structure2_chainA)
    else:
        rmsd = calculate_rmsd(structure2_chainA, structure1_chainA)

    print(f"RMSD: {rmsd:.3f} Å")

if __name__ == "__main__":

    main()
