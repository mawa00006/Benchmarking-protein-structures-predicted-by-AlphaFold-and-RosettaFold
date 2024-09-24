"""
This script compares the RMSD of loop regions in protein structures
from experimental, AlphaFold, and Rosetta models.
"""

import logging
import os
import pickle
import traceback

import biotite.structure as struc
import biotite.structure.io as strucio
import numpy as np
import pandas as pd
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBParser import PDBParser
from Bio.SVDSuperimposer import SVDSuperimposer
from tqdm import tqdm
from biotite.structure import AtomArray

# Set up logging configuration
logging.basicConfig(filename="log.log",
                    filemode='w',
                    level=logging.ERROR,
                    format='%(asctime)s - %(levelname)s - %(message)s')


def calculate_rmsd(s1: np.ndarray, s2: np.ndarray) -> float:
    """
    Calculate the RMSD between two sets of coordinates.

    :param s1: First set of coordinates.
    :param s2: Second set of coordinates.
    :returns: RMSD value.
    """
    total_distance = 0
    for l in range(len(s1)):
        c1 = s1[l]
        c2 = s2[l]
        # Calculate the distance between corresponding atoms
        distance = np.sqrt((c1[0] - c2[0]) ** 2 +
                           (c1[1] - c2[1]) ** 2 +
                           (c1[2] - c1[2]) ** 2)
        total_distance += distance

    rmsd = np.sqrt(total_distance / len(s1))
    return rmsd


def get_loop_annotation(filepath: str) -> np.ndarray:
    """
    Load loop annotation from a file.

    :param filepath: Path to the file containing loop annotation.
    :returns: A numpy array with the loop annotation.
    """
    try:
        with open(filepath, "rb") as fp:
            annotation = pickle.load(fp)
        return annotation
    except Exception as e:
        logging.error(f"Error loading loop annotation from {filepath}: {e}")
        raise


def filter_amino_acids_by_mask(atom_array: AtomArray, mask: np.ndarray) -> list[AtomArray]:
    """
    Filters the amino acids from an AtomArray based on a binary mask.

    :param atom_array: The AtomArray containing the atoms of the protein,
                       with atoms grouped by amino acids in order.
    :param mask: A binary mask where `1` indicates keeping the corresponding amino acid.
    :returns: A list of AtomArray objects, where each AtomArray contains atoms
              of the contiguous amino acids marked by `1` in the mask.
    :raises ValueError: If the mask length doesn't match the residue count.
    """
    residue_ids = np.unique(atom_array.res_id)  # Get unique residue IDs in sequence
    if len(mask) != len(residue_ids):
        raise ValueError(f"Mask length must match the number of unique residues in the AtomArray. {len(mask)},"
                         f" {len(residue_ids)}")

    blocks = group_residues_by_mask(mask, residue_ids)
    return [extract_atoms_for_block(atom_array, block) for block in blocks]


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


def superimpose_structures(x_coords: np.ndarray, y_coords: np.ndarray) -> np.ndarray:
    """
    Superimpose two structures using SVD.

    :param x_coords: Coordinates of the reference structure.
    :param y_coords: Coordinates of the structure to be superimposed.
    :returns: Transformed coordinates of y structure after superimposition.
    """
    sup = SVDSuperimposer()
    sup.set(x_coords, y_coords)
    sup.run()
    return sup.get_transformed()


def main() -> None:
    """
    Main function to process structures, filter loops, and calculate RMSD for protein structures.
    """
    loop_path = "../../data/loop_annotations"
    loop_annotation_files = os.listdir(loop_path)

    out_df = pd.DataFrame(columns=['structure_id', 'model', 'loop_len', 'rmsd'])

    for loop_annotation_file in tqdm(loop_annotation_files):

        compound_id = loop_annotation_file.split(".")[0]

        try:
            # Load the structures
            experimental_structure = strucio.load_structure(f"../../data/mmcif_files_longest_chain/{compound_id}.cif")
            alphafold_structure = strucio.load_structure(
                f"../../data/alphafold_extracted/{compound_id}/FOLD_{compound_id}_MODEL_0.cif")
            rosetta_structure = strucio.load_structure(f"../../data/rosetta/robetta_models_{compound_id}.pdb")[0]
            rosetta_structure = rosetta_structure[np.isin(rosetta_structure.atom_name, ["C", "CA", "N", "O"])]  # Only keep backbone atoms
            experimental_structure = experimental_structure[np.isin(experimental_structure.atom_name, ["C", "CA", "N", "O"])]
            alphafold_structure = alphafold_structure[np.isin(alphafold_structure.atom_name, ["C", "CA", "N", "O"])]
            loop_annotation = get_loop_annotation(os.path.join(loop_path, loop_annotation_file))

        except FileNotFoundError as fnf_error:
            logging.error(f"File not found for compound {compound_id}: {fnf_error}")
            continue
        except Exception as e:
            logging.error(f"Error loading structures or annotations for compound {compound_id}: {e}")
            logging.debug(traceback.format_exc())
            continue

        try:
            # Filter loops
            experimental_loops = filter_amino_acids_by_mask(experimental_structure, loop_annotation)
            alphafold_loops = filter_amino_acids_by_mask(alphafold_structure, loop_annotation)
            rosetta_loops = filter_amino_acids_by_mask(rosetta_structure, loop_annotation)

        except KeyError as ke:
            logging.error(f"Key error while filtering loops for compound {compound_id}: {ke}")
            logging.debug(traceback.format_exc())
            continue
        except Exception as e:
            logging.error(f"Error filtering loops for compound {compound_id}: {e}")
            logging.debug(traceback.format_exc())
            continue

        for i, (exp_loop, alpha_loop, ros_loop) in enumerate(zip(experimental_loops, alphafold_loops, rosetta_loops)):
            loop_len = len(np.unique(exp_loop.res_id))
            x_experimental = exp_loop.coord
            y_alpha = alpha_loop.coord
            y_ros = ros_loop.coord

            # Superimpose structures
            y_on_x_alpha = superimpose_structures(x_experimental, y_alpha)
            y_on_x_ros = superimpose_structures(x_experimental, y_ros)

            rmsd_alpha = calculate_rmsd(x_experimental, y_on_x_alpha)
            rmsd_ros = calculate_rmsd(x_experimental, y_on_x_ros)

            alpha_row = {'structure_id': compound_id, 'model': "alpha", 'loop_len': loop_len,
                         'rmsd': round(rmsd_alpha, 3)}
            rosetta_row = {'structure_id': compound_id, 'model': "rosetta", 'loop_len': loop_len,
                           'rmsd': round(rmsd_ros, 3)}
            # Append the row to the DataFrame
            out_df = pd.concat([out_df, pd.DataFrame([alpha_row])], ignore_index=True)
            out_df = pd.concat([out_df, pd.DataFrame([rosetta_row])], ignore_index=True)

    out_df.to_csv("results.csv")


if __name__ == "__main__":
    main()
