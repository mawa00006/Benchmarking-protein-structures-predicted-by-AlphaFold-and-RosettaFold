""" Extracts the alignment of a compound from the AlphaFold and Rosetta models. """
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

compound_id = "8TLE"


try:
    # Load the structures
    experimental_structure = strucio.load_structure(f"C:/Users/flori/Desktop/alignmentFigure/pdb chains/{compound_id}.cif")
    alphafold_structure = strucio.load_structure(
        f"C:/Users/flori/Desktop/alignmentFigure/alphafold extracted/FOLD_{compound_id}_MODEL.cif")
    rosetta_structure = strucio.load_structure(f"C:/Users/flori/Desktop/alignmentFigure/rosetta/robetta_models_{compound_id}.pdb")[0]
    rosetta_structure = rosetta_structure[np.isin(rosetta_structure.atom_name, ["C", "CA", "N", "O"])]  # Only keep backbone atoms
    experimental_structure = experimental_structure[np.isin(experimental_structure.atom_name, ["C", "CA", "N", "O"])]
    alphafold_structure = alphafold_structure[np.isin(alphafold_structure.atom_name, ["C", "CA", "N", "O"])]
except FileNotFoundError as fnf_error:
    logging.error(f"File not found for compound {compound_id}: {fnf_error}")
except Exception as e:
    logging.error(f"Error loading structures or annotations for compound {compound_id}: {e}")
    logging.debug(traceback.format_exc())


# Superimpose structures
x_experimental = experimental_structure.coord
y_alpha = alphafold_structure.coord
y_ros = rosetta_structure.coord

# Superimpose structures
y_on_x_alpha = superimpose_structures(x_experimental, y_alpha)
y_on_x_ros = superimpose_structures(x_experimental, y_ros)

def createStruc(coordis):
    num_atoms = len(coordis)

    # Create an AtomArray with the number of atoms
    atom_array = struc.AtomArray(num_atoms)

    # Set the coordinates
    atom_array.coord = coordis
    """
    # Set atom names and residue IDs (optional: chain_id, atom_type, etc.)
    atom_array.atom_name = np.array(atom_names)
    atom_array.res_id = np.array(res_ids)
    atom_array.chain_id = np.array([chain_id] * num_atoms)
    """
    return atom_array

alpha_struc = createStruc(y_on_x_alpha)
rosetta_struc = createStruc(y_on_x_ros)

# Save the superimposed structures
strucio.save_structure(f"{compound_id}_alpha_superimposed.cif", alpha_struc)
strucio.save_structure(f"{compound_id}_rosetta_superimposed.cif", rosetta_struc)
