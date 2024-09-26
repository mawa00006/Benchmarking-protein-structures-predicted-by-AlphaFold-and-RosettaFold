import logging
import os
import pickle
import traceback

import biotite.structure as struc
import biotite.structure.io as strucio
import numpy as np
from Bio.SVDSuperimposer import SVDSuperimposer


compound_id = "8IDW"


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



try:
    # Load the structures
    experimental_structure = strucio.load_structure(f"C:/Users/flori/Desktop/alignmentFigure/pdb/{compound_id}.cif")
    alphafold_structure = strucio.load_structure(
        f"C:/Users/flori/Desktop/alignmentFigure/alpha/FOLD_{compound_id}_MODEL_0.cif")
    rosetta_structure = strucio.load_structure(f"C:/Users/flori/Desktop/alignmentFigure/rosetta/robetta_models_{compound_id}.pdb")[0]
    rosetta_structure = rosetta_structure[np.isin(rosetta_structure.atom_name, ["C", "CA", "N", "O"])]  # Only keep backbone atoms
    experimental_structure = experimental_structure[np.isin(experimental_structure.atom_name, ["C", "CA", "N", "O"])]
    alphafold_structure = alphafold_structure[np.isin(alphafold_structure.atom_name, ["C", "CA", "N", "O"])]
except FileNotFoundError as fnf_error:
    logging.error(f"File not found for compound {compound_id}: {fnf_error}")
except Exception as e:
    logging.error(f"Error loading structures or annotations for compound {compound_id}: {e}")
    logging.debug(traceback.format_exc())

#strucio.save_structure(f"{compound_id}_rosetta_backbone.cif", rosetta_structure)
#strucio.save_structure(f"{compound_id}_alpha_backbone.cif", alphafold_structure)



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
    
    # Set atom names and residue IDs (optional: chain_id, atom_type, etc.)

    
    return atom_array


rosetta_structure.coord = y_on_x_ros
alphafold_structure.coord = y_on_x_alpha

strucio.save_structure(f"C:/Users/flori/Desktop/alignmentFigure/{compound_id}_alpha_superimposed.cif", alphafold_structure)
strucio.save_structure(f"C:/Users/flori/Desktop/alignmentFigure/{compound_id}_rosetta_superimposed.cif", rosetta_structure)

with open(f'annotations/{compound_id}.pickle', 'rb') as f:
    x = pickle.load(f)

    loopcounter = 1
    inloop = False
    for i in range(len(x)):


        if x[i] == 1 and inloop == False:
            inloop = True


            print("Loop ", loopcounter, " starts at ", i+1, end="")
            loopcounter += 1
        elif x[i] == 0 and inloop == True:
            inloop = False
            print(", ends at ", i)

