
import biotite.structure as struc
import biotite.structure.io as strucio
import numpy as np


def calculate_rmsd(s1, s2):
    all_dist = 0
    for l in range(len(s1)):
        c1 = s1[l].coord
        c2 = s2[l].coord
        # Calculate the RMSD between two sets of coordinates
        distance = np.sqrt((c1[0] - c2[0])** 2 + (c1[1] - c2[1])** 2 + (c1[2] - c1[2])** 2)
        all_dist += distance

    rmsd = np.sqrt((all_dist/len(s1)))
    return rmsd

def main():

    file1 = "testrmsd/robetta_models_633493.pdb"
    #file1 = "testrmsd/fold_8ton_longest_chain_model_1.cif"
    file2 = "testrmsd/fold_8ton_longest_chain_model_1.cif"


    # Parse the PDB files and get structure
    structure1 = strucio.load_structure(file1)[0]
    structure2 = strucio.load_structure(file2)


    # only compare chain A of the two input files
    chain_id = 'A'
    structure1_chainA = structure1[structure1.chain_id == chain_id]
    structure2_chainA = structure2[structure2.chain_id == chain_id]

    structure_aligned = struc.superimpose(structure1_chainA, structure2_chainA)
    structure2_chainA = structure_aligned[0]

    if len(structure1_chainA) < len(structure2_chainA):
        rmsd = calculate_rmsd(structure1_chainA, structure2_chainA)
    else:
        rmsd = calculate_rmsd(structure2_chainA, structure1_chainA)

    print(f"RMSD: {rmsd:.3f} Å")

#if __name__ == "__main__":

main()
