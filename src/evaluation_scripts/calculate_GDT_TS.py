
import biotite.structure as struc
import biotite.structure.io as strucio

import numpy as np

def get_GDT(s1, s2, d):
    gdt = 0
    num_of_residues = 0
    for l in range(len(s1)):
        if s1[l].atom_name == "CA" and s1[l].atom_name == "CA":
            c1 = s1[l].coord
            c2 = s2[l].coord

            # Calculate the gdt between two sets of coordinates
            distance = np.sqrt((c1[0] - c2[0])** 2 + (c1[1] - c2[1])** 2 + (c1[2] - c2[2]) ** 2)
            if distance < d:
                gdt += 1
            num_of_residues += 1
        elif s1[l].atom_name != s2[l].atom_name:
            if s1[l].atom_name == "CA" or s1[l].atom_name == "CA":
                print(f"ERROR: {s1[l].atom_name} != {s2[l].atom_name}")
                #pass

    gdt = (gdt / num_of_residues) * 100

    return gdt


def get_GDT_TS(s1, s2):

    gdt_ts = 0

    #GDT_TS = (GDT_P1 + GDT_P2 + GDT_P4 + GDT_P8)/4
    distances = [1, 2, 4, 8]

    for d in distances:
        gdt_ts += get_GDT(s1, s2, d)

    gdt_ts = gdt_ts / 4

    return gdt_ts

def main():

    # Parse the cif files and get structure
    structure1 = strucio.load_structure(r"C:\Users\flori\Desktop\InfoZeug\Python_Projects\SSBIGroupProject\testrmsd\test_robetta.cif")
    structure2 = strucio.load_structure(r"C:\Users\flori\Desktop\InfoZeug\Python_Projects\SSBIGroupProject\testrmsd\fold_8ton_longest_chain_model_1.cif")


    structure_aligned = struc.superimpose(structure1, structure2)
    structure2 = structure_aligned[0]

    gdt_ts = get_GDT_TS(structure1, structure2)

    print(f"GDT_TS: {gdt_ts:.3f} %")

main()