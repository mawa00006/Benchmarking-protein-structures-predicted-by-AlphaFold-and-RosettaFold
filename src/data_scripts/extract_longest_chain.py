"""Extracts the longest chain from pdbx/mmcif files and saves them as FASTA files."""
import os
from Bio import PDB, SeqIO
from Bio.Seq import Seq
from Bio.PDB import is_aa, MMCIFIO
from Bio.PDB.Chain import Chain
from Bio.PDB.Structure import Structure
from Bio.SeqRecord import SeqRecord
from Bio.PDB.MMCIFParser import MMCIFParser
from typing import Tuple, Optional, Union

# Map 3-letter amino-acid code to 1-letter
dict_aa = {'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 'GLY': 'G', 'HIS': 'H',
           'ILE': 'I', 'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q',
           'ARG': 'R', 'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'}


def parse_structure(parser: MMCIFParser, pdbx_file: str) -> Union[Structure, None]:
    """
    Parse and return the structure from an mmCIF (pdbx) file.

    :param parser: The parser instance, e.g., Bio.PDB.MMCIFParser, for parsing pdbx/mmCIF files.
    :param pdbx_file: Path to the pdbx/mmCIF file to be parsed.
    :returns: The parsed structure object or None if parsing fails.
    :raises FileNotFoundError: If the pdbx/mmCIF file does not exist.
    :raises ValueError: If the pdbx/mmCIF file cannot be parsed due to format errors.
    """
    try:
        structure = parser.get_structure('protein', pdbx_file)
    except FileNotFoundError:
        raise FileNotFoundError(f"File not found: {pdbx_file}")
    except Exception as e:

        raise ValueError(f"Failed to parse the pdbx/mmCIF file: {e}")

    return structure


def get_longest_chain(pdbx_file: str) -> Tuple[Optional[Chain], int]:
    """
    Retrieve the longest chain from an mmCIF (pdbx) file.

    :param pdbx_file: Path to the mmCIF file.
    :returns: A tuple containing the longest chain and its length.
    """

    # Initialize the MMCIFParser to parse pdbx/mmCIF files
    parser = MMCIFParser(QUIET=True)

    # Get structure
    structure = parse_structure(parser, pdbx_file)

    longest_chain = None
    longest_chain_length = 0

    # Loop through all chains in the structure
    for model in structure:
        for chain in model:
            # Get all residues in the chain and filter only standard amino acids
            residues = [residue for residue in chain if PDB.is_aa(residue, standard=True)]
            chain_length = len(residues)

            # If this chain is the longest, store it
            if chain_length > longest_chain_length:
                longest_chain_length = chain_length
                longest_chain = chain

    return longest_chain, longest_chain_length


def three_to_one(three_letter_code: str) -> str:
    """
    Convert a 3-letter amino acid code to a 1-letter code.
    If the code is unknown, it returns 'A' for Alanine.

    :param three_letter_code: The 3-letter code of the amino acid (e.g., "ALA").
    :returns: The 1-letter amino acid code (e.g., "A").
    """
    return dict_aa.get(three_letter_code, 'A')


def get_chain_sequence(chain: Chain) -> str:
    """
    Extract the amino acid sequence of the given chain in FASTA format.

    :param chain: The chain object from which to extract the sequence (Bio.PDB.Chain).
    :returns: The amino acid sequence as a string, using 1-letter codes.
    """
    sequence = "".join(
        [three_to_one(residue.resname) for residue in chain.get_residues() if is_aa(residue, standard=True)]
    )
    return sequence


def save_longest_chain_as_fasta(pdbx_file: str, output_dir_pdbx: str, longest_chain) -> None:
    """
    Extract the longest chain from a structure and save it as a FASTA file.

    :param pdbx_file: The name of the input pdbx (mmCIF) file.
    :param output_dir_pdbx: The directory where the output FASTA file will be saved.
    :param longest_chain: The chain object representing the longest chain in the structure.
    :returns: The path to the saved FASTA file or None if no valid chain is found.
    :raises ValueError: If no valid chain is found in the pdbx file
    """
    if longest_chain:
        # Get the sequence of the longest chain
        sequence = get_chain_sequence(longest_chain)
        record = SeqRecord(Seq(sequence), id=longest_chain.id, description=f"Longest chain from {pdbx_file}")

        # Write the sequence to a FASTA file
        name = pdbx_file.split(".")[0]
        fasta_file = os.path.join(output_dir_pdbx, f"{name}.fasta")
        SeqIO.write(record, fasta_file, "fasta")
        print(f"Saved FASTA: {fasta_file}")
        return None
    else:
        raise ValueError(f"No valid chain found in {pdbx_file}")


def save_longest_chain_as_mmcif(chain: Chain, output_dir: str, structure_id:str) -> None:
    """
    Save a given chain object to an MMCIF file.

    :param chain: The chain object to save in MMCIF format.
    :returns: The path to the saved MMCIF file.
    :raises ValueError: If the chain object is invalid.
    """
    if not isinstance(chain, Chain):
        raise ValueError("The provided chain is not a valid Chain object.")

    # Initialize MMCIFIO for saving in mmcif format
    io = MMCIFIO()

    # Set the chain structure to MMCIF
    io.set_structure(chain)

    # Generate the output filename
    output_file = os.path.join(output_dir,f"{structure_id}.mmcif")

    # Save the chain in MMCIF format
    io.save(output_file)


# Directory where pdbx/mmcif files are stored
pdbx_dir = "../../data/pdbx_files"
# Directory to store FASTA files containing the longest chain
output_dir_fasta = "../../data/fasta_files"
# Directory to store mmcif files containing the longest chain
output_dir_mmcif = "../../data/mmcif_files_longest_chain"

# Create the output directory if it doesn't exist
if not os.path.exists(output_dir_fasta):
    os.makedirs(output_dir_fasta)

if not os.path.exists(output_dir_mmcif):
    os.makedirs(output_dir_mmcif)

for pdbx_file in os.listdir(pdbx_dir):
    if pdbx_file.endswith(".pdbx"):
        pdb_id = pdbx_file.replace(".pdbx", "")
        cif_path = os.path.join(pdbx_dir, pdbx_file)
        longest_chain, length = get_longest_chain(cif_path)
        print(f"File: {pdbx_file}, Longest chain: {longest_chain.id}, Length: {length}")

        # Save the sequence of the longest chain in FASTA format
        save_longest_chain_as_fasta(pdbx_file, output_dir_fasta, longest_chain)
        save_longest_chain_as_mmcif(longest_chain, output_dir_mmcif, pdb_id)
