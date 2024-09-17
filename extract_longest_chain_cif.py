from Bio import PDB
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import os


dict_aa = {'ALA':'A', 'CYS':'C', 'ASP':'D', 'GLU':'E', 'PHE':'F', 'GLY':'G', 'HIS':'H', 'ILE':'I',
           'LYS':'K', 'LEU':'L', 'MET':'M', 'ASN':'N', 'PRO':'P', 'GLN':'Q', 'ARG':'R', 'SER':'S', 
           'THR':'T', 'VAL':'V', 'TRP':'W', 'TYR':'Y'}


# Function to get the sequence of the longest chain in the cif file
def get_longest_chain(cif_file):
    parser = PDB.MMCIFParser(QUIET=True)  # Parser for cif files
    structure = parser.get_structure('protein', cif_file)

    longest_chain = None
    max_length = 0

    # Loop through all chains in the structure
    for model in structure:
        for chain in model:
            # Get all residues in the chain and filter only standard amino acids
            residues = [residue for residue in chain if PDB.is_aa(residue, standard=True)]
            chain_length = len(residues)

            # If this chain is the longest, store it
            if chain_length > max_length:
                max_length = chain_length
                longest_chain = chain

    return longest_chain, max_length


# Function to 3-letter amino acid code to 1-letter code
def three_to_one(three_letter_code):
    try:
        return dict_aa[three_letter_code]
    except KeyError:
        return 'X' #Unknown amino acids to be represented as 'X


# Function to extract the sequence of the longest chain in FASTA format
def get_chain_sequence(chain):
    seq = "".join([three_to_one(residue.resname) 
                   for residue in chain.get_residues() if PDB.is_aa(residue, standard=True)])
    return seq


# Function to extract and save the longest chain as a FASTA file
def save_longest_chain_as_fasta(cif_file, output_dir_cif,longest_chain):
    #longest_chain, _ = get_longest_chain(cif_file)
    if longest_chain:
        sequence = get_chain_sequence(longest_chain)
        record = SeqRecord(Seq(sequence), id=longest_chain.id, description=f"Longest chain from {cif_file}")

        # Write the sequence to a FASTA file
        fasta_file = os.path.join(output_dir_cif, f"{cif_file}_{longest_chain.id}_longest_chain.fasta")
        SeqIO.write(record, fasta_file, "fasta")
        print(f"Saved FASTA: {fasta_file}")
        return fasta_file
    else:
        print(f"No valid chain found in {cif_file}")
        return None



cif_dir = "cif_files"  # Directory where cif files are stored
output_dir_cif = "fasta_files_from_cif"  # Directory to store FASTA files

# Create the output directory if it doesn't exist
if not os.path.exists(output_dir_cif):
    os.makedirs(output_dir_cif)


for cif_file in os.listdir(cif_dir):
    if cif_file.endswith(".cif"):
        cif_path = os.path.join(cif_dir, cif_file)
        longest_chain, length = get_longest_chain(cif_path)
        print(f"File: {cif_file}, Longest chain: {longest_chain.id}, Length: {length}")
       
        # Save the sequence of the longest chain in FASTA format
        fasta_file = save_longest_chain_as_fasta(cif_file, output_dir_cif,longest_chain)

