"""Generate a JSON file from a folder of FASTA files for batch uploading to the AlphaFold3 webserver."""
import os
import json
from Bio import SeqIO


def create_JSON(output_json, folder_path):
    """
    Generate a JSON file from a folder of FASTA files for batch uploading to the AlphaFold3 webserver.

    :param output_json: The file path where the generated JSON file will be saved.
    :param folder_path: The directory containing the FASTA files to process.
    :returns: None. Writes the generated JSON content to the specified output file.
    :raises FileNotFoundError: If the specified folder or output file path is invalid.
    :raises ValueError: If the FASTA file cannot be read or is in an incorrect format.
    """

    # List for collecting the entries
    list_of_entries = []

    # Going through all files
    for file_name in os.listdir(folder_path):
        if file_name.endswith(".fasta"):
            file_path = os.path.join(folder_path, file_name)
            print(f"Processing {file_name}...")

            # Read Sequence from the FASTA file
            record = SeqIO.read(file_path, "fasta")

            # Get protein name (first 4 characters of the file name)
            name = str(file_name)[:4]
            sequence = str(record.seq)

            # Structure of the JSON file
            dictionary = {
                "name": name,
                "modelSeeds": [],
                "sequences": [
                    {
                        "proteinChain": {
                            "sequence": sequence,
                            "count": 1
                        }
                    }
                ]
            }

            list_of_entries.append(dictionary)

    # Convert the list of entries to JSON format
    json_object = json.dumps(list_of_entries, indent=4)

    # Write the JSON object to the output file
    with open(output_json, 'w') as j_out:
        j_out.write(json_object)


# Create JSON file containing 20 entries to upload to the AlphaFold3 webserver
createJSON("../../data/json_files/fastasplit_1_20.json", "../../data/fasta_files/Fasta_longest_chain_1_20")
createJSON("../../data/json_files/fastasplit_21_40.json", "../../data/fasta_files/Fasta_longest_chain_21_40")
createJSON("../../data/json_files/fastasplit_41_60.json", "../../data/fasta_files/Fasta_longest_chain_41_60")
createJSON("../../data/json_files/fastasplit_61_80.json", "../../data/fasta_files/Fasta_longest_chain_61_80")
createJSON("../../data/json_files/fastasplit_81_90.json", "../../data/fasta_files/Fasta_longest_chain_81_90")
