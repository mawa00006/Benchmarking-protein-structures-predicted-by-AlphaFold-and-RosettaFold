"""
This script creates a JSON file for the AlphaFold3 webserver for batch upload.
Input: Foldername of fasta files which you want to upload
Output: JSON file in the right formatting for upload
"""

import json
import os
from Bio import SeqIO


def createJSON(output_json, folder_path):
    # List for collecting the entries
    listOfEntries = []

    # Going through all files
    for file_name in os.listdir(folder_path):
        if file_name.endswith(".fasta"):
            file_path = os.path.join(folder_path, file_name)
            print(f"Processing {file_name}...")

            # Read Sequence
            record = SeqIO.read(file_path, "fasta")

            # Get protein name
            name = str(file_name)[0:4]
            sequence = str(record.seq)

            # structure of the JSON file.
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

            listOfEntries.append(dictionary)

    json_object = json.dumps(listOfEntries, indent=4)

    with open(output_json, 'w') as j_out:
        j_out.write(json_object)


# Create JSON file containing 20 entries to upload to Alphafold
createJSON("../../data/json_files/fastasplit_1_20.json", "../../data/fasta_files/Fasta_longest_chain_1_20")
createJSON("../../data/json_files/fastasplit_21_40.json", "../../data/fasta_files/Fasta_longest_chain_21_40")
createJSON("../../data/json_files/fastasplit_41_60.json", "../../data/fasta_files/Fasta_longest_chain_41_60")
createJSON("../../data/json_files/fastasplit_61_80.json", "../../data/fasta_files/Fasta_longest_chain_61_80")
createJSON("../../data/json_files/fastasplit_81_90.json", "../../data/fasta_files/Fasta_longest_chain_81_90")
