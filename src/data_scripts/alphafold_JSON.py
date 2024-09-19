"""Generate a JSON file from a folder of FASTA files for batch uploading to the AlphaFold3 webserver."""
import os
import json
from Bio import SeqIO


def create_JSON(folder_path: str) -> None:
    """
    Generate JSON files from a folder of FASTA files, with each JSON file containing up to 20 FASTA entries.

    :param folder_path: The directory containing the FASTA files to process.
    :returns: None. Writes the generated JSON content to separate files.
    :raises FileNotFoundError: If the specified folder path is invalid.
    :raises ValueError: If the FASTA file cannot be read or is in an incorrect format.
    """

    # List for collecting the entries
    list_of_entries = []

    # Get all FASTA files in the folder
    fasta_files = [f for f in os.listdir(folder_path) if f.endswith(".fasta")]
    fasta_files.sort()

    # Process each file
    for file_index, file_name in enumerate(fasta_files):
        file_path = os.path.join(folder_path, file_name)
        print(f"Processing {file_name}...")

        # Read Sequence from the FASTA file
        record = SeqIO.read(file_path, "fasta")

        # Get protein name (first 4 characters of the file name)
        name = str(file_name)[:4]
        sequence = str(record.seq)

        # Structure of the JSON entry
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

        # Write JSON file for each batch of 20
        if (file_index + 1) % 20 == 0 or file_index == len(fasta_files) - 1:
            output_json = os.path.join(folder_path, f"../../data/json_files/fastasplit_{file_index - 18}_{file_index + 1}.json")

            # Convert the list of entries to JSON format
            json_object = json.dumps(list_of_entries, indent=4)

            # Write the JSON object to the output file
            with open(output_json, 'w') as j_out:
                j_out.write(json_object)

            # Clear the list for the next batch
            list_of_entries = []


# Usage
create_JSON("../../data/fasta_files")

