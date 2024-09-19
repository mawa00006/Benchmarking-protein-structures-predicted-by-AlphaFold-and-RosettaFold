"""
This script creates a JSON file for the AlphaFold3 webserver for batch upload.
Input: FOldername of fasta files which you want to upload
Output: JSOn file in the right formatting for upload
"""

import json
import os
from Bio import SeqIO


def createJSON(output_json, folder_path):

    #List for collecting the entries
    listOfEntries = []

    #Going through all files
    for file_name in os.listdir(folder_path):
        if file_name.endswith(".fasta"):
            file_path = os.path.join(folder_path, file_name)
            print(f"Processing {file_name}...")

            #It is assumed that only one sequence is in each file
            record = SeqIO.read(file_path, "fasta") 

            name = str(record.id)[0:4]
            sequence = str(record.seq)

            #structure of the JSOn file.
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

#enter the file name and the folder with the fasta files
createJSON("fastasplit41_60.json", "fastasplit41-60")

