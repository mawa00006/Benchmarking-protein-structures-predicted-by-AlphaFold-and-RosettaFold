""" This script demonstrates how to find adjacent entries in a DataFrame based on specific conditions.
    It was used to create the alignment figure in the report.
"""
import pandas as pd
import numpy as np

# Load the CSV file
file_path = '../src/results/results.csv'
data = pd.read_csv(file_path)


def assign_pair_number(df):
    df = df.reset_index(drop=True)
    df['loop_id'] = (df.index // 2) + 1  # Apply the pair numbering logic
    return df

# Apply the function to each group of structure_id
data = data.groupby('structure_id', group_keys=False).apply(assign_pair_number)

# Display the updated DataFrame with pair_number
print(data.head(20))

# Sort the data by structure_id and loop_len
data_sorted = data.sort_values(by=['structure_id', 'loop_len', 'loop_id'])


def find_adjacent_entries(df,loop_len, rosettaRMSD, alphaRMSD):
    adjacent_entries = ()
    bestDistance = np.infty

    for i in range(len(df) - 1):
        current_row = df.iloc[i]
        next_row = df.iloc[i + 1]

        # Check if structure_id and loop_len match, and models are different
        if (current_row['loop_id'] == next_row['loop_id'] and
            current_row['loop_len'] == loop_len):

            alphaDistance = abs(alphaRMSD - current_row["RMSD"])
            rosettaDistance = abs(rosettaRMSD - next_row["RMSD"])

            overallDistance = alphaDistance + rosettaDistance


            # Check if one row is 'alpha' with RMSD ~ 0.33, and the other is 'rosetta' with RMSD ~ 0.49
            if current_row['Model'] != 'AlphaFold' or next_row['Model'] != 'RoseTTAFold':
                print("ERROR!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
                      f"Alpha: {current_row['Model']}, Rosetta: {next_row['Model']}")

            if overallDistance < bestDistance:
                adjacent_entries = (current_row, next_row)
                bestDistance = overallDistance

    return adjacent_entries

loop_length = [3, 5, 7, 14, 24]
alphaRMSDs = [0.84, 0.89, 0.99, 0.68, 1.10]
rosettaRMSDs = [1.28, 1.32, 1.41, 1.20, 1.13]


for i in range(len(loop_length)):
    adjacent_entries = find_adjacent_entries(data_sorted,
                                             loop_len=loop_length[i],
                                             rosettaRMSD=rosettaRMSDs[i],
                                             alphaRMSD=alphaRMSDs[i])

    # Display the results
    print("Loop Length: ", loop_length[i])
    print("Structure ID:", adjacent_entries[0]['structure_id'])
    print("Loop ID:", adjacent_entries[0]['loop_id'])
    print("Alpha RMSD:",
          adjacent_entries[0]['RMSD'] if adjacent_entries[0]['Model'] == 'AlphaFold' else adjacent_entries[1]['RMSD'])
    print("Rosetta RMSD:",
          adjacent_entries[0]['RMSD'] if adjacent_entries[0]['Model'] == 'RoseTTAFold' else adjacent_entries[1]['RMSD'])
    print("-" * 40)