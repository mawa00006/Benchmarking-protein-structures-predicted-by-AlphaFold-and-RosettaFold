import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages

# Load the CSV file into a DataFrame 
data = pd.read_csv('results_adjacentSSE.csv')

df = data[data["model"] == "alpha"]

# Create groupings based on secondary structure conditions
helix_helix = df[(
    (df["secondary_structure_before"].isin(["H", "G", "I"])) &
    (df["secondary_structure_after"].isin(["H", "G", "I"]))
)]

beta_beta = df[(
    (df["secondary_structure_before"].isin(["B", "E"])) &
    (df["secondary_structure_after"].isin(["B", "E"]))
)]

beta_end = df[(
    ((df["secondary_structure_before"].isin(["B", "E"])) &
     (df["secondary_structure_after"] == "end")) |
    ((df["secondary_structure_before"] == "end") &
     (df["secondary_structure_after"].isin(["B", "E"])))
)]

helix_end = df[(
    ((df["secondary_structure_before"].isin(["H", "G", "I"])) &
     (df["secondary_structure_after"] == "end")) |
    ((df["secondary_structure_before"] == "end") &
     (df["secondary_structure_after"].isin(["H", "G", "I"])))
)]

beta_helix = df[(
    ((df["secondary_structure_before"].isin(["H", "G", "I"])) &
     (df["secondary_structure_after"].isin(["B", "E"]))) |
    ((df["secondary_structure_before"].isin(["B", "E"])) &
     (df["secondary_structure_after"].isin(["H", "G", "I"])))
)]

# Define colors for each loop type
loop_colors = {
    'Sheet-Sheet': 'blue',
    'Sheet-End': 'orange',
    'Sheet-Helix': 'green',
    'Helix-Helix': 'red',
    'Helix-End': 'purple'
}

# Create a list of the loop types and their corresponding dataframes
loop_types = {
    'Sheet-Sheet': beta_beta,
    'Sheet-End': beta_end,
    'Sheet-Helix': beta_helix,
    'Helix-Helix': helix_helix,
    'Helix-End': helix_end,
}

# Create a PDF to save the plot
with PdfPages('fig_loop_fractions_relative_to_adjacentSSE.pdf') as pdf:
    
    # Plot settings
    plt.figure(figsize=(8, 6))
    
    # Prepare the histogram data
    hist_data = [group_df['loop_len'].dropna() for group_df in loop_types.values()]
    # Create a stacked histogram
    plt.hist(hist_data, bins=np.arange(0, 26) - 0.5, stacked=True, 
             color=[loop_colors[loop_type] for loop_type in loop_types.keys()],
             alpha=0.7)

    # Labels and title
    plt.xlabel('Loop Length')
    plt.ylabel('Frequency')

    # Add legend and grid
    plt.legend(title='Loop Type', labels=loop_types.keys())
    plt.grid(axis='y', linestyle='--', linewidth=0.7)

    # Adjust layout
    plt.tight_layout()

    # Save the plot to the PDF
    pdf.savefig()
    plt.close()  # Close the plot after saving
