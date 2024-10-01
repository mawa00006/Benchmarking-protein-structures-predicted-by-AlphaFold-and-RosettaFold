import os
import pickle
import matplotlib.pyplot as plt
import pandas as pd
from tueplots import bundles
import numpy as np
import scipy.stats as stats

plt.rcParams.update(bundles.iclr2024(nrows=1, ncols=2))

folder_path = '../../data/loop_annotations'
df = pd.read_csv("results.csv")

# Initialize an empty list to store the percentages
percentages = []

# Iterate over each pickle file in the folder
for file_name in os.listdir(folder_path):
    if file_name.endswith('.pickle'):  # Only process pickle files
        file_path = os.path.join(folder_path, file_name)

        # Load the pickle file
        with open(file_path, 'rb') as f:
            data = pickle.load(f)

        # Calculate the percentage of 1's in the loop annotation
        percentage_of_ones = sum(data) / len(data)
        percentages.append(percentage_of_ones)

percentages = np.array(percentages)
std = np.std(percentages)


fig, (ax1, ax2) = plt.subplots(1, 2)

mu = np.median(percentages)
variance = np.var(percentages)
sigma = np.std(percentages)
x = np.linspace(mu - 3*sigma, mu + 3*sigma, 100)
ax1.plot(x, stats.norm.pdf(x, mu, sigma))


num_bins = 20  # Define the number of bins
hist, bin_edges, _ = ax1.hist(percentages, bins=num_bins, range=(0, 1), edgecolor='black', color="lightgrey")


# Add labels and title
ax1.set_ylabel('Frequency')
ax1.set_xlabel('Fraction of loops in longest chain of each structure')

ax1.set_xlim(0, 1)
ax1.set_xticks(np.arange(0, 1.01, 0.1))
ax1.set_axisbelow(True)
ax1.yaxis.grid(color='gray', linestyle='dashed')

y_min, y_max = ax1.get_ylim()
unique_loops = df['loop_len'].value_counts() / 2

print(unique_loops.values)

# Plot histogram of loop lengths
ax2.bar(unique_loops.index, unique_loops.values, color='skyblue')
ax2.set_xlabel('Loop length')
ax2.set_ylabel('Frequency')

y_min, y_max = ax2.get_ylim()
ax2.text(6, y_max*0.85, 'Mean: 5')
ax2.axvline(5, color='k', linestyle='dashed', linewidth=1)
ax2.set_axisbelow(True)
ax2.yaxis.grid(color='gray', linestyle='dashed')

plt.savefig("fig_loop_fractions.pdf", dpi=300)
plt.savefig("fig_loop_fractions.png", dpi=300)
plt.show()