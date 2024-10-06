import os
import pickle
import matplotlib.pyplot as plt
import pandas as pd
from tueplots import bundles
import numpy as np
import scipy.stats as stats
import seaborn as sns

# Apply tueplots style
plt.rcParams.update(bundles.iclr2024(nrows=1, ncols=2))
palette = {"AlphaFold": (0, 0.39, 1), "RoseTTAFold": (1, 0.5, 0)}

# Load the CSV file into a DataFrame
df = pd.read_csv("../../src/results/results.csv")
# Filter the data to only include the AlphaFold and RoseTTAFold models
alpha_df = df[df["Model"] == "AlphaFold"]
rosetta_df = df[df["Model"] == "RoseTTAFold"]

# Calculate the mean RMSD for each model
mean_alpha = np.mean(alpha_df["RMSD"].values)
mean_ros = np.mean(rosetta_df["RMSD"].values)

# Create a figure with two subplots
fig, (ax1, ax2) = plt.subplots(1, 2)

# Plotting the histogram
sns.histplot(data=df, x='RMSD', hue='Model', element='step', stat='percent', common_norm=False, ax=ax1, palette=palette)
ax1.set_xlabel('RMSD')
ax1.set_ylabel('Percentage')
ax1.set_xticks(np.arange(0, 10.1, 1))
ax1.set_xlim([0,10])
ax1.axvline(mean_ros, color='orange', linestyle='dashed', linewidth=1)
ax1.axvline(mean_alpha, color='blue', linestyle='dashed', linewidth=1)

# Grouping by loop length and model
grouped = df.groupby(['loop_len', 'Model']).agg(mean_rmsd=('RMSD', 'mean'), sem_rmsd=('RMSD', 'sem')).reset_index()

# Plotting with error bars
for model in grouped['Model'].unique():
    model_data = grouped[grouped['Model'] == model]
    ax2.errorbar(model_data['loop_len'], model_data['mean_rmsd'], yerr=model_data['sem_rmsd'], label=model, fmt=':,', markersize=3, markerfacecolor="k", linewidth=1, color=palette[model])

# Add labels and title
ax2.set_xlabel('Loop Length')
ax2.set_ylabel('Average RMSD')
ax2.legend(title='Model')

# Save the figure
plt.savefig("pdf/fig_3_rmsd.pdf", dpi=300)
plt.savefig("png/fig_3_rmsd.png", dpi=300)
plt.show()
