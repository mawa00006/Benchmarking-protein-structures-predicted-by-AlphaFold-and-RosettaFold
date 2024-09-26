import os
import pickle
import matplotlib.pyplot as plt
import pandas as pd
from tueplots import bundles
import numpy as np
import scipy.stats as stats
import seaborn as sns

plt.rcParams.update(bundles.iclr2024(nrows=1, ncols=2))

df = pd.read_csv("results.csv")

alpha_df = df[df["model"] == "alpha"]
rosetta_df = df[df["model"] == "rosetta"]

mean_alpha = np.median(alpha_df["rmsd"].values)
mean_ros = np.median(rosetta_df["rmsd"].values)

fig, (ax1, ax2) = plt.subplots(1, 2)

sns.histplot(data=df, x='rmsd', hue='model', element='step', stat='percent', common_norm=False, ax=ax1)
ax1.set_xlabel('RMSD')
ax1.set_ylabel('Percentage')
ax1.set_xticks(np.arange(0, 7.1, 1))
ax1.axvline(mean_ros, color='orange', linestyle='dashed', linewidth=1)
ax1.axvline(mean_alpha, color='blue', linestyle='dashed', linewidth=1)

#ax1.set_title('Histogram of RMSD Frequencies by Model (Percentage)')

grouped = df.groupby(['loop_len', 'model']).agg(mean_rmsd=('rmsd', 'mean'), sem_rmsd=('rmsd', 'sem')).reset_index()



# Plotting with error bars
for model in grouped['model'].unique():
    model_data = grouped[grouped['model'] == model]
    ax2.errorbar(model_data['loop_len'], model_data['mean_rmsd'], yerr=model_data['sem_rmsd'], label=model, fmt=':,', markersize=3, markerfacecolor="k", linewidth=1)

ax2.set_xlabel('Loop Length')
ax2.set_ylabel('Average RMSD')
#ax2.set_title('Average RMSD vs Loop Length with Standard Error Bars')
ax2.legend(title='Model')
plt.savefig("rmsd.pdf", dpi=300)
plt.show()
