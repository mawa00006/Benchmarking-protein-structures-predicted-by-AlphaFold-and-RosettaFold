import os
import pickle
import matplotlib.pyplot as plt
import pandas as pd
from tueplots import bundles
import numpy as np
import scipy.stats as stats
import seaborn as sns
plt.rcParams.update(bundles.iclr2024(nrows=1, ncols=2))
palette = {"AlphaFold": (0, 0.39, 1), "RoseTTAFold": (1, 0.5, 0)}

df = pd.read_csv("../../src/results/results_GDT_TS.csv")

alpha_df = df[df["Model"] == "AlphaFold"]
rosetta_df = df[df["Model"] == "RoseTTAFold"]

fig, (ax1, ax2) = plt.subplots(1, 2)
sns.violinplot(data=df, x="Model", y="GDT_TS", hue="Model", split=True, inner="quart", ax=ax1, palette=palette)

ax1.set_axisbelow(True)
ax1.yaxis.grid(color='gray', linestyle='dashed')
ax1.set_ylim([0, 100])

differences = alpha_df["GDT_TS"].values - rosetta_df["GDT_TS"].values
cnts, values, bars = ax2.hist(differences, rwidth=1)

for i, (cnt, value, bar) in enumerate(zip(cnts, values, bars)):
    if value > 0:
        bar.set_facecolor(palette["AlphaFold"])
    else:
        bar.set_facecolor(palette["RoseTTAFold"])
ax2.set_axisbelow(True)
ax2.yaxis.grid(color='gray', linestyle='dashed')
ax2.set_xlabel("GDT_TS_Alpha - GDT_TS_RoseTTA")
ax2.set_ylabel("Count")

plt.savefig("pdf/fig_2_gdt_ts.pdf", dpi=300)
plt.savefig("png/fig_2_gdt_ts.png", dpi=300)
plt.show()