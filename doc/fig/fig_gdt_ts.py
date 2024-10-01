import os
import pickle
import matplotlib.pyplot as plt
import pandas as pd
from tueplots import bundles
import numpy as np
import scipy.stats as stats
import seaborn as sns
plt.rcParams.update(bundles.iclr2024())

df = pd.read_csv("results_GDT_TS.csv")

fig, ax = plt.subplots()
sns.violinplot(data=df, x="Model", y="GDT_TS", hue="Model", split=True, inner="quart")

ax.set_axisbelow(True)
ax.yaxis.grid(color='gray', linestyle='dashed')
ax.set_ylim([0, 100])
plt.savefig("fig_gdt_ts.pdf", dpi=300)
plt.savefig("fig_gdt_ts.png", dpi=300)
plt.show()

grouped = df.groupby("Model")

for group in grouped:
    print(group[1].describe())