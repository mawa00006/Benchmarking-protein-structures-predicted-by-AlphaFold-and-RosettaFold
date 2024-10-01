import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Load the CSV file into a DataFrame (assuming the file is named 'data.csv')
df = pd.read_csv('results.csv')
df_gdt = pd.read_csv('results_GDT_TS.csv')

alpha_df = df[df["Model"] == "AlphaFold"]
rosetta_df = df[df["Model"] == "RoseTTAFold"]
alpha_df_gdt = df_gdt[df_gdt["Model"] == "AlphaFold"]
rosetta_df_gdt = df_gdt[df_gdt["Model"] == "RoseTTAFold"]

print("alpha", alpha_df.describe())
print("rosetta", rosetta_df.describe())
print("95% loop length", alpha_df["loop_len"].quantile(0.95))
print("99% loop length", alpha_df["loop_len"].quantile(0.99))
print(np.median(alpha_df["RMSD"].values))
print(np.median(rosetta_df["RMSD"].values))

print("alpha", alpha_df_gdt.describe())
print("rosetta", rosetta_df_gdt.describe())


# 3. Average RMSD as a function of loop length with standard error bars for each model
# Group by loop_len and model, then calculate mean and standard error (sem)
grouped = df.groupby(['loop_len', 'Model']).agg(mean_rmsd=('RMSD', 'mean'), sem_rmsd=('RMSD', 'sem')).reset_index()


grouped['mean_rmsd'] = grouped['mean_rmsd'].round(2)
grouped['sem_rmsd'] = grouped['sem_rmsd'].round(2)

# Save the rounded data to a new CSV file
grouped.to_csv('mean_rmsd_by_loop_and_model.csv', index=False)


alpha_df = alpha_df[alpha_df["loop_len"] > 2]
rosetta_df = rosetta_df[rosetta_df["loop_len"] > 2]

print(np.median(alpha_df["RMSD"].values))
print(np.median(rosetta_df["RMSD"].values))

print("alpha", alpha_df.describe())
print("rosetta", rosetta_df.describe())

