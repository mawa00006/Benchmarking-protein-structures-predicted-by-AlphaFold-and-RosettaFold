import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Load the CSV file into a DataFrame (assuming the file is named 'data.csv')
df = pd.read_csv('results.csv')

alpha_df = df[df["model"] == "alpha"]
rosetta_df = df[df["model"] == "rosetta"]

print("alpha", alpha_df.describe())
print("rosetta", rosetta_df.describe())
print("95% loop length", alpha_df["loop_len"].quantile(0.95))
print("99% loop length", alpha_df["loop_len"].quantile(0.99))

# 1. Histogram of loop lengths (dividing loop count by 2)
# Since each loop has two entries, we will count unique loop lengths and divide by 2
unique_loops = df['loop_len'].value_counts() / 2

# Plot histogram of loop lengths
plt.figure(figsize=(8, 6))
plt.bar(unique_loops.index, unique_loops.values, color='skyblue')
plt.xlabel('Loop Length')
plt.ylabel('Frequency')
plt.title('Histogram of Loop Lengths')
plt.savefig("loop_hist.pdf", dpi=100)
plt.show()


# 2. Histogram of RMSD frequencies grouped by model
# We'll use seaborn's histplot for this with normalized frequency (percentage)
plt.figure(figsize=(8, 6))
sns.histplot(data=df, x='rmsd', hue='model', element='step', stat='percent', common_norm=False)
plt.xlabel('RMSD')
plt.ylabel('Percentage')
plt.title('Histogram of RMSD Frequencies by Model (Percentage)')
plt.savefig("rmsd_freq.pdf")
plt.show()


# 3. Average RMSD as a function of loop length with standard error bars for each model
# Group by loop_len and model, then calculate mean and standard error (sem)
grouped = df.groupby(['loop_len', 'model']).agg(mean_rmsd=('rmsd', 'mean'), sem_rmsd=('rmsd', 'sem')).reset_index()



# Plotting with error bars
plt.figure(figsize=(8, 6))
for model in grouped['model'].unique():
    model_data = grouped[grouped['model'] == model]
    plt.errorbar(model_data['loop_len'], model_data['mean_rmsd'], yerr=model_data['sem_rmsd'], label=model, fmt='-o')

plt.xlabel('Loop Length')
plt.ylabel('Average RMSD')
plt.title('Average RMSD vs Loop Length with Standard Error Bars')
plt.legend(title='Model')
plt.savefig("rmsd_by_loop_len.pdf")
plt.show()



grouped['mean_rmsd'] = grouped['mean_rmsd'].round(2)
grouped['sem_rmsd'] = grouped['sem_rmsd'].round(2)

# Save the rounded data to a new CSV file
grouped.to_csv('mean_rmsd_by_loop_and_model.csv', index=False)