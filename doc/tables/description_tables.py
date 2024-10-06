""" Create a LaTeX table from the summary statistics of the results """
import pandas as pd

# Load the results from the CSV files
df = pd.read_csv('../../src/results/results.csv')
df_gdt = pd.read_csv('../../src/results/results_GDT_TS.csv')

# Filter the dataframes by model
alpha_df = df[df["Model"] == "AlphaFold"]
rosetta_df = df[df["Model"] == "RoseTTAFold"]
alpha_df_gdt = df_gdt[df_gdt["Model"] == "AlphaFold"]
rosetta_df_gdt = df_gdt[df_gdt["Model"] == "RoseTTAFold"]

# Get the summary statistics for each model
description_alpha = alpha_df.describe()
description_ros = rosetta_df.describe()

# Manually create a LaTeX table
def df_to_latex(df, caption=""):
    latex_str = "\\begin{table}[ht]\n\\centering\n"
    if caption:
        latex_str += f"\\caption{{{caption}}}\n"
    latex_str += "\\begin{tabular}{c" + "c" * len(df.columns) + "}\n"
    latex_str += "\\toprule\n"
    latex_str += " & \\textbf{Loop Length} & \\textbf{RMSD}"
    latex_str += "\\midrule\n"
    # Iterate over each row to format the LaTeX table rows
    for idx, row in df.iterrows():
        # Escape % sign if it exists in the row index
        row_name = idx.replace("%", r"$\%$")
        latex_str += "{} & ".format(row_name) + " & ".join([f"{val:.2f}" for val in row]) + " \\\\\n"
    latex_str += "\\bottomrule\n"
    latex_str += "\\end{tabular}\n""\\end{table}\n\n"
    return latex_str


# Write the LaTeX table to a .tex file
with open('description_tables.tex', 'w') as f:
    # Write AlphaFold summary
    f.write(df_to_latex(description_alpha, caption="Summary of AlphaFold predictions"))
    # Write RoseTTAFold summary
    f.write(df_to_latex(description_ros, caption="Summary of RoseTTAFold predictions"))
