import pandas as pd


def csv_to_latex_table(csv_file, output_file='latex_table.tex'):
    # Load the CSV data into a pandas DataFrame
    df = pd.read_csv(csv_file)

    # Pivot the data to get one row per loop length, with alpha and rosetta values in separate columns
    pivot_df = df.pivot(index='loop_len', columns='model', values=['mean_rmsd', 'sem_rmsd'])

    # Create the LaTeX table as a string
    latex_table = r"""\begin{table}[ht]
\centering
\begin{tabular}{ccc}
\toprule
\textbf{Loop Length} & \textbf{Mean RMSD Alpha (± SEM)} & \textbf{Mean RMSD Rosetta (± SEM)} \\
\midrule
"""

    # Iterate over each row to format the LaTeX table rows
    for loop_len in pivot_df.index:
        mean_alpha = pivot_df.loc[loop_len, ('mean_rmsd', 'alpha')]
        sem_alpha = pivot_df.loc[loop_len, ('sem_rmsd', 'alpha')]
        mean_rosetta = pivot_df.loc[loop_len, ('mean_rmsd', 'rosetta')]
        sem_rosetta = pivot_df.loc[loop_len, ('sem_rmsd', 'rosetta')]

        # Add a row to the LaTeX table string
        latex_table += f"{loop_len} & {mean_alpha:.2f} ± {sem_alpha:.2f} & {mean_rosetta:.2f} ± {sem_rosetta:.2f} \\\\ \n"

    # Add the bottom rule and closing LaTeX code
    latex_table += r"""\bottomrule
\end{tabular}
\caption{Mean RMSD with Standard Error (SEM) for Alpha and Rosetta Models}
\label{table:rmsd_alpha_rosetta}
\end{table}
"""

    # Write the LaTeX table to a .tex file
    with open(output_file, 'w') as f:
        f.write(latex_table)

    print(f"LaTeX table written to {output_file}")


csv_to_latex_table('mean_rmsd_by_loop_and_model.csv')
