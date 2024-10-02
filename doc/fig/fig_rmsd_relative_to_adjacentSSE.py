import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def create_plot(df1, df2, df3, df4, df5, model):
    # Create the plot
    plt.figure(figsize=(10, 6))
    
    ## Plot each dataset with error bars
    #plt.errorbar(df1['loop_len'], df1['mean_rmsd'], yerr=df1['sem_rmsd'], label='Sheet-Sheet', fmt='-o')
    #plt.errorbar(df2['loop_len'], df2['mean_rmsd'], yerr=df2['sem_rmsd'], label='Sheet-End', fmt='-o')
    #plt.errorbar(df3['loop_len'], df3['mean_rmsd'], yerr=df3['sem_rmsd'], label='Sheet-Helix', fmt='-o')
    #plt.errorbar(df4['loop_len'], df4['mean_rmsd'], yerr=df4['sem_rmsd'], label='Helix-Helix', fmt='-o')
    #plt.errorbar(df5['loop_len'], df5['mean_rmsd'], yerr=df5['sem_rmsd'], label='Helix-End', fmt='-o')
    
    plt.errorbar(df1['loop_len'], df1['mean_rmsd'], label='Sheet-Sheet', fmt='-o', linestyle='-', marker='o', linewidth=3)
    plt.errorbar(df2['loop_len'], df2['mean_rmsd'], label='Sheet-End', fmt='--^', linestyle='--', marker='^', linewidth=3)
    plt.errorbar(df3['loop_len'], df3['mean_rmsd'], label='Sheet-Helix', fmt='-.s', linestyle='-.', marker='s', linewidth=3)
    plt.errorbar(df4['loop_len'], df4['mean_rmsd'], label='Helix-Helix', fmt=':D', linestyle=':', marker='D', linewidth=3)
    plt.errorbar(df5['loop_len'], df5['mean_rmsd'], label='Helix-End', fmt='-.x', linestyle='-', marker='x', linewidth=3)
    
    # Labels and title
    plt.xlabel('Loop Length', fontsize=32)
    plt.ylabel('Average RMSD', fontsize=32)
    #plt.title(f'Average RMSD vs Loop Length Relative to Adjacent Secondary Structures for {model}')
    
    plt.xticks(fontsize=26)
    plt.yticks(fontsize=26)

    # Add legend
    plt.legend(title='Secondary Structure',
           fontsize=24, title_fontsize=28, loc='upper right',
           frameon=True, #bbox_to_anchor=(1.25, 1.25),
           handlelength=1.5, handletextpad=0.1, labelspacing=0.1)
    
    # Set outer border thickness
    plt.gca().spines['top'].set_linewidth(3)
    plt.gca().spines['right'].set_linewidth(3)
    plt.gca().spines['left'].set_linewidth(3)
    plt.gca().spines['bottom'].set_linewidth(3)
    
    
# Load the CSV file into a DataFrame 
df = pd.read_csv('results_adjacentSSE.csv')

alpha_df = df[df["model"] == "alpha"]
rosetta_df = df[df["model"] == "rosetta"]

alpha_df_helix_helix = alpha_df[((alpha_df["secondary_structure_before"] == "H") | 
                                 (alpha_df["secondary_structure_before"] == "G") |
                                 (alpha_df["secondary_structure_before"] == "I")) &
                                ((alpha_df["secondary_structure_after"] == "H") | 
                                 (alpha_df["secondary_structure_after"] == "G") |
                                 (alpha_df["secondary_structure_after"] == "I"))]
alpha_df_beta_beta = alpha_df[((alpha_df["secondary_structure_before"] == "B") | 
                               (alpha_df["secondary_structure_before"] == "E")) &
                              ((alpha_df["secondary_structure_after"] == "B") | 
                               (alpha_df["secondary_structure_after"] == "E"))]
alpha_df_beta_end = alpha_df[(((alpha_df["secondary_structure_before"] == "B") | 
                               (alpha_df["secondary_structure_before"] == "E")) &
                              (alpha_df["secondary_structure_after"] == "end")) |
                             ((alpha_df["secondary_structure_before"] == "end") &
                              ((alpha_df["secondary_structure_after"] == "B") | 
                               (alpha_df["secondary_structure_after"] == "E")))]
alpha_df_helix_end = alpha_df[(((alpha_df["secondary_structure_before"] == "H") | 
                                (alpha_df["secondary_structure_before"] == "G") |
                                (alpha_df["secondary_structure_before"] == "I")) &
                               (alpha_df["secondary_structure_after"] == "end")) |
                              ((alpha_df["secondary_structure_before"] == "end") &
                               ((alpha_df["secondary_structure_after"] == "H") | 
                                (alpha_df["secondary_structure_after"] == "G") |
                                (alpha_df["secondary_structure_after"] == "I")))]
alpha_df_beta_helix = alpha_df[(((alpha_df["secondary_structure_before"] == "H") | 
                                (alpha_df["secondary_structure_before"] == "G") |
                                (alpha_df["secondary_structure_before"] == "I")) &
                               ((alpha_df["secondary_structure_after"] == "B") | 
                                (alpha_df["secondary_structure_after"] == "E"))) |
                               (((alpha_df["secondary_structure_before"] == "B") | 
                                 (alpha_df["secondary_structure_before"] == "E")) &
                               ((alpha_df["secondary_structure_after"] == "H") | 
                                (alpha_df["secondary_structure_after"] == "G") |
                                (alpha_df["secondary_structure_after"] == "I")))]

rosetta_df_helix_helix = rosetta_df[((rosetta_df["secondary_structure_before"] == "H") | 
                                 (rosetta_df["secondary_structure_before"] == "G") |
                                 (rosetta_df["secondary_structure_before"] == "I")) &
                                ((rosetta_df["secondary_structure_after"] == "H") | 
                                 (rosetta_df["secondary_structure_after"] == "G") |
                                 (rosetta_df["secondary_structure_after"] == "I"))]
rosetta_df_beta_beta = rosetta_df[((rosetta_df["secondary_structure_before"] == "B") | 
                               (rosetta_df["secondary_structure_before"] == "E")) &
                              ((rosetta_df["secondary_structure_after"] == "B") | 
                               (rosetta_df["secondary_structure_after"] == "E"))]
rosetta_df_beta_end = rosetta_df[(((rosetta_df["secondary_structure_before"] == "B") | 
                               (rosetta_df["secondary_structure_before"] == "E")) &
                              (rosetta_df["secondary_structure_after"] == "end")) |
                             ((rosetta_df["secondary_structure_before"] == "end") &
                              ((rosetta_df["secondary_structure_after"] == "B") | 
                               (rosetta_df["secondary_structure_after"] == "E")))]
rosetta_df_helix_end = rosetta_df[(((rosetta_df["secondary_structure_before"] == "H") | 
                                (rosetta_df["secondary_structure_before"] == "G") |
                                (rosetta_df["secondary_structure_before"] == "I")) &
                               (rosetta_df["secondary_structure_after"] == "end")) |
                              ((rosetta_df["secondary_structure_before"] == "end") &
                               ((rosetta_df["secondary_structure_after"] == "H") | 
                                (rosetta_df["secondary_structure_after"] == "G") |
                                (rosetta_df["secondary_structure_after"] == "I")))]
rosetta_df_beta_helix = rosetta_df[(((rosetta_df["secondary_structure_before"] == "H") | 
                                (rosetta_df["secondary_structure_before"] == "G") |
                                (rosetta_df["secondary_structure_before"] == "I")) &
                               ((rosetta_df["secondary_structure_after"] == "B") | 
                                (rosetta_df["secondary_structure_after"] == "E"))) |
                               (((rosetta_df["secondary_structure_before"] == "B") | 
                                 (rosetta_df["secondary_structure_before"] == "E")) &
                               ((rosetta_df["secondary_structure_after"] == "H") | 
                                (rosetta_df["secondary_structure_after"] == "G") |
                                (rosetta_df["secondary_structure_after"] == "I")))]


# Average RMSD as a function of loop length with standard error bars for each model
grouped_alphafold_beta_beta = alpha_df_beta_beta.groupby(['loop_len']).agg(mean_rmsd=('rmsd', 'mean'), sem_rmsd=('rmsd', 'sem')).reset_index()
grouped_alphafold_beta_end = alpha_df_beta_end.groupby(['loop_len']).agg(mean_rmsd=('rmsd', 'mean'), sem_rmsd=('rmsd', 'sem')).reset_index()
grouped_alphafold_beta_helix = alpha_df_beta_helix.groupby(['loop_len']).agg(mean_rmsd=('rmsd', 'mean'), sem_rmsd=('rmsd', 'sem')).reset_index()
grouped_alphafold_helix_helix = alpha_df_helix_helix.groupby(['loop_len']).agg(mean_rmsd=('rmsd', 'mean'), sem_rmsd=('rmsd', 'sem')).reset_index()
grouped_alphafold_helix_end = alpha_df_helix_end.groupby(['loop_len']).agg(mean_rmsd=('rmsd', 'mean'), sem_rmsd=('rmsd', 'sem')).reset_index()

grouped_rosetta_beta_beta = rosetta_df_beta_beta.groupby(['loop_len']).agg(mean_rmsd=('rmsd', 'mean'), sem_rmsd=('rmsd', 'sem')).reset_index()
grouped_rosetta_beta_end = rosetta_df_beta_end.groupby(['loop_len']).agg(mean_rmsd=('rmsd', 'mean'), sem_rmsd=('rmsd', 'sem')).reset_index()
grouped_rosetta_beta_helix = rosetta_df_beta_helix.groupby(['loop_len']).agg(mean_rmsd=('rmsd', 'mean'), sem_rmsd=('rmsd', 'sem')).reset_index()
grouped_rosetta_helix_helix = rosetta_df_helix_helix.groupby(['loop_len']).agg(mean_rmsd=('rmsd', 'mean'), sem_rmsd=('rmsd', 'sem')).reset_index()
grouped_rosetta_helix_end = rosetta_df_helix_end.groupby(['loop_len']).agg(mean_rmsd=('rmsd', 'mean'), sem_rmsd=('rmsd', 'sem')).reset_index()

for grouped in [grouped_alphafold_beta_beta,grouped_alphafold_beta_end,grouped_alphafold_beta_helix,grouped_alphafold_helix_end,grouped_alphafold_helix_helix,
                grouped_rosetta_beta_beta,grouped_rosetta_beta_end,grouped_rosetta_beta_helix,grouped_rosetta_helix_end,grouped_rosetta_helix_helix]:
    grouped['mean_rmsd'] = grouped['mean_rmsd'].round(2)
    grouped['sem_rmsd'] = grouped['sem_rmsd'].round(2)
    

# Create a PDF to save the plots
with PdfPages('../src/evaluation_scripts/fig_rmsd_relative_to_adjacentSSE.pdf') as pdf:
       
    # First plot for model 'AlphaFold'
    create_plot(grouped_alphafold_beta_beta, grouped_alphafold_beta_end, grouped_alphafold_beta_helix,
                grouped_alphafold_helix_helix, grouped_alphafold_helix_end, model='AlphaFold')
    pdf.savefig()  # Save the current figure to the PDF
    plt.close()    # Close the figure to create a new one

    # Second plot for model 'Rosettafold' 
    create_plot(grouped_rosetta_beta_beta, grouped_rosetta_beta_end, grouped_rosetta_beta_helix,
                grouped_rosetta_helix_helix, grouped_rosetta_helix_end, model='RoseTTAFold')
    pdf.savefig()  # Save the second figure to the same PDF
    plt.close()    # Close the figure after saving

BB = len(alpha_df_beta_beta)*100/(len(alpha_df_beta_beta)+len(alpha_df_beta_end)+len(alpha_df_beta_helix)+len(alpha_df_helix_end)+len(alpha_df_helix_helix))       
BE = len(alpha_df_beta_end)*100/(len(alpha_df_beta_beta)+len(alpha_df_beta_end)+len(alpha_df_beta_helix)+len(alpha_df_helix_end)+len(alpha_df_helix_helix)) 
BH = len(alpha_df_beta_helix)*100/(len(alpha_df_beta_beta)+len(alpha_df_beta_end)+len(alpha_df_beta_helix)+len(alpha_df_helix_end)+len(alpha_df_helix_helix)) 
HH = len(alpha_df_helix_helix)*100/(len(alpha_df_beta_beta)+len(alpha_df_beta_end)+len(alpha_df_beta_helix)+len(alpha_df_helix_end)+len(alpha_df_helix_helix)) 
HE = len(alpha_df_helix_end)*100/(len(alpha_df_beta_beta)+len(alpha_df_beta_end)+len(alpha_df_beta_helix)+len(alpha_df_helix_end)+len(alpha_df_helix_helix)) 

print(f'Percentage of sheet-sheet {BB}')
print(f'Percentage of sheet-end {BE}')
print(f'Percentage of sheet-helix {BH}')
print(f'Percentage of helix-end {HE}')
print(f'Percentage of helix-helix {HH}')
