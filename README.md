# Benchmarking-protein-structures-predicted-by-AlphaFold-and-RosettaFold

# Project Description
Based on the paper of [Stevens and He](https://pubmed.ncbi.nlm.nih.gov/35883541/), our central contribution is the benchmarking of
AlphaFold and RoseTTAFold structure predictions using structures not used in model training with respect to loop length
and adjacent secondary structure type. Specifically we compare the accuracy of loop region predictions between AlphaFold 
and RoseTTAFold using the Root Mean Squared Deviation (RMSD) and the Global Distance Test Total Score (GDT-TS).


This project is part of the Structure Bioinformatics SS24 Group Project @Eberhard Karls University Tuebingen by Prof. Dr. Oliver Kohlbacher.

# Results
We found AlphaFold to give superior structural predictions than RoseTTAFold.
![alt text](https://github.com/mawa00006/Benchmarking-protein-structures-predicted-by-AlphaFold-and-RosettaFold/blob/main/doc/fig/png/fig_3_rmsd.png?raw=true)


### Reproducing Results
To reproduce all results (figures/tables/numbers) mentioned in our report follow these steps:
1. Create Python environment and install required packages 
```
# Setup new conda environment
conda create -n BenchmarkingProteinStructures

# Activate environment
conda activate BenchmarkingProteinStructures

# Install python
conda install python==3.9.19

# Install packages
conda install --file requirements.txt
```
If you also intend to redo the DSSP analysis you need to manually install DSSP from https://github.com/PDB-REDO/dssp

2. To reproduce the result csv files run the scripts located in src/evaluation_scripts
3. To reproduce figures run the scripts in doc/fig
4. To reproduce tables run the scripts in doc/tables

# Credits
This repository was created by **[Mattes Warning](https://github.com/mawa00006)**.

Additional contributors include:

- **[Florian Vögele](https://github.com/FloAvis)**
- **[Thusheera Kumarasekara](https://github.com/tckumarasekara)**
- **[Vishwa Patel](https://github.com/VPatel1412)**

Thank you to all the contributors for their work and dedication to the project. Special thanks to our Tutors Tom Mueller and Alexander Röhl.