# Post Translational Modification of Proteins
Description
This repository provides R scripts for post-translational modification analysis in proteomics data, focusing on six amino acids: Arginine (R), Lysine (K), Proline (P), Threonine (T), Cysteine (C), and Methionine (M). However, modifications on Cysteine (C) and Methionine (M) are disregarded due to their alteration by chemicals used in sample processing. The analysis includes identifying peptide positions, summarizing modifications, and generating modification-to-total ratios for each amino acid.

Features

Identification of modified peptide positions within protein sequences.
Calculation of modification summaries for R, K, P, and T amino acids, with C and M modifications ignored.
Generation of modification-to-total ratios, offering insights into modification patterns across peptides.
Project Structure

Functions:

read_n_skip_fun: Reads and prepares CSV files by skipping non-relevant lines.

get_non_zero_positions: Identifies positions of non-zero modification indicators.

add_peptide_positions: Maps peptide sequences within protein sequences to locate modifications.

calculate_modification_summary: Summarizes counts of modifications for R, K, P, and T.

calculate_summary: Computes total amino acid counts for normalization.

raw_to_ratio: Executes the entire workflow, from data preparation to modification ratio calculation and CSV export.
Data Requirements

Input CSV files should include columns specifying peptides, protein sequences, variable modifications, and modification positions.

How to Use

Clone the repository and install dependencies (e.g., tidyverse).
Run raw_to_ratio(<folder_path>), where <folder_path> points to the folder containing input CSV files.
The output will be saved as CSV files in a "ratio" folder within the specified directory.
