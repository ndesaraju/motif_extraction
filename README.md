# motif_extraction

The purpose of this project is to continue the pipeline that begins with checking for the existence of a transcription factor binding site, using tools from https://github.com/DiscoveryDNA/TFBS_presence. This project takes nonthresholded files outputted from this script and retrives the DNA sequence for the region of interest. 

The program:

1. Reads in non-thresholded files and thresholds them based on a score which indicates the likelihood of there being a motif.
2. Uses motif start index indicated in thresholded file to find the DNA sequences from fasta files containing the raw DNA sequences.
3. Populates a csv file with the DNA sequences for every species starting at the previously found index.


# Motivation

We want to be able to recognize patterns in the evolution of motif sequences, and analyze TFBS spacing patterns.


# Usage

Inputs: Path to the directory containing nonthresholded csv files, path to directory containing raw data, integer length of the sequence you want pulled before and after the motif, string of motif name (e.g. "bcd"), jaspar file containing motif, and output directory.

Outputs: A directory of CSV files containing the motif sequence and the specified number of nucleotides before and after. If it doesn't exist, the output directory will be created.

Running the command: In terminal you should be able to type - 
\npython motif_extraction.py "/data/jaspar_redo_2019_09_06/bcd", "/data/3.24_species_only", 6, "bcd", "/data/jaspar_fm/modified/MA0212.1_bcd.jaspar" "/python_script_output"
