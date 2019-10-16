# motif_extraction

The purpose of this project is to continue the pipeline that begins with checking for the existence of a transcription factor binding site, using tools from https://github.com/DiscoveryDNA/TFBS_presence. This project takes nonthresholded files outputted from this script and retrives the DNA sequence for the region of interest. 

The program:

1. Reads in non-thresholded files and thresholds them based on a score which indicates the likelihood of there being a motif.
2. Uses motif start index indicated in thresholded file to find the DNA sequences from fasta files containing the raw DNA sequences.
3. Populates a csv file with the DNA sequences for every species starting at the previously found index.


# Motivation

We want to be able to recognize patterns in the evolution of motif sequences, and analyze TFBS spacing patterns.


# Features

I'm unsure what to put here.
