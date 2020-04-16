#!/usr/bin/env python

from Bio import SeqIO
from Bio import Seq
import pandas as pd
import numpy as np
from glob import glob
import csv
import re
import os
import sys

# sys inputs
try:
    path = sys.argv[1]
    path_to_raw = sys.argv[2]
    length = sys.argv[3]
    motif = sys.argv[4]
    motif_path = sys.argv[5]
    output_path = sys.argv[6]

    if os.isdir(path):
        if not os.listdir(path)[0].endswith(".csv"):
            raise Exception("Files in directory are not CSV's.")
    else:
        raise Exception("Input is not a directory")

    if os.isdir(path_to_raw):
        if not os.listdir(path_to_raw)[0].endswith(".fasta"):
            raise Exception("Files in directory are not fasta's.")
    else:
        raise Exception("path_to_raw input is not a directory")

    if not type(length) == int:
        raise Exception("Length provided is not an integer.")

    if not type(motif) == str:
        raise Exception("Motif is not a string.")

    if not motif_path.endswith(".jaspar"):
        raise Exception("Motif file is not a jaspar file.")

    if not os.isdir(output_path):
        print("Creating output directory {0}".format(output_path))
        try:
            os.makedirs(output_path)
        except:
            pass
except:
    print("Incorrect number of arguments provided. Make sure to provide: \n1. {0} \n2. {1} \n3. {2} \n4. {3} \n5. {4} \n6. {5}".format("Path to the directory containing nonthresholded csv files", "Path to directory containing raw data", "Integer length of the sequence you want pulled before and after the motif", "String of motif name (e.g. \"bcd\")", "Jaspar file containing motif", "Output directory"))




def seqToString(motif):
    """
    motif: a Seq that represents the motif 
    
    Returns the String representation of the motif
    """
    i = 0
    string = ""
    length = motif.__len__()
    while i < length:
        string += motif.__getitem__(i)
        i += 1
    return string

def getNegative(pos_seq):
    """
    pos_seq = dna sequence in the positive direction reading from the file
    
    Returns the negative counterpart of the positive sequence.
    """
    dict = {"A":'T','T':'A','G':'C','C':'G','-':'-','N':'N'}
    negative = ""
    last_index = len(pos_seq) - 1
    while last_index > -1:
        negative += dict[pos_seq[last_index].upper()]
        last_index -= 1
    return negative

def getMotifLength(motif_key, path_to_motif):
    try:
        f = open(path_to_motif, 'r')
        length = 0
        for line in f:
            a = line.split()
            length = len(a)
            return length
    except:
        raise ValueError("The motif file {0} does not exist.".format(path_to_motif))

def prep_table(no_thresh_path):
    """ 
        path: The path to the folder with the thresholded and non thresholded directories
        
        Returns table containing the columns: strand direction, aligned position, score, species id, raw postion,
        and whether or not the position is considered the beginning of a motif
        
    """
    no_thresh = pd.read_csv(no_thresh_path) 
    thresh = no_thresh.loc[no_thresh["score"] >= 7] 
    verif_motifs = thresh.drop(columns = ['score', 'motif', 'raw_position', 'Unnamed: 0', 'strand']) 
    
    thresh = thresh.drop_duplicates(subset = 'align_position')
    thresh = thresh.drop(columns = ['score', 'motif', 'raw_position', 'Unnamed: 0'])
    no_thresh = no_thresh.drop(columns = ['motif', 'Unnamed: 0']) 
    #no_thresh = no_thresh.drop_duplicates(subset = ['align_position', 'strand'])
    result = pd.merge(thresh, no_thresh, on = ['align_position', 'strand'])

    result = result.drop(columns=['species_x'])
    return result 

def get_raw_string(region_group, length, no_thresh_path, motif_key, path_to_raw, motif_path, output_path):
    """ region_group: The regex group of the region number ex: 11048
        length: The length of the extra sequence that we want
        no_thresh_path: The path to the folder with the non thresholded files
        motif_key: the motif we are trying to extract (ex:bcd)
        path_to_raw: the path to the raw data
        motif_path: the path to the motif jaspar file
        output_path: path to output directory where csv files will be saved
        
        Returns a csv saved in the location specified in output_path
        """
    motif_length = getMotifLength(motif_key, motif_path)
    result = prep_table(no_thresh_path)
    num_region = region_group[1]
    region = region_group[0]
    if not os.path.isfile(os.path.join(output_path, motif_key, "{0}_final_raw.fa.csv".format(region))):
        # print(path_to_raw)
        # print("num_region:" + num_region)
        # print("region:" + region)
        region_path = os.path.join(path_to_raw, "outlier_rm_with_length_{0}.fa".format(region))
        record_dict = SeqIO.to_dict(SeqIO.parse(region_path, "fasta"))
        
        sequences = []
        before = []
        after = []
        
        for index, row in result.iterrows():
            speci = row['species_y']
            pos = row['raw_position']
            strand = row['strand']
            seq = record_dict[speci]
            if strand == 'negative':
                sequences.append(getNegative(seqToString(seq[pos:pos + motif_length])))
                before.append(getNegative(seqToString(seq[pos - length:pos])))
                after.append(getNegative(seqToString(seq[pos + motif_length:pos + motif_length + length])))
            else:
                sequences.append(seqToString(seq[pos:pos + motif_length]))
                before.append(seqToString(seq[pos - length:pos]))
                after.append(seqToString(seq[pos + motif_length:pos + motif_length + length]))

        result['raw_seq'] = np.array(sequences)
        result['before_seq'] = np.array(before)
        result['after_seq'] = np.array(after)
        result
        if not os.path.exists(output_path):
            os.mkdir(output_path)
        if not os.path.exists(os.path.join(output_path, motif_key)):
            os.mkdir(os.path.join(output_path, motif_key))
        result.to_csv(os.path.join(output_path, motif_key, "{0}_final_raw.fa.csv".format(region)))


def extract_motifs(path, path_to_raw, length, motif_key, motif_path, output_path):
    num_files = len(os.listdir(path))
    step = num_files // 10
    count = 0
    print("Beginning motif extraction for motif {0}.".format(motif_key))
    print("There are {0} files available in this directory.".format(num_files))
    for motif_file in os.listdir(path):
        count += 1
        region_group = re.search("VT([0-9]*)", motif_file)
        no_thresh_path = os.path.join(path, motif_file)
        get_raw_string(region_group, length, no_thresh_path, motif_key, path_to_raw, motif_path, output_path)
        if count % step == 0:
            print("Motif Extraction {0}% completed!".format(round(count/num_files*100, 2)))
    print("Motif Extraction Completed!")

extract_motifs(path, path_to_raw, length, motif, motif_path, output_path)
        






 
