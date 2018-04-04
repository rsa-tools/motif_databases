#!/usr/bin/python

__author__ = "bahin"
""" Script to extract data from CISBP-RNA bulk download (http://cisbp-rna.ccbr.utoronto.ca/bulk.php) and convert it to "JASPAR" format.
Main features to deal with:
    - PWM matrices are stored in individual files
    - there might be several motifs for one RBP
"""

import argparse
import os
import pandas as pd

# Getting command-line options back
parser = argparse.ArgumentParser()
parser.add_argument("-r", dest="reference", required=True, help="Reference file, 'RBP_Information_all_motifs.txt', containing the link between RBPs and PWM IDs.")
parser.add_argument("-o", dest="output", required=True, help="Unique JASPAR matrix.")
options = parser.parse_args()

# Parsing the reference file
RBPs = {}
with open(options.reference, "r") as input_file:
    input_file.readline()  # Burning the header line
    for line in input_file:
        filename = line.split("\t")[3]
        if filename != ".":  # Only processing lines linked to a PWM matrix file
            id = line.split("\t")[3]
            prot = line.split("\t")[6]
            # Listing the PWM matrix files linked to RBPs
            if prot not in RBPs:
                RBPs[prot] = [filename]
            elif filename not in RBPs[prot]:
                RBPs[prot].append(filename)

# Writing the output file
absent = []  # Listing absent matrix files (to print the warning only once)
with open(options.output, "w") as output_file:
    for prot in RBPs:
        letter_nb = ord("A")  # If there are several motifs for one RBP, its IDs will followed by an underscore and an increasing letter
        for filename in RBPs[prot]:
            # Setting the PWM matrix file path
            prot_exp_filepath = "pwms_counts_tf/" + filename + ".txt"
            # Checking if the matrix file exists
            if os.path.isfile(prot_exp_filepath):
                if len(RBPs[prot]) > 1:
                    output_file.write(">" + prot + "_" + chr(letter_nb) + " /name=" + filename + "\n")
                    letter_nb += 1
                else:
                    output_file.write(">" + prot + " /name=" + filename + "\n")
                pfm_matrix = pd.read_csv(prot_exp_filepath, sep="\t", header=None, index_col=0, skipfooter=1, engine="python")  # The engine is specified to avoid raising an error.
                # Writing the correctly formatted matrix to the output file
                for nt in ["A", "C", "G", "T"]:
                    output_file.write(
                        nt + " [ " + " ".join(([str(pfm_matrix.loc[nt].loc[c]) for c in pfm_matrix.columns])) + " ]\n")
            else:
                if prot_exp_filepath not in absent:
                    absent.append(prot_exp_filepath)
                    print "Warning: The file '" + prot_exp_filepath + "' doesn't exist."