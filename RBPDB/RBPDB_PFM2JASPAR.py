#!/usr/bin/python

__author__ = "bahin"
""" Script to extract data from the matrix list file of RBPDB bulk download (http://rbpdb.ccbr.utoronto.ca/downloads/PFMDir.zip) and convert it to JASPAR format.
Main features to deal with:
    - PWM matrices are stored in individual files
    - there might be several motifs for one RBP
"""

import argparse
import os
import pandas as pd

# Getting command-line options back
parser = argparse.ArgumentParser()
parser.add_argument("-r", dest="reference", required=True, help="Reference file containing the link between RBPs and PWM IDs.")
parser.add_argument("-p", dest="dirpath", required=True, help="Path to directory where the PWM matrices files are stored.")
parser.add_argument("-o", dest="output", required=True)
options = parser.parse_args()

# Parsing the reference file
RBPs = {}
with open(options.reference, "r") as input_file:
    for line in input_file:
        id = line.split("\t")[0]
        prot = line.split("\t")[2]
        # Listing the matrix IDs linked to RBPs
        if prot not in RBPs:
            RBPs[prot] = [id]
        elif id not in RBPs[prot]:
            RBPs[prot].append(id)

# Writing the output file
with open(options.output, "w") as output_file:
    for prot in RBPs:
        letter_nb = ord("A")  # If there are several motifs for one RBP, its IDs will followed by an underscore and an increasing letter
        for id in RBPs[prot]:
            # Setting the matrix file path
            prot_exp_filepath = options.dirpath + id + ".correct.pfm"
            # Checking if the matrix file exists
            if os.path.isfile(prot_exp_filepath):
                # Checking whether there is one or several motifs for the processed protein
                if len(RBPs[prot]) > 1:
                    output_file.write(">" + prot + "_" + chr(letter_nb) + "\n")
                    letter_nb += 1
                else:
                    output_file.write(">" + prot + "\n")
                # Getting the pfm file into a pandas DataFrame
                pfm_matrix = pd.read_csv(prot_exp_filepath, header=None)
                # Setting the indexes
                pfm_matrix.index = ["A", "C", "G", "T"]
                # Checking if this is a frequency or a counts matrix (round because sometime the sum equals 0.99999)
                if round(pfm_matrix.ix[:,0].sum(axis=0)) == 1.0:
                    pfm_matrix = pfm_matrix.multiply(100)
                # Writing the correctly formatted matrix to the output file
                for nt in pfm_matrix.index:
                    output_file.write(nt + " [ " + " ".join(
                        ([str(pfm_matrix.loc[nt].loc[c]) for c in pfm_matrix.columns])) + " ]\n")
            else:
                print "Warning: The file '" + prot_exp_filepath + "' doesn't exist. Aborting."