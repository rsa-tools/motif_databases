#!/usr/bin/python

__author__ = "bahin"
""" Script to extract data from ATtRACT database bulk download (https://attract.cnic.es/download) and convert it to cb format.
Main feature to deal with:
    - some motif stem from mutated genes, they have to be discarded
    - the database the motif stem from has to be reported
    - they name the files PWM but they are PFM
    - motifs from different databases are mixed in the PFM matrices file
    - there might be several motifs for one RBP
    - some of the several motifs might be redundant
"""

import argparse
import pandas as pd
import copy

# Getting command-line options back
parser = argparse.ArgumentParser()
parser.add_argument("-r", dest="reference", required=True, help="Reference file containing the link between RBPs and PFM IDs.")
parser.add_argument("-p", dest="pfm", required=True, help="File containing the PFM matrices.")
parser.add_argument("-o", dest="output", required=True)
options = parser.parse_args()

# Defining databases short names
databases = {"PDB": "PDB", "C": "CISBP-RNA", "R": "RBPDB", "S": "Spliceaid-F", "AEDB": "ASD"}

# Parsing the references file
RBPs = {}
with open(options.reference, "r") as db_file:
    db_file.readline()  # Burning the header line
    for line in db_file:
        if line.split("\t")[2] == "no":  # Motif from mutated genes are discarded
            id = line.split("\t")[11]
            prot = line.split("\t")[0]
            database = line.split("\t")[7]
            # Listing the PFM IDs linked to RBPs
            if prot not in RBPs:
                RBPs[prot] = {id: database}
            elif id not in RBPs[prot]:
                RBPs[prot][id] = database

# Parsing the PFM matrices file
matrices = {}
l = 1
with open(options.pfm, "r") as pfm_file:
    for line in pfm_file:
        # Processing header lines
        if line.startswith(">"):
            id = line.rstrip().split("\t")[0].split(">")[1]
            n = int(line.rstrip().split("\t")[1])
            # Storing the n next lines as a pandas DataFrame
            mat = pd.read_csv(options.pfm, sep="\t", header=None, skiprows=l, nrows=n)
            matrices[id] = mat
        l += 1

# Removing redundancies (same matrices with different IDs for the same RBP)
RBPs_copy = copy.deepcopy(RBPs)
for prot in RBPs_copy:  # Looping on the copy dict (the original dict will have some entries removed)
    to_keep = []  # Storing the IDs to keep (because one redundancy is found twice, we only want to remove one of the 2 redundant matrices=)
    for id in RBPs_copy[prot]:
        to_keep.append(id)  # Keeping one of the matrices when a redundancy is found
        for id2 in RBPs_copy[prot]:
            if (id != id2) and (matrices[id].equals(matrices[id2])):
                if (id2 in RBPs[prot]) and (id2 not in to_keep):  # The entry might already have been removed
                    del RBPs[prot][id2]  # Removing from the original dict

# Writing the output file
with open(options.output, "w") as output_file:
    for prot in RBPs:
        # If there are several motifs for one RBP, its IDs will followed by an underscore and an increasing letter or a left letter and a right letter if more than 26 motifs
        letter_nb = ord("A")  # Right letter index
        prefix_nb = 0  # Left letter ("tens") index if more than 26 motifs
        # Processing the RBPs
        for p in RBPs[prot]:
            if len(RBPs[prot]) > 1:  # If more than one motif for the processed RBP
                if prefix_nb == 0:  # If the motif processed is before the 26th, only one letter is joined
                    output_file.write(">" + databases[RBPs[prot][p]] + "__" + p + " /name=" + prot + "_" + chr(letter_nb) + "\n")
                else:  # Otherwise a couple of letters are joined
                    output_file.write(">" + databases[RBPs[prot][p]] + "__" + p + " /name=" + prot + "_" + chr(prefix_letter_nb) + chr(letter_nb) + "\n")
                if letter_nb == 90:  # If the loop reached "Z" for the right letter, setting the left letter
                    prefix_letter_nb = ord("A") + prefix_nb
                    prefix_nb += 1
                    letter_nb = 64  # Re-initializing the right letter
                letter_nb += 1
            else:
                output_file.write(">" + databases[RBPs[prot][p]] + "__" + p + " /name=" + prot + "\n")
            matrices[p].to_csv(output_file, sep="\t", header=False, index=False)