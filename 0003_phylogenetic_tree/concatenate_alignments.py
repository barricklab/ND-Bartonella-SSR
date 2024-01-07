#!/usr/bin/env python3


import glob
import os.path
import re
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

concatenated_sequences = dict()
alignment_length = 0

input_files=glob.glob(os.path.join("03_single_copy_ortholog_protein_alignments_bmge_trimmed", "*.fa"))
for g in input_files:
    print(g)
    this_alignment_length = None
    found_in_this_alignment = dict()
    for record in SeqIO.parse(g, "fasta"):
        id = record.id.split("_")[0]

        if this_alignment_length is None:
            this_alignment_length = len(record)
        elif len(record) != this_alignment_length:
            sys.exit("Disagreement on alignment length.")

        #print(id)
        #add if not found
        if not id in concatenated_sequences:
            concatenated_sequences[id] = SeqRecord("-" * alignment_length, id=id)
        concatenated_sequences[id] = concatenated_sequences[id] + record.seq

        found_in_this_alignment[id] = True

    #add sequence for missing ones - but there aen't any!
    for id in concatenated_sequences.keys():
        if not id in found_in_this_alignment:
            print("Missing padded: " + id)
            concatenated_sequences[id] = concatenated_sequences[id] + "-" * this_alignment_length

    alignment_length = alignment_length + this_alignment_length

print("Concatenated alignment total sequences: " + str(len(concatenated_sequences)))
print("Concatenated alignment total length: " + str(alignment_length))

# write file of all sequences
output_file = "concatenated.fasta"
with open(output_file, "w") as f:
    SeqIO.write(concatenated_sequences.values(), f, "fasta")