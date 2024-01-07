#!/usr/bin/env python3

import argparse

from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, required=True)
parser.add_argument("-o", "--output", type=str, required=True)
args = parser.parse_args()

file_name = args.input
all_entries = []

with open(file_name, 'r') as GBFile:

    GBcds = SeqIO.InsdcIO.GenBankCdsFeatureIterator(GBFile)

    for cds in GBcds:
        #print(cds)
        if cds.seq is not None:
            cds.id = cds.name
            cds.description = ''
            all_entries.append(cds)


SeqIO.write(all_entries, args.output, 'fasta')