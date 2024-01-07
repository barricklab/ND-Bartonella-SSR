#!/usr/bin/env python3

import argparse
import pandas
import glob
import os.path
import re
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SimpleLocation
from Bio.SeqFeature import ExactPosition

parser = argparse.ArgumentParser(description="Survey homopolymer CL candidates across genomes")
parser.add_argument("-g", "--genbank", help="Input folder of genbanks", action="store", required=True)
parser.add_argument("-o", "--output", help="Output CSV file with gene descriptions", action="store", default="output.csv")

args = parser.parse_args()


full_list = []

genbank_input_files=glob.glob(os.path.join(args.genbank, "*.gbk"))
for g in genbank_input_files:
    p = re.compile("/(.+).gbk")
    genome_id = p.search(g).group(1)
    print(genome_id)
    this_list = []
    for record in SeqIO.parse(g, "genbank"):
        seq_id = record.id

        for feature in record.features:
           if feature.type == 'CDS':
                locus_tag = feature.qualifiers.get("locus_tag", [""])[0]
                gene = feature.qualifiers.get("gene", [""])[0]
                product = feature.qualifiers.get("product", [""])[0]
                bp_length = len(feature)

                new_entry = { 'locus_tag' : locus_tag, 'gene' : gene, 'product' : product, "bp_length" : bp_length}

                full_list.append(pandas.DataFrame(new_entry, index=[0]))

full_df = pandas.concat(full_list, ignore_index=True)
full_df.to_csv(args.output, index=False)