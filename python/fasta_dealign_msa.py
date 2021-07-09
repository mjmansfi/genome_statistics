#!/bin/python3

import argparse
from argparse import RawTextHelpFormatter
import os

parser = argparse.ArgumentParser(description='Rename the headers of a fasta file to make them play nice with different programs. Writes the new headers to an easily parsed tsv file for when you need this information back.', formatter_class=RawTextHelpFormatter)
parser.add_argument('-f','--fasta', help='FASTA file to be dealigned.', required=True)
parser.add_argument('-o','--outfasta', help='Optional: dealigned FASTA output file.\nDefault: standard output.', required=False)
args = parser.parse_args()

if args.outfasta:
	out_fasta = args.outfasta
else:
	out_fasta = False

fasta_d = {}
with open(args.fasta) as i:
	for line in i:
		line = line.strip()
		if line.startswith('>'):
			header = line.split('>')[1]
			fasta_d[header] = ''
		else:
			fasta_d[header] += line

if not out_fasta:
	for fasta_header in fasta_d:
		print('>%s\n%s' % ( fasta_header, fasta_d[fasta_header].replace('-','')))
else:
	with open(out_fasta, 'w') as output_handler:
		for fasta_header in fasta_d:
			output_handler.write('>%s\n%s\n' % (fasta_header, fasta_d[fasta_header].replace('-','')))
