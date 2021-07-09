#!/bin/python3

"""

(C) Mike Mansfield, PhD 2021

This script selects N random positions from a FASTA file which are not
soft-masked. The reason this is needed is for NanoSV, which requires
a .bed file containing non-soft-masked positions in a BAM file for
depth support of variants.

"""

import binascii
import gzip
import os
import sys
import argparse
from argparse import RawTextHelpFormatter
import random

parser = argparse.ArgumentParser(description="This script selects N random positions from a FASTA file which are not soft-masked. The reason this is needed is for NanoSV, which requires a .bed file containing non-soft-masked positions in a BAM file for depth support of variants.", formatter_class=RawTextHelpFormatter)
parser.add_argument("-f", "--fasta", help="Input soft-masked FASTA file.\nDefault: none", required=True)
parser.add_argument("-N", "--Npositions", help="Number of positions to sample.\nDefault: 1000000", required=False, default=1000000, type=int)
parser.add_argument("-m", "--minContigLength", help="Minimum contig length for sampling.\nDefault: 0", required=False, default=0)
parser.add_argument("-o", "--output", help="Output bed file.\nDefault: [stdout]", required=False, default=False)

args = parser.parse_args()


"""
Simply checks if a file is GZIP compressed or not. Instead of hacking away with file
extensions, instead this checks the first two bytes of a file after opening. For all
GZIP files, the first two bytes of the file are "1f8b", which will not likely be the
case if it's not a compressed file.
"""
def is_gzip_file(filePath):
    with open(filePath, "rb") as test_file:
        return binascii.hexlify(test_file.read(2)) == b'1f8b'

"""
Simple FASTA reader. A FASTA parser can be written as a generator instead, but since
it passes through a few if/etc. logical separators, so this is easier to read.
"""
def read_fasta_file(fastaFile, minContigLength=0):
	fasta_dictionary = {}
	if is_gzip_file(fastaFile):
		with gzip.open(fastaFile, "rt") as fasta_handler:
			for line in fasta_handler:
				line = line.strip()
				if line.startswith(">"):
					header = line.split(">")[1]
					fasta_dictionary[header] = []
				else:
					fasta_dictionary[header].append(line)
	else:
		with open(fastaFile, "r") as fasta_handler:
			for line in fasta_handler:
				line = line.strip()
				if line.startswith(">"):
					header = line.split(">")[1]
					fasta_dictionary[header] = []
				else:
					fasta_dictionary[header].append(line)
	fasta_dictionary = {header : "".join(fasta_dictionary[header]) for header in fasta_dictionary if len(fasta_dictionary[header]) > minContigLength}
	return(fasta_dictionary)


"""
Select N random positions that are not soft-masked.

Note that the bed file is 1-indexed, not 0-indexed.
"""
def select_N_nonmasked_positions(fastaDict, N_positions=1000000):
	nonmasked_pos_list = []
	for fasta_header in fastaDict:
		sequence = fastaDict[fasta_header]
		sequence_length = len(sequence)
		for position in range(0, sequence_length):
			sequence_position = sequence[position]
			start = position + 1
			end = position + 2
			if sequence_position not in ['n', 'N', '-']:
				if not sequence_position.islower():
					nonmasked_pos_list.append([fasta_header, str(start), str(end)])

	random_positions = random.sample(nonmasked_pos_list, N_positions)
	
	return(random_positions)


"""
Print (or write) the bed file.
"""
def write_bed_file(nonmasked_pos_list, output=None):
	nonmasked_pos_list = sorted(nonmasked_pos_list)
	if output:
		with open(output, "w") as o:
			[o.write('\t'.join(pos) + '\n') for pos in nonmasked_pos_list]
	else:
		[print('\t'.join(pos)) for pos in nonmasked_pos_list]


"""
Wrapping up everything in one function.
"""
def workflow():
	fasta_dict = read_fasta_file(args.fasta)
	position_list = select_N_nonmasked_positions(fasta_dict, N_positions=args.Npositions)
	write_bed_file(position_list, output=args.output)

if __name__ == "__main__":
	workflow()
else:
	print("Why are you importing this? Don't.")
