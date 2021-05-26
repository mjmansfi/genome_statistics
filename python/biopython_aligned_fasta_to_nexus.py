#!/bin/python3

import sys
import argparse
from argparse import RawTextHelpFormatter
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

parser = argparse.ArgumentParser(description="This script re-formats an aligned FASTA file, (the output for many multiple alignment programs) and returns the same alignment in NEXUS format.", formatter_class=RawTextHelpFormatter)
parser.add_argument("-f", "--fasta", help="Input FASTA file to be converted.\nDefault: none", required=True)
parser.add_argument("-m", "--molecule", help="Input molecule type/sequence alphabet.\nAvailable options: DNA, RNA, protein \nDefault: DNA", required=False, default="DNA")
parser.add_argument("-n", "--nexus", help="Output nexus file name.\nDefault: [input fasta].nexus", required=False)

args = parser.parse_args()

if args.fasta:
	fasta_file = args.fasta

if args.nexus:
	nexus_file = args.nexus
else:
	nexus_file = args.fasta + ".nexus"

molecule_types = ["DNA", "RNA", "protein"]
if args.molecule:
	molecule_type = args.molecule
	if molecule_type not in molecule_types:
		raise ValueError("Unexpected molecule type: %r" % args.alphabet)
	else:
		molecule_type = args.molecule

SeqIO.convert(fasta_file, "fasta", nexus_file, "nexus", molecule_type)

