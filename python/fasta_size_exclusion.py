#!/bin/python3

import argparse
from argparse import RawTextHelpFormatter
import os
import collections

def check_positive_integer(value):
    ivalue = int(value)
    if ivalue <= 0:
        raise argparse.ArgumentTypeError("%s is not 0 or a positive integer" % value)
    return ivalue

parser = argparse.ArgumentParser(description='Rename the headers of a fasta file to make them play nice with different programs. Writes the new headers to an easily parsed tsv file for when you need this information back.', formatter_class=RawTextHelpFormatter)
parser.add_argument('-f','--fasta', help='FASTA file to be renamed.\nRequired.', required=True)
parser.add_argument('-o', '--outfasta', help='Renamed FASTA output file.\nDefault: standard output.', required=False)
parser.add_argument('-m', '--minlen', help='Minimum length to be included in output file.\nDefault: none',  required=False, default = 0, type=check_positive_integer)
parser.add_argument('-std', '--standard', help='Remove sequences containing non-standard characters.\nWrites non-standard sequences to [input fasta].warnings.txt.\nOptions:\n    [aa] - amino acids including X\n    [aax] - amino acids excluding X\n    [nt] - nucleotides including N\n    [ntn] - nucleotides excluding N', required=False, choices=['aa','aax','nt','ntn'])
args = parser.parse_args()

if args.fasta:
	fasta = args.fasta

check_seqs_for_standard = False
if args.standard:
	check_seqs_for_standard = True
	alphabet = args.standard
	out_warnings = os.path.basename(args.fasta) + '.warnings.txt'

if args.minlen:
    min_contig_length = args.minlen
else:
    min_contig_length = 0

if args.outfasta:
	outfasta = args.outfasta
else:
	outfasta = False

def check_standard(sequence, alphabet):
	check = True
	if alphabet == 'aa':
		standard_aas = 'ACDEFGHIKLMNPQRSTVWXY'
		for letter in sequence:
			if letter.upper() not in standard_aas:
				check = False
	if alphabet == 'aax':
		standard_aas = 'ACDEFGHIKLMNPQRSTVWY'
		for letter in sequence:
			if letter.upper() not in standard_aas:
				check = False
	if alphabet == 'nt':
		standard_nts = 'ATCGN'
		for letter in sequence:
			if letter.upper() not in standard_nt:
				check = False
	if alphabet == 'ntn':
		standard_nts = 'ATCG'
		for letter in sequence:
			if letter.upper() not in standard_nt:
				check = False
	return check


"""
Simple FASTA reader. A FASTA parser can be written as a generator instead, but since
it passes through a few if/etc. logical separators, so this is easier to read.
"""
def read_fasta_file(fastaFile, minContigLength=0):
	fasta_dictionary = {}
	num_lines = sum(1 for line in open(fastaFile, "r"))
	with open(fastaFile, "r") as fasta_handler:
		for line in fasta_handler:
			line = line.strip()
			if line.startswith(">"):
				header = line.split(">")[1]
				fasta_dictionary[header] = []
			else:
				fasta_dictionary[header].append(line)
	# Trim the whitespace
	fasta_dictionary = {header : "".join(fasta_dictionary[header]) for header in fasta_dictionary}
	# Then exclude by length
	fasta_dictionary = {header : fasta_dictionary[header] for header in fasta_dictionary if len(fasta_dictionary[header]) > minContigLength}
	return(fasta_dictionary)

fasta_d = read_fasta_file(fasta, minContigLength=min_contig_length)

if check_seqs_for_standard:
	pop_list = []
	with open(out_warnings, 'w') as warnings_handler:
		warnings_handler.write('fasta_header\tnonstandard_seq\n')
		for seq in fasta_d:
			seq_status = check_standard(fasta_d[seq], alphabet)
			if not seq_status:
				#write the suspicious fasta to warnings.txt
				warnings_handler.write('%s\t%s\n' % (seq, fasta_d[seq]))
				pop_list.append(seq)
	#pop all suspicious sequences from fasta dict
	for seq in pop_list:
		fasta_d.pop(seq)

fasta_d = collections.OrderedDict(sorted(fasta_d.items(), key = lambda x : len(x[1]), reverse=True))

if outfasta:
	with open(outfasta, 'w') as out:
		for fasta_header in fasta_d:
			out.write(">%s\n%s\n" % (fasta_header, fasta_d[fasta_header]))
else:
	for fasta_header in fasta_d:
		print(">%s\n%s" % (fasta_header, fasta_d[fasta_header]))
