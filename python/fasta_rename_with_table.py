#!/bin/python3

import argparse
from argparse import RawTextHelpFormatter
from tqdm import tqdm
import os

def check_positive_integer(value):
    ivalue = int(value)
    if ivalue <= 0:
        raise argparse.ArgumentTypeError("%s is not 0 or a positive integer" % value)
    return ivalue

parser = argparse.ArgumentParser(description='Rename the headers of a fasta file to make them play nice with different programs. Writes the new headers to an easily parsed tsv file for when you need this information back.', formatter_class=RawTextHelpFormatter)
parser.add_argument('-f','--fasta', help='FASTA file to be renamed.', required=True)
parser.add_argument('-p','--prefix', help='Which prefix to prepend to renamed FASTA sequences.\nDefault: s (each sequence becomes s1, s2, ... sN)', default='s', required=False)
parser.add_argument('-o','--outfasta', help='Renamed FASTA output file.\nDefault: [input fasta].renamed.fa', required=False)
parser.add_argument('-tsv','--outtsv', help='Renamed TSV output file.\nDefault: [input fasta].renamed.tsv', required=False)
parser.add_argument('-seq', '--sequence', help='Boolean. Include sequence information in output .tsv?', required=False, action='store_true')
parser.add_argument('-std', '--standard', help='Remove sequences containing non-standard characters.\nWrites non-standard sequences to .warnings.txt.\nOptions:\n    [aa] - amino acids including X\n    [aax] - amino acids excluding X\n    [nt] - nucleotides including N\n    [ntn] - nucleotides excluding N', required=False, choices=['aa','aax','nt','ntn'])
args = parser.parse_args()

if args.outfasta:
	out_fasta = args.outfasta
else:
	out_fasta = os.path.basename(args.fasta) + '.renamed.fa'

prefix = args.prefix

if args.outtsv:
	out_tsv = args.outtsv
else:
	out_tsv = os.path.basename(args.fasta) + '.renamed.tsv'
	
if args.sequence:
	include_sequence = True

check_seqs_for_standard = False
if args.standard:
	check_seqs_for_standard = True
	alphabet = args.standard
	out_warnings = os.path.basename(args.fasta) + '.warnings.txt'

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

fasta_d = {}
with open(args.fasta) as i:
	print('Reading in fasta file %s...' % args.fasta)
	for line in tqdm(i):
		line = line.strip()
		if line.startswith('>'):
			header = line.split('>')[1]
			fasta_d[header] = ''
		else:
			fasta_d[header] += line

if check_seqs_for_standard:
	#make list of suspicious sequences
	pop_list = []
	print('Checking fasta sequences for non-standard characters...')
	with open(out_warnings, 'w') as warnings_handler:
		warnings_handler.write('fasta_header\tnonstandard_seq\n')
		for seq in tqdm(fasta_d):
			seq_status = check_standard(fasta_d[seq], alphabet)
			if not seq_status:
				#write the suspicious fasta to warnings.txt
				warnings_handler.write('%s\t%s\n' % (seq, fasta_d[seq]))
				pop_list.append(seq)
	#pop all suspicious sequences from fasta dict
	for seq in pop_list:
		fasta_d.pop(seq)

with open(out_fasta, 'w') as fasta_handler:
	with open(out_tsv, 'w') as tsv_handler:
		if args.sequence:
			tsv_header = 'renamed\toriginal\tseqlen\tsequence\n'
		else:
			tsv_header = 'renamed\toriginal\tseqlen\n'
		tsv_handler.write(tsv_header)
		print('Writing outputs to %s and %s...' % (out_fasta, out_tsv))
		seq_count = 0
		for seq in tqdm(fasta_d.keys()):
			seq_count += 1
			seq_len = len(fasta_d[seq])
			orig_header = seq
			seq_string = fasta_d[seq]
			fasta_handler.write('>%s%i\n%s\n' % (prefix, seq_count, seq_string))
			if args.sequence:
				tsv_out = ('%s%i\t%s\t%i\t%s\n' % (prefix, seq_count, orig_header, seq_len, seq_string))
			else:
				tsv_out = ('%s%i\t%s\t%i\n' % (prefix, seq_count, orig_header, seq_len))
			tsv_handler.write(tsv_out)
print("Done!")
print("    FASTA written to: %s\n    TSV written to: %s" % (out_fasta, out_tsv) )
if check_seqs_for_standard:
	print('    Warnings written to: %s' % (out_warnings))
