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
parser.add_argument('-p','--prefix', help='Which prefix to prepend to renamed FASTA sequences.\nDefault: None (each sequence becomes 1, 2, ... N)', default="", required=False)
parser.add_argument('-ap','--add-prefix', help='If True, then prepend option provided by --prefix to all names.\nDefault: False', default=False, required=False, action='store_true')
parser.add_argument('-k', '--keep-old-names', help='Keep old name of the sequence. Can still prepend a prefix if --prefix is provided and --add-prefix is used. Useful to filter only for nonstandard characters without renaming.\nDefault: False, will do the renaming.', default=False, required=False, action='store_true')
parser.add_argument('-o','--outfasta', help='Renamed FASTA output file.\nDefault: [input fasta].renamed.fasta', required=False)
parser.add_argument('-tsv','--outtsv', help='Renamed TSV output file.\nDefault: [input fasta].renamed.tsv', required=False)
parser.add_argument('-seq', '--sequence', help='Boolean. Include sequence information in output .tsv?', required=False, action='store_true')
parser.add_argument('-std', '--standard', help='Remove sequences containing non-standard characters.\nWrites non-standard sequences to .warnings.txt.\nOptions:\n    [aa] - 20 canonical amino acids \n    [aax] - canonical amino acids including X\n    [aastar] - canonical amino acids including "*"\n    [aaxstar] - amino acids including "*" and X\n    [nt] - 4 canonical nucleotides\n    [ntn] - canonical nucleotides including N', required=False, choices=['aa','aax','aastar','aaxstar','nt','ntn',])
args = parser.parse_args()

if args.outfasta:
	out_fasta = args.outfasta
else:
	out_fasta = os.path.basename(args.fasta) + '.renamed.fasta'



if args.outtsv:
	out_tsv = args.outtsv
else:
	out_tsv = os.path.basename(args.fasta) + '.renamed.fasta.tsv'
	
if args.sequence:
	include_sequence = True

check_seqs_for_standard = False
if args.standard:
	check_seqs_for_standard = True
	alphabet = args.standard
	out_warnings = out_fasta + '.warnings.txt'

dont_rename = False
if args.keep_old_names:
	dont_rename = True

prefix = args.prefix
add_prefix = False
if args.add_prefix:
	if prefix == "":
		print('Error: prefix cannot be blank if --addprefix was used.', file=sys.stderr)
		sys.exit(1)
	else:
		add_prefix = True
	

def check_standard(sequence, alphabet):
	check = True
	if alphabet == 'aa':
		standard_aas = 'ACDEFGHIKLMNPQRSTVWY'
		for letter in sequence:
			if letter.upper() not in standard_aas:
				check = False
	if alphabet == 'aax':
		standard_aas = 'ACDEFGHIKLMNPQRSTVWYX'
		for letter in sequence:
			if letter.upper() not in standard_aas:
				check = False
	if alphabet == 'aastar':
		standard_aas = 'ACDEFGHIKLMNPQRSTVWY*'
		for letter in sequence:
			if letter.upper() not in standard_aas:
				check = False
	if alphabet == 'aaxstar':
		standard_aas = 'ACDEFGHIKLMNPQRSTVWXY*'
		for letter in sequence:
			if letter.upper() not in standard_aas:
				check = False
	if alphabet == 'nt':
		standard_nts = 'ATCG'
		for letter in sequence:
			if letter.upper() not in standard_nts:
				check = False
	if alphabet == 'ntn':
		standard_nts = 'ATCGN'
		for letter in sequence:
			if letter.upper() not in standard_nts:
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
			tsv_header = 'seq_num\toriginal\tnew_name\tseqlen\tsequence\n'
		else:
			tsv_header = 'seq_num\toriginal\tnew_name\tseqlen\n'
		tsv_handler.write(tsv_header)
		print('Writing outputs to %s and %s...' % (out_fasta, out_tsv))
		seq_count = 0
		for seq in tqdm(fasta_d.keys()):
			seq_count += 1
			seq_len = len(fasta_d[seq])
			orig_header = seq
			seq_string = fasta_d[seq]
			if dont_rename:
				if add_prefix:
					output_name = '%s%s' % (prefix, orig_header)
				else:
					output_name = orig_header
			else:
				output_name = '%s%s' % (prefix, seq_count)
			fasta_handler.write('>%s\n%s\n' % (output_name, seq_string))
			if args.sequence:
				tsv_out = ('%s\t%s\t%s\t%i\t%s\n' % (seq_count, orig_header, output_name, seq_len, seq_string))
			else:
				tsv_out = ('%s\t%s\t%s\t%i\n' % (seq_count, orig_header, output_name, seq_len))
			tsv_handler.write(tsv_out)
print("Done!")
print("    FASTA written to: %s\n    TSV written to: %s" % (out_fasta, out_tsv) )
if check_seqs_for_standard:
	print('    Warnings written to: %s' % (out_warnings))
