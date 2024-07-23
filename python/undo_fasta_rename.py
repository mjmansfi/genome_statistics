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
parser.add_argument('-f','--fasta', help='FASTA file containing renamed sequences. Default: none', required=True)
parser.add_argument('-tsv','--tsv', help='TSV table containing sequences that were at one point renamed.\nDefault: none', required=True)
parser.add_argument('-o','--outfasta', help='Renamed FASTA output file.\nDefault: [input fasta].unrenamed.fasta', required=False)

args = parser.parse_args()

if args.outfasta:
	out_fasta = args.outfasta
else:
	out_fasta = os.path.basename(args.fasta) + '.unrenamed.fasta'

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

tsv_d = {}
with open(args.tsv) as i:
	print('Reading in TSV file %s...' % args.tsv)
	for line in tqdm(i):
		line = line.strip()
		if not line.startswith('seq_num'):
			split_line = line.split('\t')
			seq_num = split_line[0]
			original = split_line[1]
			renamed = split_line[2]
			seqlen = split_line[3]
			tsv_d[renamed] = original

unrenamed_d = {}
for seq in fasta_d:
	if seq in tsv_d:
		original_name = tsv_d[seq]
		unrenamed_d[original_name] = fasta_d[seq]

with open(out_fasta, 'w') as fasta_handler:
	print('Writing outputs to %s...' % (out_fasta))
	seq_count = 0
	for seq in tqdm(unrenamed_d.keys()):
		seq_string = unrenamed_d[seq]
		fasta_handler.write('>%s\n%s\n' % (seq, seq_string))

print("Done!")
print("    FASTA written to: %s\n" % (out_fasta) )


