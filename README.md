# genome_statistics
A collection of useful scripts for analyzing genomic data.

Hey, it's better than random bash scripts, no?

##plot_samtools_depth.R

A script to plot different parameters based on the output of *samtools depth -a [sorted, indexed bamfile]*. The **-a** flag is very important because ALL positions for ALL sequences need to be represented in the data. There is a rudimentary check for this included in the script, so the script should complain if you do not use this flag.

Required options:

###-d, --depthFile

Path to the file containing the output of samtools depth -a. Can be gzip compressed. The output should look like:
```bash
seq1	1	5
seq1	2	6
seq1	3	0
seq1	4	0
...
seqN	X	Y
```
Note that sequence names may be truncated in the process of creating the depth file. To avoid this, ensure that the sequence names in the FASTA file are simple.

###-f, --filter

Comma-separated list of query sequences to analze. Useful if you are interested in one or a few sequences represented in the depth file. Usage: `--filter seq1,seq5`

###-o, --outFile

Output file name. If not specified, the default is `[depth file input name].pdf`. If `--outFile` is specified, the script appends the file extension, if it is not already specified.
