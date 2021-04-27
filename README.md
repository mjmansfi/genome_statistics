# genome_statistics #
A collection of useful scripts for analyzing genomic data.

Hey, it's better than random bash scripts, no?

## plot_samtools_depth.R

A script to plot different parameters based on the output of *samtools depth -a [sorted, indexed bamfile]*. The **-a** flag is very important because ALL positions for ALL sequences need to be represented in the data. There is a rudimentary check for this included in the script, so the script should complain if you do not use this flag.

Non-base library dependencies:
- `optparse` - to parse input arguments
- `data.table` - to read input files
- `R.utils` - to read compressed input files
- `viridis` - to allow more colourful plots

### Required arguments:
#### -d, --depthFile

Path to a file containing the output of `samtools depth -a`. Can be gzip compressed. The output should look like:
```bash
seq1	1	5
seq1	2	6
seq1	3	0
seq1	4	0
...
seqN	X	Y
```
Note that sequence names may be truncated in the process of creating the depth file. To avoid this, ensure that the sequence names in the FASTA file are simple (e.g., by using `python/fasta_rename_with_table.py`!)


### Optional arguments:
#### -o, --outFile
Base name for output files. If not specified, the default is `[depth file input name]`. The script automatically appends the strings `_depthPlot.pdf`, `_depthPerSeqBarplot.pdf`, `_`

#### -w, --windowSize
Size of the window to calculate average depth values. By default, the script does not calculate any windows (because by default `--windowSize` equals 1), but this may result in a noisy coverage plot. Any non-negative number is acceptable; if a window's end is greater than the length of a sequence in the FASTA file, the end is the last position.

#### -s, --stepSize
Size of the steps between window intervals, necessary for calculating sliding window averages. By default, `--windowSize` is 1 and no windowed averages are calculated. If `--windowSize` is specified but `--stepSize` is not, by default, `--stepSize` equals `--windowSize` (no sliding windows). Sliding windows are only calculated when `--windowSize` is greater than 1 and `--stepSize` is manually specified number less than `--windowSize`.

#### -f, --filter
Comma-separated list of query sequences to analyze. Useful if you are interested in one or a few sequences represented in the depth file. Note that this filter statement is sensitive to strange characters and whitespace, so if you want to use this feature make sure that the sequence names in the depth file are simple and easy to specify. For example: `--filter seq1,seq5`

#### --scaleLenPerSeq, --scaleDepthPerSeq
Scaling options that determine the limits of each plot on the X (length) and Y (Depth) axes. The following image (which, for now, is manually coloured) shows what these options do more clearly:

![Scaling options](https://raw.githubusercontent.com/mjmansfi/genomics_scripts/main/assets/plot_samtools_depth_scaling.png)

#### --seqsPerPage
Option to include multiple sequences on a single PDF page. By default, one sequence is shown per page in the depth plot. Increasing this is mostly useful if long, thin plots is preferable to the square ones produced by default. However, this number cannot be more than 5, because base R will start to complain about margin sizes. I recommend post-processing in [Inkscape](https://inkscape.org/) or similar.

Do note that this affects the way barplots are distributed if `--plotDepthPerSeq` is specified.

#### --plotDepthPerSeq
Optional barplot output displaying average coverage per sequence in the depth file, coloured by sequence length. Example:
![Scaling options](https://raw.githubusercontent.com/mjmansfi/genomics_scripts/main/assets/plot_samtools_depth_plotDepthPerSeq.png)

Colours can be changed with the `--viridisPalette` argument.

#### --useSplines
Fit smoothing splines to coverage distributions instead of raw counts (or window counts). Example:
![Scaling options](https://raw.githubusercontent.com/mjmansfi/genomics_scripts/main/assets/plot_samtools_depth_useSplines.png)

The splines sometimes make it easier to spot trends in the data, and also drastically reduces the number of vertices in the polygon, making it easier to manipulate in external image editing programs. But, sometimes it looks a bit silly, especially for low-coverage sequences.

#### -t, --threads
Number of threads to use for parallel processing.

#### --writeDepthPerSeqTable, --writeDepthWindowTable
Optionally write tables of per-sequence coverage or per-window coverage to external `tsv` files.

#### --verbose
Print nice summary messages.

#### --force
Force overwrite of existing output files.
