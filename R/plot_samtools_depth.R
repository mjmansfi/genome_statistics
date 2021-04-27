#!/bin/Rscript

######################################################################################
#
# Mike Mansfield 2021
# The input to this script is the output of samtools depth -a.
# It assumes that you have run samtools depth -a [sam/bam file] (which ensures that
# positions with 0 depth are reported in the output).
#
# If you haven't used the -a flag, basically everything is going to be wrong because
# the script assumes there is one row per position per sequence.
#
# What the script does is calculate the mean coverage over windows of size X.
# You can specify a step size if you want to do sliding windows; by default,
# it calculates coverage over chunks instead of over a sliding window.
# Meaning, with the defaults, it calculates coverage for every sequence from positions
# 1 to 5000, then 5001 to 10000, ... etc.
#
# Note that sequences with 0 depth at all positions are excluded from output plots.
#
######################################################################################

suppressPackageStartupMessages(require(tools))
suppressPackageStartupMessages(require(utils))
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(data.table))
suppressPackageStartupMessages(require(R.utils))
suppressPackageStartupMessages(require(splines))
suppressPackageStartupMessages(require(parallel))
suppressPackageStartupMessages(require(viridis))

####################################################################
#
# OPTPARSE INPUT OPTIONS
#
####################################################################
option_list = list(
	# Input options
	make_option(c('-d', '--depthFile'), action='store', default=NA, type='character',
		help='Path to file containing the output of samtools depth -a. Can be compressed. Default: none'),
	make_option(c('-w', '--windowSize'), action='store', default=1, type='integer',
		help='Window size to calculate mean coverage. Default: 1 (no window averaging)'),
	make_option(c('-s', '--stepSize'), action='store', default=NA, type='integer',
		help='An integer, which must be equal to or less than --windowSize, used to calculate sliding window averages. If equal (default), averages are calculated in segments of length --windowSize. If unequal, sliding window averages of size --windowSize are calculated for each --stepSize interval. Default: equal to --windowSize, no sliding window'),
	make_option(c('-f', '--filter'), action='store', default=NA, type='character',
		help='Comma-separated list of query sequences to analyze. Note that this is sensitive to white space, so the sequence names in the depth file must be simple. Default: all sequences analyzed'),

	# Output options
	make_option(c('-o', '--outFile'), action='store', default=NA, type='character',
		help='Base name to use for output files. Default: [depth file base name]'),
	make_option(c('-e', '--extension'), action='store', default='pdf', type='character',
		help='Output file type. Default: pdf'),
	make_option(c('--force'), action='store_true', default=FALSE, type='logical',
		help='Force overwrite of existing output files. Default: False'),

	# Plot options
	make_option(c('--scaleLenPerSeq'), action='store_true', default=FALSE, type='logical',
		help='Flag specifying whether sequence lengths should be scaled independently. Default: False. By default the X limit of each plot is scaled to the largest sequence in the data set.'),
	make_option(c('--scaleDepthPerSeq'), action='store_true', default=FALSE, type='logical',
		help='Flag specifying whether sequence depths should be scaled independently. Default: False. By default the Y limit of each plot is scaled to the maximum depth in the data set.'),
	make_option(c('--seqsPerPage'), action='store', default=1, type='integer',
		help='An integer between 1 and 5 specifying the number of sequences to include on each page of the coverage plot. Default: 1'),
	make_option(c('--plotDepthPerSeq'), action='store_true', default=FALSE, type='logical',
		help='Flag specifying whether or not to show plot summarizing mean depth for each sequence. Default: False.'),
	make_option(c('--plotDepthWindowDistPerSeq'), action='store_true', default=FALSE, type='logical',
		help='Flag specifying whether or not to show box plots with mean depth window distributions for each sequence. Default: False. By default coverage windows are not calculated; if they are, showPerSeqWindowPlot is FALSE unless specified.'),
	make_option(c('--viridisPalette'), action='store', default='viridis', type='character',
		help='A string specifying a viridis colour palette. This is only used in DepthPerSeq and DepthWindowDistPerSeq plots. Default: viridis'),
	make_option(c('--useSplines'), action='store_true', default=FALSE, type='logical',
		help='Flag specifying whether or not to use a smoothing spline in the coverage plot. Note: can look silly when combined with large values for --windowSize. Default: False'),

	# Runtime options
	make_option(c('-t', '--threads'), action='store', default=1, type='numeric',
		help='Number of threads for parallel processing. Default: 1'),
	make_option(c('-v', '--verbose'), action='store_true', default=FALSE, type='logical',
		help='Flag specifying whether to print verbose output to stdout. Prints summary stats and such. Default: False')
)

opt = parse_args(OptionParser(option_list=option_list))

# Make sure the depth file exists. If not, exit.
if(is.na(opt$depthFile)){
	cat('Error: an input samtools depth file must be specified with -d or --depthfile!\n', file=stderr())
	quit(status=1)
}

####################################################################
#
# INPUT SANITIZATION AND VARIABLE CREATION
#
####################################################################
############################################################
# Multithreading check. Needed for fread.
############################################################
num_threads = opt$threads

############################################################
# Depth file checking.
############################################################
# Read in the samtools depth file. Complain if it doesn't exist.
if(!file.exists(opt$depthFile)){
	cat(paste('Error: the input file specified by', opt$depthFile, 'does not exist!\n'), file=stderr())
	quit(status=1)
}

# This function makes sure that ALL positions are represented in the depth
# file. Which is to say, it ensures the -a flag was used with samtools depth.
ensure_depth_file_used_minus_a_flag <- function(depthDF){
	depthDF_by_seq = split(depthDF, depthDF$Sequence)
	test_sequential = lapply(depthDF_by_seq, function(seq){
		necessary_positions = seq(from=1, to=max(seq$Position), by=1)
		if(all(necessary_positions == seq$Position)){
			return(TRUE)
		} else {
			return(FALSE)
		}
	})
	if(all(unlist(test_sequential))){
		return(TRUE)
	} else {
		return(FALSE)
	}
}


############################################################
# Window and step size checking.
############################################################
# Check the window size.
if(opt$windowSize < 1){
	cat('Error: the --windowSize parameter must be a positive integer!\n', file=stderr())
	quit(status=1)
} else {
	window_size = opt$windowSize
}

# Check if a step size is specified.
# If it isn't set, then step size equals window size and sliding window mode is off.
# If it is set, make sure that the step size makes sense, or else complain.
if(is.na(opt$stepSize)){
	step_size = window_size
} else {
	if(opt$stepSize <= window_size){
		step_size = opt$stepSize
	} else {
		cat('Error: the --stepSize parameter cannot be larger than --windowSize! That would be silly.\n', file=stderr())
		quit(status=1)
	}
}

############################################################
# Filter statement checking.
############################################################
# Takes the comma-separated filter statement and turn it into a vector for later.
turn_filter_statement_into_vector <- function(filter_statement){
	filter_vector = unlist(strsplit(filter_statement, ','))
	return(filter_vector)
}

############################################################
# Output file checking.
############################################################
# If --outfile is not specified, use the depth file name as output.
if(is.na(opt$outFile)){
	output_base = basename(opt$depthFile)
} else{
	output_base = opt$outFile
}

# Check file extension options
if(opt$extension %in% c('pdf', 'png', 'jpeg', 'ps', 'bmp')){
	output_file_extension = opt$extension
} else{
	cat(paste('Error: the option provided to --extension is not recognized!\nAvailable options:    pdf, png, jpeg, ps, bmp\n'), sep='')
	quit(status=1)
}

# A simple check to overwrite existing output files
check_file_overwrite = function(fileName, forceOverwrite=F) {
	if(file.exists(fileName)){
		if(forceOverwrite){
			return(TRUE)
		} else {
			force_overwrite = menu(c('Yes', 'No'), title=paste('Warning: The output file', opt$outfile, 'already exists. Overwrite?'))
			if (force_overwrite == 1) {
				return(TRUE)
			} else {
				return(FALSE)
			}
		}
	}
}

# Check all output files to see if they exist.
depth_plot_file_name = paste(output_base, '_depthPlot', output_file_extension, sep='')
check_file_overwrite(depth_plot_file_name, forceOverwrite=opt$force)

if(opt$plotDepthPerSeq){
	depth_per_seq_plot_file_name = paste(output_base, '_depthPerSeqBarplot', output_file_extension, sep='')
	check_file_overwrite(depth_plot_file_name, forceOverwrite=opt$force)
}

if(opt$plotDepthPerSeq){
	depth_per_seq_plot_file_name = paste(output_base, '_depthPerSeqBarplot', output_file_extension, sep='')
	check_file_overwrite(depth_plot_file_name, forceOverwrite=opt$force)
}

############################################################
# Plot parameter checking.
############################################################
# Whether or not to apply global sequence length scaling vs. per-sequence.
if(opt$scaleLenPerSeq){
	scale_len_per_seq = TRUE
} else {
	scale_len_per_seq = FALSE
}

# Whether or not to apply global sequence length scaling vs. per-sequence.
if(opt$scaleDepthPerSeq){
	scale_depth_per_seq = TRUE
} else {
	scale_depth_per_seq = FALSE
}

# Number of seqs to plot on each page.
if(opt$seqsPerPage == 1){
	seqs_per_page = 1
} else if (opt$seqsPerPage < 6) {
	seqs_per_page = opt$seqsPerPage
} else {
	cat('Error: the --seqsPerPage parameter must be an integer less than 5! More than 3 might break things. More than 5 almost certainly will.\n', file=stderr())
	quit(status=1)
}

# Which viridis palette to use.
if(!opt$viridisPalette %in% c('viridis', 'magma', 'inferno', 'plasma', 'cividis', 'rocket', 'mako', 'turbo')){
	cat(paste('Error: the argument provided to --viridisPalette is not recognized!\n', 'Available options: ', sep=''), file=stderr())
	quit(status=1)
} else {
	viridis_palette = opt$viridisPalette
}

# Whether or not to use splines.
if(opt$useSplines){
	use_splines = TRUE
} else {
	use_splines = FALSE
}

############################################################
# Summarize inputs and outputs.
############################################################
print_summary_stats <- function(depthDF) {
	number_of_rows = nrow(depthDF)
	number_of_rows_above_0 = nrow(subset(depthDF, depthDF$Depth > 0))
	number_of_sequences = length(unique(depthDF$Sequence))
	rows_per_sequence = sort(table(depthDF$Sequence), decreasing=T)
	rows_per_sequence <- apply(data.frame(head(rows_per_sequence, n=5)), 1, function(x) gsub('[[:space:]]', '', paste(x[1], x[2], sep=':')))
	depth_mean = format(mean(depthDF$Depth), digits=2)
	depth_median = format(median(depthDF$Depth), digits=2)
	depth_sd = format(sd(depthDF$Depth), digits=2)

	stat_summary = paste(
		'INPUT OPTIONS:', '\n',
		'  Input depth file:                  ', opt$depthFile, '\n',
		'  Filter statement:                  ', turn_filter_statement_into_vector(opt$filter), '\n',
		'  Window size:                       ', window_size, '\n',
		'  Step size:                         ', step_size, '\n',
		'  Number of threads:                 ', num_threads, '\n',
		'OUTPUT OPTIONS:', '\n',
		'  Output base name:                  ', output_base, '\n',
		'  Output file type:                  ', output_file_extension, '\n',
		'  Plot per-seq depth barplot?        ', opt$plotDepthPerSeq, '\n',
		'  Plot per-seq window distribution?  ', opt$plotDepthWindowDistPerSeq, '\n',
		'  Plot coverage using splines?       ', opt$useSplines, '\n',
		'SUMMARY OF DEPTH FILE:', '\n',
		'  Total number of positions:         ', format(number_of_rows, big.mark=','), '\n',
		'  Num. positions with coverage >0:   ', format(number_of_rows_above_0, big.mark=','), '\n',
		'  Overall cov. mean, median, and SD: ', paste(depth_mean, depth_median, depth_sd, sep=', '), '\n',
		'  Top 5 sequences by length:         ', paste(rows_per_sequence[1:5], collapse=', '), '\n',
		sep='')
	cat(stat_summary)
}


####################################################################
#
# CALCULATION AND PLOTTING FUNCTIONS
#
####################################################################
############################################################
# A helper function to remove sequences that have a max
# of 0 coverage (i.e., 0 depth at all positions).
############################################################
remove_sequences_with_0_depth = function(depthDF){
	depthDF_by_seq = split(depthDF, depthDF$Sequence)
	depthDF_by_seq = depthDF_by_seq[sapply(depthDF_by_seq, function(x) max(x$Depth)>0)]
	return(depthDF_by_seq)
}

############################################################
# A helper function to open output files in the right
# format.
############################################################
open_plot_output_file = function(fileExtension=NA, fileName=NA) {
	if(!is.na(fileName)){
		if(fileExtension == 'pdf'){
			file_function = pdf(fileName)
		} else if (fileExtension == 'png') {
			file_function = png(fileName)
		} else if (fileExtension == 'jpeg') {
			file_function = jpeg(fileName)
		} else if (extension == 'ps'){
			file_function = postscript(fileName)
		} else if (fileExtension == 'bmp'){
			file_function = bmp(fileName)
		}
	}
}

############################################################
# This function calculates depth statistics for each
# sequence in a depth file, without fancy sliding windows
# or anything.
############################################################
calculate_depth_per_seq <- function(depthDF, numThreads=1){
	# Remove sequences with 0 depth at all positions.
	depthDF_by_seq = remove_sequences_with_0_depth(depthDF=depthDF)
	depthDF_by_seq = split(depthDF, depthDF$Sequence)

	# For each sequence...
	list_of_chunked_sequences = mclapply(mc.cores=numThreads, depthDF_by_seq, function(seq){
		# Calculate statistics...
		seq_len = nrow(seq)
		seq_name = unique(seq$Sequence)

		df_with_0s = data.frame(Position=1:seq_len)
		df_with_0s = merge(df_with_0s, seq, by='Position', all=T)
		df_with_0s$Depth[is.na(df_with_0s$Depth)] = 0

		depth_mean = mean(df_with_0s$Depth)
		depth_median = median(df_with_0s$Depth)
		depth_sd = sd(df_with_0s$Depth)

		depth_list = list(Sequence=seq_name, Length=seq_len, Depth_mean=depth_mean, Depth_median=depth_median, Depth_sd=depth_sd)
		return(depth_list)
	})

	joined_chunks = as.data.frame(do.call(rbind, list_of_chunked_sequences))

	joined_chunks$Sequence = as.character(unlist(joined_chunks$Sequence))
	joined_chunks$Length = as.numeric(unlist(joined_chunks$Length))
	joined_chunks$Depth_mean = as.numeric(unlist(joined_chunks$Depth_mean))
	joined_chunks$Depth_sd = as.numeric(unlist(joined_chunks$Depth_sd))
	joined_chunks$Depth_median = as.numeric(unlist(joined_chunks$Depth_median))
	joined_chunks = joined_chunks[order(joined_chunks$Length, decreasing=T),]
	rownames(joined_chunks) <- NULL

	# Return a formatted dataframe
	return(joined_chunks)
}


############################################################
# This is a helper function that takes a numeric vector
# and splits it into X breaks, returning an annotated
# dataframe of the same length with $colour and $label
# columns. This makes it easy to pass colour annotations to
# to plotting functions.
############################################################
numeric_vector_to_labels_and_colours = function(numeric_vector, breaks=8, viridis_palette=viridis, prefix='', suffix='', digits=2){
	colour_palette = viridis_palette(breaks)
	cut_vector = cut(numeric_vector, breaks=breaks)
	cut_vector_levels = levels(cut_vector)

	bounds = as.data.frame(cbind(lower = as.numeric(sub('\\((.+),.*', '\\1', cut_vector_levels)), upper = as.numeric(sub('[^,]*,([^]]*)\\]', '\\1', cut_vector_levels))))
	bounds$lower[1] = min(numeric_vector)
	bounds$upper[nrow(bounds)] = max(numeric_vector)
	bounds$label = trimws(paste(prefix, ' ', trimws(format(round(bounds$lower, digits=digits))), '-', trimws(format(round(bounds$upper, digits=digits))), ' ', suffix, sep=''))
	bounds$colour = colour_palette

	vector_with_colours_and_labels = data.frame(vec=numeric_vector)
	vector_with_colours_and_labels$cut_factor = cut_vector
	vector_with_colours_and_labels$colour = bounds$colour[as.numeric(cut(numeric_vector, breaks=breaks))]
	vector_with_colours_and_labels$label = bounds$label[as.numeric(cut(numeric_vector, breaks=breaks))]
	return(vector_with_colours_and_labels)
}


############################################################
# This function plots coverage per sequence.
############################################################
plot_depth_per_seq <- function(chunkedDF, seqs_per_plot=25, plots_per_page=1, viridis_palette=viridis, breaks=8){
	# Call helper function to get colours and labels according to a numeric vector (length).
	labels_and_colours = numeric_vector_to_labels_and_colours(chunkedDF$Length/1000, breaks=breaks, viridis_palette=viridis_palette, prefix='', suffix='kbp', digits=2)
	chunkedDF$colour = labels_and_colours$colour
	chunkedDF$label = labels_and_colours$label
	# Split and reorder by length
	split_by_seq = split(chunkedDF, chunkedDF$Sequence)
	split_by_seq = split_by_seq[order(sapply(split_by_seq, function(x) max(x$Length)), decreasing=T)]
	split_start_points = seq(from=1, to=length(split_by_seq), by=seqs_per_plot)

	# Making plots
	max_depth = max(chunkedDF$Depth_mean)
	par(mfrow=c(plots_per_page, 1))
	A = lapply(split_start_points, function(x){
		list_start = x
		if((list_start+seqs_per_plot) < length(split_by_seq)){
			list_end = x + seqs_per_plot
		} else {
			list_end = length(split_by_seq)
		}

		list_subsection = split_by_seq[c(list_start:list_end)]
		joined_df = do.call(rbind, list_subsection)

		barplot(joined_df$Depth_mean ~ joined_df$Sequence, las=2, xlab='Sequence', ylab='Depth', ylim=c(0, max_depth), col=joined_df$colour)
		legend('top', horiz=T, pch=15, col=unique(chunkedDF$colour), legend=unique(chunkedDF$label), bty='n')

	})
	par(mfrow=c(1,1))
}


############################################################
# This function splits a depth data frame into a list of
# data frames based on sequence name. Then, for each
# sequence, it calculates mean depth over a window of
# size --windowsize. Optionally calculates mean depth over
# a sliding window instead of static chunks.
# A data frame including these data is returned.
#
# Note that the function assumes that the columns are in
# exactly the same order and named as if they were parsed
# by this script.
# The samtools depth column names should be:
#     c('Sequence', 'Position', 'Depth')
############################################################
calculate_depth_windows <- function(depthDF, windowSize=5000, stepSize=5000, numThreads=1){
	# Remove sequences with 0 depth at all positions.
	depthDF_by_seq = remove_sequences_with_0_depth(depthDF=depthDF)
	depthDF_by_seq = split(depthDF, depthDF$Sequence)

	# For each sequence...
	list_of_chunked_sequences = lapply(depthDF_by_seq, function(seq){
		seq_len = nrow(seq)
		seq_name = unique(seq$Sequence)

		# Break sequence into chunks...
		chunk_start_points = seq(from=1, to=seq_len, by=stepSize)

		chunk_data_frames = mclapply(mc.cores = numThreads, chunk_start_points, function(chunk_start){
			# Then calculate depth statistics.
			if(chunk_start + windowSize <= seq_len){
				chunk_end = chunk_start + windowSize
			} else {
				chunk_end = seq_len
			}
			chunk_mid = (chunk_start+chunk_end)/ 2

			current_chunk = subset(seq, chunk_start <= seq$Position & seq$Position < chunk_end)

			df_with_0s = data.frame(Position=chunk_start:chunk_end)

			df_with_0s = merge(df_with_0s, current_chunk, by='Position', all=T)
			df_with_0s$Depth[is.na(df_with_0s$Depth)] = 0

			depth_mean = mean(df_with_0s$Depth)
			depth_median = median(df_with_0s$Depth)
			depth_sd = sd(df_with_0s$Depth)

			chunk_list = list(Sequence=seq_name, Start=chunk_start, End=chunk_end, Midpoint=chunk_mid, Depth_mean=depth_mean, Depth_median=depth_median, Depth_sd=depth_sd)
			return(chunk_list)
		})

		joined_chunks = as.data.frame(do.call(rbind, chunk_data_frames))

		joined_chunks$Sequence = as.character(unlist(joined_chunks$Sequence))
		joined_chunks$Start = as.numeric(unlist(joined_chunks$Start))
		joined_chunks$End = as.numeric(unlist(joined_chunks$End))
		joined_chunks$Midpoint = as.numeric(unlist(joined_chunks$Midpoint))
		joined_chunks$Depth_mean = as.numeric(unlist(joined_chunks$Depth_mean))
		joined_chunks$Depth_sd = as.numeric(unlist(joined_chunks$Depth_sd))
		joined_chunks$Depth_median = as.numeric(unlist(joined_chunks$Depth_median))
		return(joined_chunks)
	})
	concatenated_df_of_chunked_sequences = do.call(rbind, list_of_chunked_sequences)
	rownames(concatenated_df_of_chunked_sequences) = NULL
	return(concatenated_df_of_chunked_sequences)
}


############################################################
# This function plots depth window distributions
# across all sequences in the depthfile.
############################################################
plot_depth_window_distribution_per_seq <- function(chunkedDF, seqs_per_plot=20, plots_per_page=1, viridis_palette=viridis, breaks=8){
	# Call helper function to get colours and labels according to a numeric vector (length).
	labels_and_colours = numeric_vector_to_labels_and_colours(chunkedDF$Length/1000, breaks=breaks, viridis_palette=viridis_palette, prefix='', suffix='kbp', digits=2)
	chunkedDF$colour = labels_and_colours$colour
	chunkedDF$label = labels_and_colours$label

	split_by_seq = split(chunkedDF, chunkedDF$Sequence)
	# Order by sequence size (maximum end size)
	split_by_seq = split_by_seq[order(sapply(split_by_seq, function(x) max(x$End)), decreasing=T)]
	split_start_points = seq(from=1, to=length(split_by_seq), by=seqs_per_plot)

	# Plot parameters
	min_depth = min(chunkedDF$Depth_mean)

	max_depth = max(chunkedDF$Depth_mean)

	par(mfrow=c(plots_per_page, 1))
	A = lapply(split_start_points, function(x){
		list_start = x
		if( (list_start+seqs_per_plot) < length(split_by_seq)){
			list_end = x + seqs_per_plot
		} else {
			list_end = length(split_by_seq)
		}

		list_subsection = split_by_seq[c(list_start:list_end)]
		joined_df = do.call(rbind, list_subsection)

		boxplot( joined_df$Depth_mean ~ joined_df$Sequence, las=2, xlab='Sequence', ylab='Depth', ylim=c(0, max_depth), col=joined_df$colour)
		legend('topright', horiz=F, pch=15, col=unique(chunkedDF$colour), legend=unique(chunkedDF$label), bty='n')
	})
	par(mfrow=c(1,1))
}


############################################################
# This function plots coverage across the length of each
# sequence in the depth file.
############################################################
plot_depth_windows_over_seqs <- function(chunkedDF, plots_per_page=1, windows=F, splines=F, per_seq_depth_scaling=T, per_seq_length_scaling=F){
	split_by_seq = split(chunkedDF, chunkedDF$Sequence)
	# Order by sequence size (maximum end size)
	split_by_seq = split_by_seq[order(sapply(split_by_seq, function(x) max(x$End)), decreasing=T)]

	max_depth_across_seqs = max(chunkedDF$Depth_mean)
	max_length_across_seqs = max(chunkedDF$End)

	par(mfrow=c(plots_per_page, 1))
	# For every sequence...
	A = lapply(split_by_seq, function(seq){
			# Plot coverage. First, determine X and Y limits based on function arguments.
			if(per_seq_length_scaling){
				max_x = max(seq$End)
			} else {
				max_x = max_length_across_seqs
			}
			if(per_seq_depth_scaling){
				max_y = max(seq$Depth_mean)
			} else {
				max_y = max_depth_across_seqs
			}

			# If averaging over a window, plot the coverage as blocks.
			if(windows){
				if(splines) {
					plot(type='n', x=1, y=1, xlim=c(0, max_x), ylim=c(0, max_y), ylab='Depth', xlab='Window start position', main=unique(seq$Sequence))
					spline = smooth.spline(seq$Depth_mean ~ seq$Start)
					x_values = c(0, spline$x, max(spline$x))
					y_values = c(0, spline$y, 0)
					polygon(x=x_values, y=y_values, col='#595959')
				} else {
					plot(type='n', x=1, y=1, xlim=c(0, max_x), ylim=c(0, max_y), ylab='Depth', xlab='Window start position', main=unique(seq$Sequence))
					rect(xleft=seq$Start, xright=seq$End, ybottom=0, ytop=seq$Depth_mean, col='#595959', border=NA)
				}
			# If not averaging over a window, plot either the lines themselves or their spline.
			} else {
				if(splines){
					spline = smooth.spline(seq$Depth_mean ~ seq$Start)
					x_values = c(0, spline$x, max(spline$x))
					y_values = c(0, spline$y, 0)
				} else {
					x_values = c(0, seq$Start, max(seq$Start))
					y_values = c(0, seq$Depth_mean, 0)
				}
				plot(type='n', x=1, y=1, xlim=c(0, max_x), ylim=c(0, max_y), ylab='Depth', xlab='Position', main=unique(seq$Sequence))
				polygon(x=x_values, y=y_values, col='#595959')
			}
		})
	par(mfrow=c(1,1))
}


####################################################################
#
# WORKFLOW
#
####################################################################
depth_df <- fread(opt$depthFile, nThread=num_threads)
# Now that the depth DF is loaded, give the columns some names.
colnames(depth_df) <- c('Sequence', 'Position', 'Depth')

if(!ensure_depth_file_used_minus_a_flag(depthDF=depth_df)){
	cat(paste('Error: not all positions are represented for every sequence in the depth file!\n'), file=stderr())
	quit(status=1)
}

# If there is a filter statement, subset the depth file to include only matching sequences.
if(!is.na(opt$filter)){
	seqs_of_interest <- turn_filter_statement_into_vector(opt$filter)
	depth_df = subset(depth_df, depth_df$Sequence %in% seqs_of_interest)
}
# Make sure that there are still sequences left after filtering.
if(nrow(depth_df) == 0){
	cat(paste('Error: No sequences are left after filtering based on your filter statement!\nMake sure that the sequence(s) provided to --filter are in the depth file.\nYour --filter argument:    ', opt$filter), file=stderr())
	quit(status=1)
}

if(opt$verbose){
	print_summary_stats(depthDF=depth_df)
}
