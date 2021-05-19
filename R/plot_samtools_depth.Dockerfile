FROM rocker/r-base:4.0.4

LABEL authors="michaeljamesmansfield@gmail.com" \
	description="Docker image with all requirements to run the plot_samtools_depth.R script."

# Update apt repositories
RUN apt-get update \
	&& apt-get clean \ 
	&& rm -rf /var/lib/apt/lists/*

RUN R -e "install.packages(c('optparse', 'data.table', 'R.utils', 'viridis'), dependencies=TRUE, repos='https://cloud.r-project.org')"
