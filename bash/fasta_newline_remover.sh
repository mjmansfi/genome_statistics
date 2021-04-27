#!/usr/bin/env bash
# 	Note, sometimes you need to remove stupid DOS newlines! use dos2unix utility to do this

if [[ $# -eq 0 ]] ; then
        echo "Error: you must specify a fasta file with -f."
        exit 1
fi

case "$1" in
	-f)
		FILE="$2"
		shift
	;;
	*)
		echo "Error: unknown option '$1'. Please use the -f flag only."
		exit 1
	;;
esac

if [[ -z "$FILE" ]] ; then
	echo "Error: please specify a file with -f [filename]."
	exit 1
fi

awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' $FILE
