###########################################################
#
#	Splits a multiple fasta into many individual fasta
#	files. File names are based on the fasta headers, and
#	there is some minor reformatting in the script (e.g.,
#	spaces, slashes and colons become underscores), but if 
#	there are other characters that result in non-UNIX
#	compliant file names (e.g., brackets or length too long),
#	you must edit this first.
#
###########################################################

#!/usr/bin/env/bash
if [[ $# -eq 0 ]] ; then
        echo "Error: you must specify a fasta file with -f."
        exit 1
fi

#collects input flags
case "$1" in
	-f)
		INFILE="$2"
		shift
	;;
	*)
		echo "Error: unknown option '$1'. Please use the -f flag only."
		exit 1
	;;
esac

#checks if input is a file
if [ -e "$INFILE" ] ;
then
	echo "Splitting the multiple fasta file $INFILE into individual fasta files..."	
else
	echo "Error: please specify a file with -f [filename]."
	exit 1
fi

while read LINE
do
	if [[ ${LINE:0:1} == '>' ]]
	then
		OUTHANDLE=`echo ${LINE#>} | sed 's/ /_/g' | sed 's/\//_/g' | sed 's/\:/_/g'`
		#echo $OUTHANDLE
		OUTFILE="${OUTHANDLE#>}.fa"
		#echo "    "$OUTFILE
		echo $LINE > "$OUTFILE"
	else
		echo $LINE >> "$OUTFILE"
	fi
done < $INFILE
