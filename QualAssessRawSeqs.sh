#!/bin/bash
#  @/Author: Christopher A. Gulvik, Ph.D.
#  +/Version 1.1
#  +/6 February 2015
#  Dependencies:  awk, bc, cut, GNU sed, grep, gunzip, mkdir, readlink, paste, perl, printf, R, wc
#    script $HOME/scripts/FastQ2fastaQual.pl
#    script $HOME/scripts/qscore_graphs_illumina.pl

function usage() { 
	echo "
	USAGE: $0 Input_R1.fastq Input_R2.fastq Outdir

	The raw reads must be in Illumina 1.8+ FastQ
	format, and filenames must be formatted:
	<name>.fastq.gz or <name>.fastq, where
	<name> is the name of the sample.

	Absolute paths are required. Relative paths
	are not accepted.
	"
	}

#Help and usage information
if [[ "$1" == "" || "$1" == "--help" || "$1" == "-h" ]]; then
	usage
	exit 1
fi

#Require 3 arguments
if [ $# -ne 3 ]; then
	usage
	exit 1
fi

#Make outdir if absent
if [ ! -d "$3" ]; then
	mkdir -p "$3"
	echo "created output directory path:  $3"
fi

###Filename handling
dIFS=$IFS  #default IFS
IFS=$''  #change for newlines
#takes raw unpaired R1 and R2 FastQ files as input
#tests file input extensions
if [[ $1 == *.fastq.gz && $2 == *.fastq.gz ]]; then  #both files are gunzipped
	echo "extracting the gunzipped file $1..."
	gunzip -k $1
	echo "$1 was extracted..."
	SEQ1NAME=$(ls $1 | cut -d . -f 1-2)
	SEQ1=$(echo $SEQ1NAME)
	echo "extracting the gunzipped file $2..."
	gunzip -k $2
	echo "$2 was extracted..."
	SEQ2NAME=$(ls $2 | cut -d . -f 1-2)
	SEQ2=$(echo $SEQ2NAME)
elif [[ $1 == *.fastq && $2 == *.fastq ]]; then  #both files are fastq format
	echo "both are fastq..."
	SEQ1=$(echo $1)
	SEQ2=$(echo $2)
elif [[ $1 == *.fastq.gz && $2 == *.fastq ]]; then
	echo "extracting the gunzipped file $1..."
	gunzip -k $1
	echo "$1 was extracted..."
	SEQ1NAME=$(ls $1 | cut -d . -f 1-2)
	SEQ1=$(echo $SEQ1NAME)
	SEQ2="$2"  #$2 is still in fastq format
elif [[ $1 == *.fastq && $2 == *.fastq.gz ]]; then
	SEQ1="$1"  #$1 is still in fastq format
	echo "extracting the gunzipped file $2..."
	gunzip -k $2
	echo "$2 was extracted..."
	SEQ2NAME=$(ls $2 | cut -d . -f 1-2)
	SEQ2=$(echo $SEQ2NAME)
else
	usage
	exit 1
fi
IFS=$dIFS  #back to default IFS
wait

###File input sequence information
#file regex from Illumina 1.8+ (gunzip extracted):
#  ([0-9A-Za-z]{1,19})_S([0-9A-Za-z]{1,3})_L([0-9]{3})_R[12]_([0-9]{3}).fastq
#trim filenames to something meaningful but succinct
SAMPLE=$(echo $SEQ1 | awk -F / '{print $(NF +0)}' | cut -d _ -f 1-2)  #grabs all characters before the 1st and underscore
SAMPLEA=$(echo $SEQ1 | awk -F / '{print $(NF +0)}'| cut -d _ -f 1,2,4)  #grabs full filename regardless of directories in path
SAMPLEB=$(echo $SEQ2 | awk -F / '{print $(NF +0)}'| cut -d _ -f 1,2,4)  #grabs full filename regardless of directories in path
wait
echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
echo " $SAMPLEA and $SAMPLEB sequences were read in..."
echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

###Count total nucleotides
TN1=$(awk 'BEGIN {SUM=0;} {if(NR%4==2) {SUM+=length($0);}} END {print SUM;}' $SEQ1)
TN2=$(awk 'BEGIN {SUM=0;} {if(NR%4==2) {SUM+=length($0);}} END {print SUM;}' $SEQ2)
BP_TOT="$(($TN1 + $TN2))"
wait
echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
echo " counted the total nucleotides..."
echo "        Total: $BP_TOT"
echo "        R1: $TN1"
echo "        R2: $TN2"
echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

OPSYS=$(uname -s)

###Count nucleotides >=Q20 and >=Q30 for Illumina 1.8+
if [[ "$OPSYS" == Linux ]]; then  #should be linux where sed default is GNU sed
	Q201=$(sed -n '4~4'p $SEQ1 | grep -o '[0-9:\;\<=\>\?@A-J]' | paste -s -d"\0" - | wc -m)
	Q202=$(sed -n '4~4'p $SEQ2 | grep -o '[0-9:\;\<=\>\?@A-J]' | paste -s -d"\0" - | wc -m) 
	Q301=$(sed -n '4~4'p $SEQ1 | grep -o '[\?@A-J]' | paste -s -d"\0" - | wc -m)
	Q302=$(sed -n '4~4'p $SEQ2 | grep -o '[\?@A-J]' | paste -s -d"\0" - | wc -m)
elif [[ "$OPSYS" == Darwin ]]; then  #use gsed instead of sed in Mac OS X here
	Q201=$(gsed -n '4~4'p $SEQ1 | grep -o '[0-9:\;\<=\>\?@A-J]' | paste -s -d"\0" - | wc -m)
	Q202=$(gsed -n '4~4'p $SEQ2 | grep -o '[0-9:\;\<=\>\?@A-J]' | paste -s -d"\0" - | wc -m) 
	Q301=$(gsed -n '4~4'p $SEQ1 | grep -o '[\?@A-J]' | paste -s -d"\0" - | wc -m)
	Q302=$(gsed -n '4~4'p $SEQ2 | grep -o '[\?@A-J]' | paste -s -d"\0" - | wc -m)
else
	echo 'ERROR: Your operating system does not appear to be supported.' >&2
	exit 1
fi
wait
Q20_TOT="$(($Q201 + $Q202))" 
Q30_TOT="$(($Q301 + $Q302))" 
wait
echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
echo ' counted the Q20 and Q30 nucleotides...'
echo "        ${Q201} nucleotides are >=Q20 in $SAMPLEA" 
echo "        ${Q301} nucleotides are >=Q30 in $SAMPLEA" 
echo "        ${Q202} nucleotides are >=Q20 in $SAMPLEB" 
echo "        ${Q302} nucleotides are >=Q30 in $SAMPLEB" 
echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

###Calculate percentages of >=Q20 and >=Q30
PQ201="$(echo "scale=2;($Q201/$TN1)*100" | bc)"
PQ202="$(echo "scale=2;($Q202/$TN2)*100" | bc)"
PQ301="$(echo "scale=2;($Q301/$TN1)*100" | bc)"
PQ302="$(echo "scale=2;($Q302/$TN2)*100" | bc)"
wait
echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
echo ' calculated the percentages of Q20 and Q30...'
echo "        ${PQ201}% nucleotides are >=Q20 in $SAMPLEA" 
echo "        ${PQ301}% nucleotides are >=Q30 in $SAMPLEA" 
echo "        ${PQ202}% nucleotides are >=Q20 in $SAMPLEB" 
echo "        ${PQ302}% nucleotides are >=Q30 in $SAMPLEB" 
echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

###Count number of reads
NUMLN_L=$(wc -l $SEQ1| awk '{print $1}')
NUMLN_R=$(wc -l $SEQ2 | awk '{print $1}')
wait
NUM_L="$(($NUMLN_L / 4))"
NUM_R="$(($NUMLN_R / 4))"
READ_TOT="$(($NUM_L + NUM_R))"
wait
echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
echo ' counted the number of reads...'
echo "        There are $NUM_L reads in $SAMPLEA" 
echo "        There are $NUM_R reads in $SAMPLEB" 
echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

########################################################
#unfixed bug in one of the perl dependencies, so this section is commented out for now
###Boxplot visualization###
###TOM###
#perl $HOME/scripts/FastQ2fastaQual.pl "$SEQ1"
#perl $HOME/scripts/FastQ2fastaQual.pl "$SEQ2"
#rm *.fasta
#echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
#echo ' created qual files...'
#echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
#QUALFILE1=$(echo "$SEQ1")
#QUALFILE2=$(echo "$SEQ2")
#perl $HOME/scripts/qscore_graphs_illumina.pl $QUALFILE1.qual
#perl $HOME/scripts/qscore_graphs_illumina.pl $QUALFILE2.qual
#rm *.qual
#wait
#echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
#echo ' created quality graphs...'
#echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
########################################################

#Create results output file
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
Sample_Name Q20_Total_[bp] Q30_Total_[bp] Q20_R1_[bp] Q20_R2_[bp] Q20_R1_[%] Q20_R2_[%] Q30_R1_[bp] Q30_R2_[bp] Q30_R1_[%] Q30_R2_[%] Total_Sequenced_[bp] Total_Sequenced_[reads] NUMR1_[bp] NUMR2_[bp] \
"$SAMPLE" "$Q20_TOT" "$Q30_TOT" "$Q201" "$Q202" "$PQ201%" "$PQ202%" "$Q301" "$Q302" "$PQ301%" "$PQ302%" "$BP_TOT" "$READ_TOT" "$NUM_L" "$NUM_R" > "$3"/QualAssessRawSeqs_"$SAMPLE"_results.tsv ;
wait
echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
echo ' created results output file...'
echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

if [[ "$OPSYS" == Linux ]]; then
	R1FILEPATH=$(readlink -m $1)
	R2FILEPATH=$(readlink -m $2)
elif [[ "$OPSYS" == Darwin ]]; then
	R1FILEPATH=$(echo `pwd`/`ls "$1"`)
	R2FILEPATH=$(echo `pwd`/`ls "$2"`)
fi

#Create log file
printf "`date`\n$USER\n%s\n%s\n\n" \
"$R1FILEPATH" "$R2FILEPATH" \
> "$3"/QualAssessRawSeqs_"$SAMPLE"_results.log ;
wait
echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
echo ' created log output file...'
echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
echo " Quality assessment of $SAMPLE raw sequences completed" 
echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
