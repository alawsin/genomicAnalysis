#!/bin/bash
#  @/Author: Christopher A. Gulvik, Ph.D.
#  +/Version 1.1
#  +/6 February 2015
#  Dependencies:  awk, cat, chmod, cut, grep, ls, paste, rm, sed, tail, tr, uname, uniq, wc
#    script must also be present $HOME/scripts/QualAssessRawSeqs.sh

function usage() { 
	echo "
	USAGE: $0 -d InputDirectory || -p /Input/Path -o /Output/Directory [-T threads] [-m email@host.com]
	The raw reads must be in Illumina 1.8+ FastQ format.
    
	InputDirectory must be formatted:
	$HOME/MiSeq_Raw_Data/<directory>/Data/Intensities/BaseCalls/, where
	<directory> is the name of the run's directory
	(e.g., 140514_M02901_0015_000000000-ABEED).
	
	InputPath must be the complete filepath to the .fastq.gz files
	(e.g., $HOME/MiSeq_Data/project/raw/)

	OutputDirectory must also be a complete filepath.
	"
	}

#Script dependency check
QASCRIPT_DEPEND="$HOME/scripts/QualAssessRawSeqs.sh"
if [[ -f "$QASCRIPT_DEPEND" ]]; then
	echo "$QASCRIPT_DEPEND located..."
	if [[ -x "$QASCRIPT_DEPEND" ]]; then
		echo "$QASCRIPT_DEPEND is executable..."
	else
		echo "$QASCRIPT_DEPEND is not executable by $USER..." >&2
		exit 1
	fi
else
	echo "$SCRIPT_DEPEND absent"
	echo 'This script depends on $HOME/scripts/QualAssessRawSeqs.sh' >&2
	exit 1
fi

#Require arguments
if [[ $# -lt 4  || $# -gt 9 ]]; then
	echo "ERROR: incorrect number ("$#") of arguments were provided." >&2
	usage
	exit 1
fi

#Replaced getopts while loop with for statement
nopts=$#
for ((i=1 ; i <= nopts ; i++)); do
	case "$1" in
		-d | --in-dir)
			RUNFOLDER="$2"
			INDATADIR="$HOME"/MiSeq_Raw_Data/"$RUNFOLDER"/Data/Intensities/BaseCalls
			shift 2
			;;
		-p | --in-path)
			PPATH="$2"
			INDATADIR="$PPATH"
			shift 2
			;;
		-m | --mailto)
			EMAIL="$2"
			echo "${BOLD}-m${NORM} argument invoked:  $2"
			shift 2
			;;
		-o | --out-path)
			OUTDATADIR="$2"
			shift 2
			;;
		-T | --threads)
			NUMPROC="$2"
			echo "${BOLD}-T${BOLD} argument invoked:  $2"
			shift 2
			;;
		-v | --verbose)
			VERBOSE=1
			echo "${BOLD}verbose${NORM} mode invoked"
			shift
			;;
		-h | --help)
			usage
			exit 1
			;;
		\?)
			echo "ERROR: ${BOLD}$2${NORM}is not a valid argument" >&2
			usage
			exit 1
			;;
	esac
done

#I/O handling
if [[ -z "$RUNFOLDER$PPATH" ]]; then  #require at least one arg exists
    echo 'ERROR: Option -d or -p is required and needs an argument.' >&2
    exit 1
fi

if [ "${INDATADIR:0:1}" != "/" ] && [ "${INDATADIR:0:1}" != "$" ]; then
	echo 'ERROR: The full path was not specfied.' >&2
	exit 1
fi

if [[ -z "$OUTDATADIR" ]]; then  #require
    echo "ERROR: Option -o is required and needs an argument." >&2
    exit 1
fi

if [ ! -d "$OUTDATADIR" ]; then  #create outdir if absent
	mkdir -p "$OUTDATADIR"
fi

echo "Input Sequence Data directory path:  $INDATADIR"
echo "Output Data Analysis directory path:  $OUTDATADIR"

#Determine cluster or desktop environment
OPSYS=$(uname -s)
if [ "$OPSYS" != Linux ] && [ "$OPSYS" != Darwin ]; then
	echo 'ERROR: Your operating system does not appear to be supported.' >&2
	exit 1
elif [[ NSLOTS -gt 0 ]]; then  #cluster
	NUMPROC="$NSLOTS"
elif [[ "$NUMPROC" != "" ]]; then  #noncluster, threads arg used
	echo "$NUMPROC processors will be used..."
elif [[ "$NUMPROC" == "" ]]; then  #noncluster, threads arg unused
	#OS and CPU information
	NUMPROC=1  #Default
	if [[ "$OPSYS" == Linux ]]; then
		NUMPROC=$(grep -c ^processor /proc/cpuinfo)
		if [[ "$NUMPROC" -gt 2 ]]; then
			echo "found $NUMPROC processors..."
			NUMPROC="$(($NUMPROC - 1))"  #leave 1 process open
		fi
	elif [[ "$OPSYS" == Darwin ]]; then
		NUMPROC=$(sysctl hw.physicalcpu | awk '{print $2}')
		if [[ "$NUMPROC" -gt 2 ]]; then
			echo "found $NUMPROC processors..."
			NUMPROC="$(($NUMPROC - 1))"  #leave 1 process open
		fi
	fi
else
	echo 'ERROR' >&2
	exit 1
fi

#Test barcodes/indices were entered correctly on the machine
# if either undetermined.fastq.gz file is > any other sample.fastq.gz file, report a warning
#for undfile in Undetermined_*.fastq.gz; do
#	UNDT1SZ=$(du -b Undetermined_S0_L001_R1_001.fastq.gz | cut -f 1)
#	UNDT2SZ=$(du -b Undetermined_S0_L001_R2_001.fastq.gz | cut -f 1)
#verify two smallest files in directory are 'Undetermined' reads	
UNDTESTCOUNT=$(ls "$INDATADIR"/*.fastq.gz | grep -o 'Undetermined' | wc -l)  #should return 2 if undetermined files R1 and R2 exist
UNDTESTSTR=$(ls -S "$INDATADIR"/*.fastq.gz | tail -n 2 | awk -F / '{print $(NF +0)}' | cut -d _ -f 1 | tr '\n' '_')  #concats two smallest sample names with underscore
if [[ UNDTESTCOUNT -gt 0 ]]; then  #at least 1 Undetermined file found
	if [[ UNDTESTSTR == Undetermined_Undetermined_\n ]]; then
		echo 'Undetermined_*.fastq.gz files are the smallest two sequence files'
		echo 'for this run. This is usually an indication that barcodes were'
		echo 'properly entered into the MiSeq.'
	elif [[ UNDTESTSTR != Undetermined_Undetermined_ ]]; then
		echo 'WARNING: undetermined file(s) larger than at least one sample file'
		echo '         Consider verifying barcodes were properly entered.'
	else
		echo 'WARNING: Unable to determine whether Undetermined_*.fastq.gz'
		echo '         files are the smallest two sequence files for this run.'
	fi
else
	echo 'No Undetermined_*.fastq.gz sequence files were found in this run.'
fi

#Test if even number of sequence files (for paired end reads)
SEQFILES=$(ls "$INDATADIR"/*.fastq.gz | uniq -u | wc -l)
if [[ $(($SEQFILES % 2)) -ne 0 ]]; then  #value is odd
	echo "ERROR: uneven number ($SEQFILES) of sequence files" >&2
	exit 1
fi
echo "$SEQFILES *.fastq.gz files will be processed..."

PAIREDFILES="$(($SEQFILES / 2))"
if [[ PAIREDFILES -lt NUMPROC ]]; then
	NUMPROC="$PAIREDFILES"
	echo 'Number of paired files were less than number of processes (-T) requested...'
fi
echo "using $NUMPROC processors..."

#Generate bash script to perform QA calculations
ls "$INDATADIR"/*.fastq.gz | uniq -u | paste -d " " - - | awk '{print "$HOME/scripts/QualAssessRawSeqs.sh " $0}' | awk -v awk_var="${OUTDATADIR} &" '{print $0,awk_var}' | awk -v cpuvar="${NUMPROC}" '1; NR%cpuvar==0{print "wait"}' | awk 'BEGIN{print "#!/bin/bash"}{print}' > "$OUTDATADIR"/Pairs.sh
echo 'wait' >> "$OUTDATADIR"/Pairs.sh
chmod u+x "$OUTDATADIR"/Pairs.sh
bash "$OUTDATADIR"/Pairs.sh
wait

#Summarize and clean up individual data and log files
echo 'Summarizing and cleaning up files...'
cat "$OUTDATADIR"/QualAssessRawSeqs*_results.tsv > "$OUTDATADIR"/QualAssessRawSeqs_Sum.tsv
if [[ "$OPSYS" == Linux ]]; then
	sed -n '3~2!p' "$OUTDATADIR"/QualAssessRawSeqs_Sum.tsv > "$OUTDATADIR"/QualAssessRawSeqs_Summary.tsv
elif [[ "$OPSYS" == Darwin ]]; then
	gsed -n '3~2!p' "$OUTDATADIR"/QualAssessRawSeqs_Sum.tsv > "$OUTDATADIR"/QualAssessRawSeqs_Summary.tsv
fi
rm "$OUTDATADIR"/QualAssessRawSeqs_Sum.tsv "$OUTDATADIR"/QualAssessRawSeqs*_results.tsv
cat "$OUTDATADIR"/QualAssessRawSeqs*_results.log > "$OUTDATADIR"/QualAssessRawSeqs_Summary.log
echo -e '\n\nThe script to generate these data was' >> "$OUTDATADIR"/QualAssessRawSeqs_Summary.log
cat "$OUTDATADIR"/Pairs.sh >> "$OUTDATADIR"/QualAssessRawSeqs_Summary.log
rm "$OUTDATADIR"/Pairs.sh
echo -e '\n' >> "$OUTDATADIR"/QualAssessRawSeqs_Summary.log 
rm "$OUTDATADIR"/QualAssessRawSeqs*_results.log
echo 'Batch quality assessment of raw sequences has completed.'

#Email data
#use mail or mailx to send a plain text file
#mail -s "Subject: QA of Raw Reads (attached)" "$EMAIL" < "$OUTDATADIR"/QualAssessRawSeqs_Summary.tsv
#echo "Data were just sent to $EMAIL."
