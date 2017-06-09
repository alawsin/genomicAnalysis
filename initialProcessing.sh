#!/bin/sh -l

# This program takes in a path to a directory containing reads, and then trims and removes contaminants in them with  
# BBduk, assembles them with SPAdes, identifies the rna sequences with RNAmmer and then pulls the 16S sequences and  
# runs Kraken and BLAST on them to identify their taxonomy.

# Dependencies: BBduk, SPAdes, RNAmmer, Kraken

# Usage: initialProcessing.sh /path/to/reads 

procs=16 # Number of processors
#mem=24 # Gigabytes of memory - no limit on cluster

kraken_db="$HOME/MiniKrakenDB/"


if [ "${BASH_VERSINFO}" -lt 4 ];
then
	echo "Sorry, you need at least bash-4.0 to run this script." >&2; 
	exit 1; 
fi

DIR="$1"

if [ $# -ne 1 ]
then
        echo "Usage: initialProcessing.sh /path/to/reads" 1>&2
        exit 3

elif [ ! -e "$DIR" ]
then
        echo "Directory does not exist!"
        exit 4

elif [ ! -d "$DIR" ]
then
        echo "Not a Directory"
        exit 5
else
        echo "Path is okay: $DIR"

fi

cd $DIR

for file in *_R1_001.fastq;
do
	filename=$(basename $file _R1_001.fastq);
	echo $filename;

	if [ ! -f ${file}.nt.RemoteBLASTN ]; then
		# ######Trimming Using BBDUK######
		#Removing Contaminants in $HOME/UniVec.fasta#
		echo "Removing Contaminants";
		bbduk.sh -Xmx20g threads=${procs} in=${filename}_R1_001.fastq in2=${filename}_R2_001.fastq out=${filename}-noPhiX-R1.fsq out2=${filename}-noPhiX-R2.fsq ref=$HOME/phiX.fasta k=31 hdist=1
		#Quality and Adapter Trimming in $HOME/adapters.fasta#
		echo "Quality and Adapter Trimming";
		trimmomatic PE -phred33 -threads $procs ${filename}-noPhiX-R1.fsq ${filename}-noPhiX-R2.fsq ${filename}_R1.paired.fq ${filename}_R1.unpaired.fq ${filename}_R2.paired.fq ${filename}_R2.unpaired.fq ILLUMINACLIP:$HOME/adapters.fasta:2:20:10:8:TRUE SLIDINGWINDOW:20:30 LEADING:20 TRAILING:20 MINLEN:50
		# Merge and cleanup
		cat ${filename}_R1.unpaired.fq ${filename}_R2.unpaired.fq > ${filename}.single.fq
		rm *_R[12].unpaired.fq *.fsq

		######Run Kraken on cleaned runs ###### 

		bash $HOME/kraken_krona_fastq_wUnclassified.sh ${filename}_R1.paired.fq ${filename}_R2.paired.fq;

		######Assembling Using SPAdes######
		echo "Assembling Using SPAdes";
		spades.py --only-assembler --careful --pe1-1 ${filename}_R1.paired.fq --pe1-2 ${filename}_R2.paired.fq --pe1-s ${filename}.single.fq -o ${filename}_Assembled --phred-offset 33 -t ${procs} -m 26 --cov-cutoff auto;
		#spades.py --plasmid --only-assembler --careful --pe1-1 ${filename}_R1.paired.fq --pe1-2 ${filename}_R2.paired.fq --pe1-s ${filename}.single.fq -o ${filename}_AssembledPlasmid --phred-offset 33 -t ${procs} -m 26 --cov-cutoff auto;
		
		perl $HOME/removeShortContigs.pl ${filename}_Assembled/scaffolds.fasta

		######Getting 16S RNA Sequences Using RNAmmer######
		echo "Getting 16S RNA Sequences Using RNAmmer";
		rnammer -S bac -m ssu -xml ${filename}_Assembled/scaffolds.fasta.TRIMMED.fasta_rRNA_seqs.xml -gff ${filename}_Assembled/scaffolds.fasta.TRIMMED.fasta_rRNA_seqs.gff -h ${filename}_Assembled/scaffolds.fasta.TRIMMED.fasta_rRNA_seqs.hmmreport -f ${filename}_Assembled/scaffolds.fasta.TRIMMED.fasta_rRNA_seqs.fasta < ${filename}_Assembled/scaffolds.fasta.TRIMMED.fasta;


		#####Blast 16S Against NT######
		#makeblastdb -in $HOME/ntDB/nt -out $HOME/ntDB/nt -dbtype nucl; ##Now using remoteBlast
		blastn -word_size 10 -task blastn -remote -db nt -max_hsps 1 -max_target_seqs 1 -query ${filename}_Assembled/scaffolds.fasta.TRIMMED.fasta_rRNA_seqs.fasta -out ${file}.nt.RemoteBLASTN -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen sscinames";
		kraken --threads $procs --db $kraken_db ${filename}_Assembled/scaffolds.fasta.TRIMMED.fasta_rRNA_seqs.fasta > ${filename}.16s.kraken; 
		kraken-mpa-report --db $kraken_db ${filename}.16s.kraken > ${filename}.16s.kraken.mpa; 
		python ~/Metaphlan_to_krona.py -p ${filename}.16s.kraken.mpa -k ${filename}.16s.kraken.krona 
		ktImportText -o ${filename}.16s.kraken.onlyClassified.html ${filename}.16s.kraken.krona 

	fi

done

######Gather all 16S results
for i in *.RemoteBLASTN; do blastname=$(basename $i _L001_R1_001.fastq.nt.RemoteBLASTN);output=`xargs < $i`;echo -e "${blastname}\t${output}" > ${blastname}.16S.tab;done
cat *.16S.tab > all.16s.tsv
perl ~/16sSummary.pl all.16s.tsv > 16sSummary.tsv

######Summarize all Kraken on cleaned runs like Chris Gulvik's AR Bank Summaries###### 
for F in *kraken; do kraken-report --db ~/MiniKrakenDB $F > "$F".tab; done
for F in *.paired.fq.kraken.tab; do bash ~/summarize_kraken_data/summarize_kraken.sh "$F" >> Summary_reads.kraken.tab; done
for F in *.16s.kraken.tab; do bash ~/summarize_kraken_data/summarize_kraken.sh "$F" >> Summary_16s.kraken.tab; done
sed -i -E 's/_S[0-9]{1,2}_L001_R[1-2]\.paired\.fq\.kraken\.tab//g' Summary_reads.kraken.tab
sed -i -E 's/_S[0-9]{1,2}_L001\.16s\.kraken\.tab//g' Summary_16s.kraken.tab
# sort on genus name, then on CSID
sort -k 7,7 -k 1,1 Summary_reads.kraken.tab > Summary_reads.kraken.sorted.tab
sort -k 7,7 -k 1,1 Summary_16s.kraken.tab > Summary_16s.kraken.sorted.tab
