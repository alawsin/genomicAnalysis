#!/bin/sh 

# This will pull data from MiSeq.  Must have config file
# Usage: miseqRun.sh config.txt

server=`grep "^miseqLocation" $1 | cut -d ":" -f 2-`
localRawDir=`grep "^localRawDirName" $1 | cut -d ":" -f 2-`
date=`grep "^date" $1 | cut -d ":" -f 2-`
user=`grep "^username" $1 | cut -d ":" -f 2-`

echo "Getting raw reads from $server"
echo "Placing raw reads in $HOME/MiSeq_Raw_Data/$localRawDir"
echo "Fastq and Fasta files will be placed in $HOME/Analysis/$date"


#rsync -zarv --prune-empty-dirs --include="*/" --include="*R*_001.fastq.gz" --exclude="*" ${user}@aspen.biotech.cdc.gov:$server $HOME/MiSeq_Raw_Data/$localRawDir

time rsync -aP ${user}@aspen.biotech.cdc.gov:$server $HOME/MiSeq_Raw_Data/$localRawDir

rename 's/I1_001.fastq.gz$/I1_001.fastq.gz.bk/' $HOME/MiSeq_Raw_Data/$localRawDir/Intensities/BaseCalls/*I1_001.fastq.gz
rename 's/I1_001.fastq.gz$/I1_001.fastq.gz.bk/' $HOME/MiSeq_Raw_Data/$localRawDir/Data/Intensities/BaseCalls/*I1_001.fastq.gz

mkdir -p $HOME/Analysis/$date/Fastq/

bash $HOME/scripts/Batch_QualAssessRawSeqs.sh -d $localRawDir -o $HOME/Analysis/$date/Fastq/

mv $HOME/MiSeq_Raw_Data/$localRawDir/Data/Intensities/BaseCalls/*.fastq $HOME/Analysis/$date/Fastq/

cd $HOME/Analysis/$date/Fastq/

rm Undetermined*

bash $HOME/scripts/Batch_QualAssessRawSeqs.sh -p $HOME/Analysis/$date/Fastq/ -o $HOME/Analysis/$date/Fastq/
bash $HOME/scripts/initialProcessing.sh $HOME/Analysis/$date/Fastq/
python $HOME/scripts/extract_TRIMMED_SPAdes.IDBA_genomes_V2.py $HOME/Analysis/$date/Fastq/

mkdir ../Fasta

mv *.fna ../Fasta/

cp 16sSummary.tsv ..
cp Summary_reads.kraken.sorted.tab ..
cp Summary_16s.kraken.sorted.tab ..

cd ../Fasta/

mlst *.fna > allMLST.tsv

for file in *.fna; 
do 
	filename=$(basename $file .fna);
	quast -o ../quast/$filename $file --min-contig 500 --gene-finding --ambiguity-usage one --no-html --no-snps --threads 24;
 	makeblastdb -in $file -out $file -dbtype nucl;
 	blastn -word_size 11 -task blastn -evalue 1e-10 -perc_identity 95 -max_target_seqs 1 -culling_limit 1 -query ~/MCR1gene.fna -db $file -out ${filename}.MCR1.BLASTN \
 	-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen";
 	blastn -word_size 11 -task blastn -evalue 1e-10 -perc_identity 95 -max_target_seqs 1 -culling_limit 1 \
 	-query ~/plasmid_database.fsa -db $file -out ${filename}.PlasmidFinder.BLASTN \
 	-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen";
	blastn -word_size 11 -task blastn -evalue 1e-10 -perc_identity 95 -max_target_seqs 1 -culling_limit 1 \
 	-query ~/all.qac.fasta -db $file -out ${filename}.qac.BLASTN \
 	-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen";
	python ~/c-SSTAR/c-SSTAR100 -s 95 -g $file -d $HOME/ARGannot.r1.fasta > ${file}.ARG-ANNOT.sstar;
	python ~/c-SSTAR/c-SSTAR100 -s 95 -g $file -d $HOME/ResFinder_12-14-2015.srst2.fasta > ${file}.ResFinder.sstar;
	cat ${file}.*.sstar > ${file}.separated.txt;
	python ~/c-SSTAR/c-SSTAR-porin -s 80 -g $file -d $HOME/porins.fasta > ${file}.porins.sstar;

done


#Run SSTAR on all samples, running ARG_ANNOT and Resfinder on all, naming the found gene list <sample>-argRes and translation <sample>-trans###

for file in *.separated.txt;
do
	echo >> $file;
done

for file in *.porins.sstar;
do
	echo >> $file;
done


cat ../quast/*/transposed_report.tsv > ../quast/all.Quast.tsv
sed -i '3~2d' ../quast/all.Quast.tsv
cat *.qac.BLASTN > all.qac.BLASTN
perl ~/scripts/qacSort.pl all.qac.BLASTN > Summary_qac.tsv 
perl ~/scripts/aniAndBuscoFromKraken.pl ../Fastq/Summary_reads.kraken.tab
cat run*.BUSCO/*.BUSCO.tab > Summary_BUSCO.tsv
sed -i '3~2d' Summary_BUSCO.tsv
cat ANI*/ANIm/BestAni*.tab > Summary_ANI.tsv
sed -i '3~2d' Summary_ANI.tsv
rm *.fa
cat *separated.txt > allArgRes.txt
cat *porins.sstar > allPorins.txt
perl $HOME/scripts/plasmidFinderSummary.pl .
perl $HOME/scripts/plasmidSummaryToMugsi.pl summaryPlasmidFinder.blast $HOME/plasmidCategories.txt > plasmidTypes.tsv
perl $HOME/blastnResultsToMugsi.pl allArgRes.txt $HOME/Categories.txt > categories.tsv
perl $HOME/Genes.pl allArgRes.txt $HOME/CPEgenes.txt > cpe.txt
perl $HOME/Genes.pl allArgRes.txt $HOME/betalactamases.txt > betaLacts.txt
perl $HOME/PorinResultsToMugsi.pl allPorins.txt $HOME/porins.txt > porins.tsv
perl $HOME/scripts/sortMMBsheet.pl

#For MUGSI - make sure logLookup.txt in main folder and HOME has MUGSI_MASTER
if [ -s ../logLookup.txt ]; then
	awk -F'\t' '{ print $3 }' ../logLookup.txt > ../aliquotId.txt;
	awk -F'\t' '{ print $4 }' ../logLookup.txt > ../csid.txt;
	cut -f 6  ../logLookup.txt > MALDI.tsv
	awk -F'\t' '{ print $5,"\t",$4,"\t",$2,"\t",$3}' ../logLookup.txt > ../temp1.tsv;
	perl ~/reOrder.pl allMLST.tsv ../aliquotId.txt > allMLST.reOrdered.tsv;
	awk -F'\t' '{ print $2,"\t",$3}' allMLST.reOrdered.tsv > ../temp2.tsv;
	cut -f 3 allMLST.reOrdered.tsv > justMLST.tsv
	perl ~/reOrder.pl Summary_qac.tsv ../aliquotId.txt > Summary_qac.reOrdered.txt
	awk -F'\t' '{ print $2}' Summary_qac.reOrdered.txt > ../tempQac.tsv;
	perl ~/reOrder.pl cpe.txt ../aliquotId.txt > cpe.reOrdered.txt
	awk -F'\t' '{ print $2,"\t","\t"}' cpe.reOrdered.txt > ../temp3.tsv;
	perl ~/reOrder.pl betaLacts.txt ../aliquotId.txt > betaLacts.reOrdered.tsv
	awk -F'\t' '{ print $2}' betaLacts.reOrdered.tsv> ../temp4.tsv;
	perl ~/reOrder.pl categories.tsv ../aliquotId.txt 110 > categories.reOrdered.tsv;
	cut -f 2-111 categories.reOrdered.tsv > ../temp5.tsv;
	#awk '{for(i=2;i<=NF-1;i++) print $i}' categories.reOrdered.tsv > ../temp5.tsv;
	perl ~/reOrder.pl ../16sSummary.tsv ../aliquotId.txt > ../16sSummary.reOrdered.tsv;
	cut -f 2 ../16sSummary.reOrdered.tsv > ../temp6.tsv
	perl ~/reOrder.pl ~/MUGSI_MASTER.txt ../csid.txt > ../temp7.tsv 
	paste ../testOutput.tsv ../csid.txt > testOutput1.tsv
	paste ../temp1.tsv ../temp2.tsv ../tempQac.tsv ../temp3.tsv ../temp4.tsv ../temp5.tsv ../temp6.tsv justMLST.tsv ../aliquotId.txt ../temp7.tsv > ../mugsi1.output.tsv

	####Paste 1st part into MUGSI spreadsheet - different column lengths####
	####Check SHV and TEM and any new betalactams####

	perl ~/reOrder.pl plasmidTypes.tsv ../aliquotId.txt 35 > plasmidTypes.reOrdered.tsv
	cut -f 2-36  plasmidTypes.reOrdered.tsv > ../temp8.tsv;
	perl ~/reOrder.pl porins.tsv ../aliquotId.txt 10 > porins.reOrdered.tsv
	cut -f 2- porins.reOrdered.tsv > ../temp9.tsv;
	paste MALDI.tsv ../csid.txt ../aliquotId.txt ../temp8.tsv MALDI.tsv ../csid.txt ../aliquotId.txt ../temp9.tsv > ../mugsi2.output.tsv
	rm ../temp[1-9].tsv

	####Paste 2nd part into MUGSI spreadsheet - check porins - possible repeats push into extra columns####


#If not MUGSI, and has aliquotId.txt, to just reOrder#
else
	if [ -s ../aliquotId.txt ] || [ -s ../aliquotID.txt ]; then
		if [ -s ../aliquotId.txt ]; then
			aliquotFile="../aliquotId.txt"
		else aliquotFile="../aliquotId.txt"
		fi
			perl ~/reOrder.pl allMLST.tsv $aliquotFile > ../allMLST.reOrdered.tsv;
			perl ~/reOrder.pl cpe.txt $aliquotFile > ../cpe.reOrdered.txt
			perl ~/reOrder.pl betaLacts.txt $aliquotFile > ../betaLacts.reOrdered.tsv
			perl ~/reOrder.pl categories.tsv $aliquotFile 110 > ../categories.reOrdered.tsv;
			perl ~/reOrder.pl plasmidTypes.tsv $aliquotFile 35 > ../plasmidTypes.reOrdered.tsv
			perl ~/reOrder.pl porins.tsv $aliquotFile 10 > ../porins.reOrdered.tsv
			perl ~/reOrder.pl ../Fastq/16sSummary.tsv $aliquotFile > ../16sSummary.reOrdered.tsv;
			perl ~/reOrder.pl ../Fastq/Summary_16s.kraken.sorted.tab $aliquotFile > ../Summary_16s.kraken.reOrdered.tsv
			perl ~/reOrder.pl ../Fastq/Summary_reads.kraken.sorted.tab $aliquotFile > ../Summary_reads.kraken.reOrdered.tsv
			perl ~/reOrder.pl ../MMBoutput.tsv $aliquotFile > ../MMBoutput.tsv.reOrdered.tsv
		
	fi
	if [ -s ../csid.txt ]; then
		perl ~/reOrder.pl ~/MUGSI_MASTER.txt ../csid.txt > ../temp7.tsv 
	fi
fi
