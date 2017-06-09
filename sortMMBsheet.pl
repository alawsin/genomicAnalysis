#!/usr/bin/perl -w

# Assumes Kraken and QC summary is in Run/Fastq folder and all other necessary files
# are in the Run/Fasta folder.


use strict;

open(OUTPUT,">../MMBoutput.tsv");
my %outputHash;
my %q30Hash;

###Kraken ID###
open(KRAKEN, "../Fastq/Summary_reads.kraken.tab");
my @krakenArray=<KRAKEN>;

foreach (@krakenArray){
	my @krakenFields=split("\t",$_);
	$outputHash{$krakenFields[0]}=$krakenFields[15];
}
close KRAKEN;

###QC Scores###
open(QC, "../Fastq/QualAssessRawSeqs_Summary.tsv");
my @qcArray=<QC>;
foreach (%outputHash){
	for (my $x=1;$x<@qcArray;$x++){
		if($qcArray[$x] =~ /^$_/){
			my @qcFields=split("\t",$qcArray[$x]);
			$q30Hash{$_}=$qcFields[2]; # Needed to calculate Estimated coverage
			for (my $y=1;$y<=12;$y++){
				$outputHash{$_}.="\t$qcFields[$y]";
			}
		}
	}
}
close QC;

###Quast Scores###
open(QUAST, "../quast/all.Quast.tsv");
my @quastArray=<QUAST>;
foreach (%outputHash){
	for (my $x=1;$x<@quastArray;$x++){
		if($quastArray[$x] =~ /^$_/){
			my @quastFields=split("\t",$quastArray[$x]);
			my $coverage = $q30Hash{$_} / $quastFields[15]; # Coverage = Q30 total / Total length
			$outputHash{$_}.="\t$coverage\t$quastFields[13]\t$quastFields[15]";
		}
	}
}
close QUAST;

###BUSCO Scores###
open(BUSCO, "Summary_BUSCO.tsv");
my @buscoArray=<BUSCO>;
foreach (%outputHash){
	for (my $x=1;$x<@buscoArray;$x++){
		if($buscoArray[$x] =~ /^$_/){
			my @buscoFields=split("\t",$buscoArray[$x]);
			$outputHash{$_}.="\t$buscoFields[1]";
		}
	}
}
close QUAST;

###ANI Scores###
open(ANI, "Summary_ANI.tsv");
my @aniArray=<ANI>;
foreach (%outputHash){
	for (my $x=1;$x<@aniArray;$x++){
		if($aniArray[$x] =~ /^$_/){
			chomp $aniArray[$x];
			my @aniFields=split("\t",$aniArray[$x]);
			my $aniScore=int($aniFields[2]*10000)/100;
			my $aniMatch=$aniFields[1];
			$outputHash{$_}.="\t$aniScore%\-$aniMatch";
		}
	}
}
close ANI;

###Print OUTPUT###
print OUTPUT "Sample\tKRAKEN ID\tQ20_Total_[bp]\tQ30_Total_[bp]\tQ20_R1_[bp]\tQ20_R2_[bp]\t";
print OUTPUT "Q20_R1_[%]\tQ20_R2_[%]\tQ30_R1_[bp]\tQ30_R2_[bp]\tQ30_R1_[%]\tQ30_R2_[%]\t";
print OUTPUT "Total_Sequenced_[bp]\tTotal_Sequenced_[reads]\tEstimated coverage\t";
print OUTPUT "# contigs\tCumulative_Length_Assembly (bp)\tBUSCO\tANI\n";
foreach (sort keys %outputHash){
	print OUTPUT "$_\t$outputHash{$_}\n";
}
close OUTPUT;