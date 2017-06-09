#!/usr/bin/perl -w

# This script takes in a kraken output summary and gets the top genus for every sample.
# With the genus identified, ANI is ran on every sample against a collection of RefSeq
# genomes for every genus.  Also, depending on the genus, BUSCO is ran with the most 
# appropriate Bacteria dataset.

# Assumes kraken output summary is in Run/Fastq folder and fasta samples are in Run/Fasta
# Assumes ANI db is at ~/aniDB 
# Assumes BUSCO lineages are at ~/busco/lineages

# Dependencies: bbduk, prokka, BUSCO, pyaniANI.py (*executable* and in PATH)

# Usage: perl aniAndBuscoFromKraken.pl Summary_reads.kraken.tab

use File::Basename;
use File::Copy;
use Cwd;
use strict;

my $aniDB="~/aniDB";
my $buscoDB="~/busco/lineages";

my $dirname=dirname($ARGV[0]);
print "Directory of kraken summary is $dirname\n";
chdir $dirname;
chdir "../Fasta";
my $currentDir=cwd();
print "Current directory is $currentDir\n";
opendir(CDH, $currentDir) || die "Can't opendir $currentDir: $!";
my @fileList=grep {/\.fna$/} readdir(CDH);
closedir CDH;
#print $_,"\n" for (@fileList);

open(FILE,$ARGV[0]);
my @file=<FILE>;

my %buscoHash=(
	"Klebsiella"=>"enterobacteriales_odb9",
	"Escherichia"=>"enterobacteriales_odb9",
	"Acinetobacter"=>"gammaproteobacteria_odb9",
	"Burkholderia"=>"betaproteobacteria_odb9",
	"Enterobacter"=>"betaproteobacteria_odb9",
	"Enterococcus"=>"lactobacillales_odb9",
	"Flavobacterium"=>"bacteroidetes_odb9",
	"Morganella"=>"enterobacteriales_odb9",
	"Mycobacterium"=>"actinobacteria_odb9",
	"Pseudomonas"=>"gammaproteobacteria_odb9",
	"Staphylococcus"=>"bacillales_odb9",
	"Serratia"=>"enterobacteriales_odb9",
	"Proteus"=>"enterobacteriales_odb9",
);


my $aniBase=glob($aniDB);
my $buscoBase=glob($buscoDB);
#mkdir "ANI";
for(my $count=0;$count<@file;$count++){
	if ($file[$count] !~ /^\s*$/){
		my @fields=split("\t",$file[$count]);
		my $sample=$fields[0];
		my $genus=$fields[6];
		my @file=grep {/^$sample/} @fileList;

		###BUSCO RUN###
		print "BUSCO\t$sample\t$genus\t$file[0]\n";
		print "${buscoBase}/$buscoHash{$genus}\n";
		my $buscoDir=$buscoBase."/".$buscoHash{$genus};
		if (-e $buscoDir and -d $buscoDir) {
			`rename.sh in=$file[0] out=$file[0].fa prefix=contig`;
			my $prokkaDir=$sample.".prokka";
			`prokka --outdir $prokkaDir --prefix $sample $file[0].fa`;
			`BUSCO.py -i ${prokkaDir}/${sample}.faa -l ${buscoDB}/$buscoHash{$genus} -o ${sample}.BUSCO -m prot`;
			open(BPRINT,">run_${sample}.BUSCO/$sample.BUSCO.tab");
			open(BOUTPUT, "run_${sample}.BUSCO/short_summary_${sample}.BUSCO.txt");
			my @bfile=<BOUTPUT>;
			my $bnumerator;
			my $bdenominator;
			for (@bfile){
				chomp $_;
				if ($_ =~ /Complete BUSCOs/){
					my @bfields=split("\t",$_);
					$bnumerator=$bfields[1];
				}
				if ($_ =~ /Total BUSCO groups searched/){
					my @bfields=split("\t",$_);
					$bdenominator=$bfields[1];
				}
			}
			my $bpercent=${bnumerator}/${bdenominator}*100;
			print BPRINT "Sample\tBUSCO Score\tBUSCO Percentage\n";
			print BPRINT "$sample\t${bnumerator}/${bdenominator}\t${bpercent}%\n";
			close BPRINT;
		}
		else{
			print "BUSCO directory does not exist\n";
		}

		###ANI RUN###
		print "ANI\t$sample\t$genus\t$file[0]\n";
		my $aniDir=$aniBase."/".$genus;
		print "$aniDir\n";
		if (-e $aniDir and -d $aniDir) {
			opendir(ADH, $aniDir)|| die "Can't opendir $aniDir: $!";
			my @refList=grep {/\.fna$/} readdir(ADH);
			mkdir "ANI-$file[0]-$genus";
			print "ANI-$file[0]-$genus\n";
			copy ($file[0], "ANI-$file[0]-$genus");
			print "$file[0]\n";
			for (@refList){
				copy ("${aniDir}/$_",  "ANI-$file[0]-$genus");
			}
			`average_nucleotide_identity.py -i "ANI-$file[0]-$genus" -o "ANI-$file[0]-$genus/ANIm" -m ANIm -g`;
			open(ANIRESULT, "ANI-$file[0]-${genus}/ANIm/ANIm_percentage_identity.tab");
			my @aniResults=<ANIRESULT>;
			close ANIRESULT;
			my $max=-1;
			my $maxCol=1;
			my $heading=$aniResults[0];
			chomp $heading;
			my @headingFields=split("\t",$heading);
			for (my $x=1;$x<@aniResults;$x++){
				if ($aniResults[$x] =~ /^${sample}/){
					chomp $aniResults[$x];
					my @aniFields=split("\t",$aniResults[$x]);
					for (my $y=1;$y<@aniFields;$y++){
						if (($aniFields[$y]>$max) && ($headingFields[$y] !~ /^${sample}/)){
							$max=$aniFields[$y];
							$maxCol=$y;
						}
					}
				}
			}
			open(BESTHIT,"${aniDir}/$headingFields[$maxCol].fna");
			my $firstline = <BESTHIT>;
			chomp $firstline;
			$firstline =~ s/>//;
			my @firstFields=split(" ",$firstline);
			my $species="$firstFields[1] $firstFields[2]";
			open(BESTANI, ">ANI-$file[0]-${genus}/ANIm/BestAni-${sample}.tab");
			print BESTANI "Sample\tBest Match\tANIm Score\n";
			print BESTANI "${sample}\t$species($headingFields[$maxCol].fna)\t$max\n";
			close BESTANI;
			close BESTHIT;
		}
	}
	else{
		next;
	}
}

close CDH;
close FILE;