#!/usr/bin/perl -w

# This script takes in the concatenated output from SSTAR and compares it to a tsv of categories 
# that will go in the MUGSI spreadsheet.  The SSTAR output file must have a blank line in-between
# the different samples and not start with a blank line or headings line.

# USAGE: perl blastnResultsToMugsi.pl allPorinSSTARresults.txt Categories.txt > categoriesResults.tsv

use strict;

open(FILE,$ARGV[0]);
my @file=<FILE>;
close FILE;
open(CATEGORIES,$ARGV[1]);
my $categories=<CATEGORIES>;
close CATEGORIES;
chomp $categories;
my $flag;
my @separatedCategories=split("\t",$categories);
my $currentStart=0;
my $currentEnd;
for (my $x=0;$x<@file;$x++){
	if ($file[$x] !~ /^\s*$/){
		next;
	}
	else{
		if (($file[$x] =~ /^\s*$/) && ($file[$x-1] =~ /^\s*$/)){
			print "\n";
			$currentStart=$x+1;
			next;
		}
		$currentEnd=$x-1;
		my @foundArray;
		my $printOutput="";
		my @fields=split("\t",$file[$currentStart]);	
		my @contigs=split('_',$fields[2]);
		if ($contigs[0] eq "MUGSI"){
			print "$contigs[0]_$contigs[1]\t";
		}
		else{
			print "$contigs[0]\t";
		}
		for (my $y=0;$y<@separatedCategories;$y++){	
			my $flag=0;
			my $truncation=0;
			for (my $z=$currentStart;$z<=$currentEnd;$z++){
				if ($file[$z] =~ /\Q$separatedCategories[$y]/i){
					my @parts = split("\t",$file[$z]);
					if (($parts[1] =~ /TR$/) && ($parts[1] =~ /^omp/i)){
						chop $parts[3];
						my $idLimit;
						$idLimit = 80;
						if ($parts[3] >= $idLimit){
							if ($parts[4] >= $parts[5]){
									$flag=1;
									push @foundArray, $separatedCategories[$y];
									$printOutput.="$parts[1]\t";
									#print "$parts[1]";
									#print "\t";
							}
						}
					}
				}
			}
			if($flag != 1){
				$printOutput.="\t";
				#print "\t";
			}

		}
		print join(', ', @foundArray);
		print "\t",$printOutput;
		$currentStart=$x+1;
		print "\n";
	}
}
