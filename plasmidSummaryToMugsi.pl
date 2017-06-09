#!/usr/bin/perl -w

# This script takes in the concatenated output from SSTAR and compares it to a tsv of categories 
# that will go in the MUGSI spreadsheet.  The SSTAR output file must have a blank line in-between
# the different samples and not start with a blank line or headings line.

# USAGE: perl plasmidSummaryToMugsi.pl allArgAnnotResfinderSSTARresults.txt Categories.txt > categoriesResults.tsv

use strict;

open(FILE,$ARGV[0]);
my @file=<FILE>;
close FILE;
open(CATEGORIES,$ARGV[1]);
my $categories=<CATEGORIES>;
close CATEGORIES;
chomp $categories;

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
		my @foundArray;
		my $printOutput="";
		$currentEnd=$x-1;
		my @fields=split("\t",$file[$currentStart]);	
		my @contigs=split('_',$fields[1]);
		if ($contigs[0] eq "MUGSI"){
			print "$contigs[0]_$contigs[1]\t";
		}
		else{
			print "$contigs[0]\t";
		}
		for (my $y=0;$y<@separatedCategories;$y++){	
			my $flag=0;
			for (my $z=$currentStart;$z<=$currentEnd;$z++){
				if ($file[$z] =~ /\Q$separatedCategories[$y]/i){
					my @parts = split("\t",$file[$z]);
					my @gene = split("_", $parts[0]);
					$gene[0].=" ";
					if($gene[0] =~ /\Q$separatedCategories[$y] /i){
						$flag=1;
					}
				 }
			}
			if($flag == 1){
				$printOutput.="$separatedCategories[$y]\t";
				#print "$separatedCategories[$y]\t";
				push @foundArray, $separatedCategories[$y];
			}
			else{
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
