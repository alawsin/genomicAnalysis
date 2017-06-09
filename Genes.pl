#!/usr/bin/perl -w

# This script takes in the concatenated output from SSTAR and compares it to a tsv of categories 
# that will go in the MUGSI spreadsheet.  The SSTAR output file must have a blank line in-between
# the different samples and not start with a blank line or headings line.

# USAGE: perl Genes.pl allArgAnnotResfinderSSTARresults.txt geneList.txt > foundGeneResults.tsv

use strict;
use List::MoreUtils qw(uniq);

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
		my @genes;
		$currentEnd=$x-1;
		my @fields=split("\t",$file[$currentStart]);	
		my @contigs=split('_',$fields[2]);
		print "$contigs[0]\t";
		my $flag=0;
		for (my $y=0;$y<@separatedCategories;$y++){	
			for (my $z=$currentStart;$z<$currentEnd;$z++){
				if($file[$z] =~ /\Q$separatedCategories[$y]/i){
					my @parts = split("\t",$file[$z]);
					if ($parts[1] !~ /TR$/){
						$parts[1] =~ s/^bla//;
						my $idLimit;
						chop $parts[3];
						if ($parts[1] =~ /$separatedCategories[$y](_|\-|\*|\?|T|$)/i){
							my @alleleParts = split ("_", $parts[1]);
							my $variant = $alleleParts[0];
							if (($parts[1] =~ /^TEM/i) || ($parts[1] =~ /^SHV/i)){
								$idLimit = 95; # 95% ID
								if ($parts[3] >= $idLimit){
									#my $lengthLimit = .95 * $parts[5]; # 95% Coverage
									my $lengthLimit = $parts[5]; # 100% Coverage
									if ($parts[4] >= $lengthLimit){
										$flag=1;
										#push @genes, $separatedCategories[$y];
										push @genes, $variant;
									}
								}
							}
							else{ # if not a TEM or SHV, must be 100%ID
								$idLimit = 100;
								if ($parts[3] == $idLimit){
									if ($parts[4] == $parts[5]){
										$flag=1;
										#push @genes, $separatedCategories[$y];
										push @genes, $variant;
									}
								}
							}
						}
					}					
				}
			}
		}
		my @unique_genes = uniq @genes;
		if($flag == 1){
			print "Yes - ", "@unique_genes";
		}
		else{
			print "No";
		}
		print "\n";
		$currentStart=$x+1;
	}
}
