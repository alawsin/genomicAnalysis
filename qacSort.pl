#!/usr/bin/perl -w

# Sort qac blast output.  Only print out 100% coverage hits per sample.

# Usage: perl qacSort.pl qacBlast.txt > qacBlastSorted.tsv

use strict;

my %qacHash=(
	"NZ_CP006659.2:2298981-2299820"=>"qacDeltaD",
	"NC_001735.4:33886-34229"=>"qacE",
	"NC_019081.1:5454-5786"=>"qacH",
	"NC_019286.1:3347-3679"=>"qacF",
	"NG_052493.1"=>"qacE",
	"qacI"=>"qacI",
	"qacL"=>"qacL",
);
open (FILE, $ARGV[0]);
my @file=<FILE>;
my $name;
my $current;
my $currentHits="";
for (my $x=0;$x<@file;$x++){
	if ($file[$x] !~ /^\s*$/){
		my @fields=split("\t",$file[$x]);
		my @nameFields=split("_",$fields[1]);
		if ($nameFields[0] =~ /^MUGSI/){
			$name="$nameFields[0]_$nameFields[1]";
		}
		else{
			$name=$nameFields[0];
		}
		if ($name ne $current){
			if ($currentHits !~ /^\s*$/){
				print "$current\t$currentHits\n";
				$currentHits="";
			}
			$current = $name;
		}
		if ($fields[3] == $fields[12]){
			if ($currentHits !~ /^\s*$/){
				$currentHits.=",";
			}
			$currentHits.="$qacHash{$fields[0]}($fields[2])";
		}
	
	}
}
if ($currentHits !~ /^\s*$/){
				print "$name\t$currentHits\n";
}