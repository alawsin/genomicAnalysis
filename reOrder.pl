#!/usr/bin/perl -w

use strict;

open(FILE,$ARGV[0]);
my @file=<FILE>;
close FILE;
open(LIST,$ARGV[1]);
my @list=<LIST>;
close LIST;
my $numTabs=0;
if (defined $ARGV[2]){
	$numTabs=$ARGV[2];
}
else{
	$numTabs=0;
}


for (my $x=0;$x<@list;$x++){
	my $flag=0;
	chomp $list[$x];
	for (my $y=0;$y<@file;$y++){
		if ($file[$y] =~ m/$list[$x]/i){
			$flag=1;
			print $file[$y];
		}
	}
	if ($flag == 0){
		print "$list[$x]";
		if ($numTabs != 0){
			for(my $count=0;$count<$numTabs;$count++){
				print "\t";
			}
		}
		print "\n";
	}

}
