#! usr/bin/perl -w

use strict;

open (INPUT, $ARGV[0]);
open (OUTPUT, ">$ARGV[0].TRIMMED.fasta");

my $printFlag = 0;
while (my $line=<INPUT>){
	if ($line =~ /^>/){
		my @fields = split('_',$line);
		if (($fields[3] > 500) && ($fields[5] >= 2)){
			$printFlag = 1;
		}
		else{
			$printFlag = 0;
		}
	}
	if ($printFlag == 1){
		print OUTPUT $line;
	}
}
close INPUT;
close OUTPUT;
