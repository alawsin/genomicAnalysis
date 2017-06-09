#!usr/bin/perl -w

# This program will take all the plasmidFinder Blast output files in a given directory
# and print the filename, bacteria, and ST-type.

use strict; 

opendir (DH, $ARGV[0]);
my @allFiles = readdir (DH);
my @sortedFiles = sort @allFiles;
closedir(DH);

open (OUTPUT, ">$ARGV[0]/summaryPlasmidFinder.blast");
foreach my $file (@sortedFiles){
	if ($file =~ /PlasmidFinder.BLASTN/){
		#print $file;
		open (EACHBLAST, "$ARGV[0]/$file");
		while (my $line = <EACHBLAST>){
			chomp $line;
			my @fields = split("\t", $line); 
			if (($fields[2] >= 95.00) && ($fields[3] == $fields[12])){ # At least 95% ID & 100% Coverage
				print OUTPUT "$line\n";
				#print OUTPUT "$file\t$line\n";
				#print "$fields[0]\t$fields[2]% ID\t$fields[3] Matching Length\t$fields[12] Query Length\n";
			}
		}
		print OUTPUT "\n";
	}
	
}
close OUTPUT;
