#!/usr/bin/env python

# python ~/scripts/filter_contigs.py 1000 INPUT.fna
# for FILE in *.fna; do python ~/scripts/filter_contigs.py 1000 $FILE; done

import sys
from Bio import SeqIO
from Bio.SeqIO import FastaIO

min_length, fasta_file_path = sys.argv[1:]

def filtered_contigs_generator(min):
	for contig in SeqIO.parse(input_fasta, 'fasta'):
		if len(contig) >= min:
			yield contig

#with open(fasta_file_path.replace('fna', 'fa'.format(min_length)), 'w') as filtered_fasta:
with open(fasta_file_path.replace('fna', 'filtered{}.fa'.format(min_length)), 'w') as filtered_fasta:
	with open(fasta_file_path, 'rU') as input_fasta:
		SeqIO.write(filtered_contigs_generator(int(min_length)), filtered_fasta, 'fasta')
