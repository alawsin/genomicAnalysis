#!/usr/bin/env python
#  @ Author: Christopher A. Gulvik, Ph.D.
#  + Version 1.1
#  + 18 May 2015
#  Dependencies:  none


from __future__ import print_function
import os
import re
import string
import sys


if sys.argv[1].startswith('/'):
	print ('Finding all contig.fa and scaffolds.fasta files in: ' + (sys.argv[1]))
elif sys.argv[1].startswith('$'):
	print ('Finding all contig.fa and scaffolds.fasta files in: ' + (sys.argv[1]))
elif sys.argv[1].startswith('~'):
	print ('Finding all contig.fa and scaffolds.fasta files in: ' + (sys.argv[1]))
else:
	print ('Usage: python script.py /full/filepath')
	sys.exit (1)


basepath = ''
for root, dirs, files in os.walk(sys.argv[1], topdown=False):  #recursive searching
	path = os.path.relpath(root, basepath).split('/')
	#print ((len(path) - 1) * '~~~' , os.path.basename(root))
	for file in files:
	 	if file == 'contig.fa' or file == 'scaffolds.fasta.TRIMMED.fasta':  #'contig.fa' for IDBA and 'scaffolds.fasta' for SPAdes
			#print (len(path) * '==>' , file)
			#print ('found ' + os.path.join(root, file))
			filepath, filename = os.path.split(os.path.abspath(os.path.join(root, file)))
			spades_kmerdir = re.compile(r"K[0-9]{2}")
			if not spades_kmerdir.search(str(filepath)):  #SPAdes puts 'scaffolds.fasta' also in highest kmer directory; ignore this
				parentdir = os.path.basename(filepath)
				newheader = os.path.splitext(parentdir)[0]
				fasta = open(os.path.join(root, file))
				newfasta = open(str(filepath) + '.fna', 'w')
				i = 1
				for line in fasta:
					line = re.sub('\n', '', line)
					if line.startswith('>'):
						newfasta.write('>' + newheader + '_' + str(i) + '\n')
						i += 1
					else:
						newfasta.write(line + '\n')
				print (' renamed ' + os.path.join(root, file))
				print ('  to ' + os.path.join(str(filepath) + '.fna') + '\n')
				fasta.close()
				newfasta.close()

##AWK to count and add underscore to each header line:  awk '/^>/{$0=$0"_"(++i)}1' in.fa > out.fa
