# Randomizes bases in each entry of an input FASTA
# Usage:
# cat in.fa | python scambleFasta.py > out.fa

import random
import sys
from Bio import SeqIO

seqs = SeqIO.parse(sys.stdin, 'fasta')

i = 1
for entry in seqs:
	seq = list(entry.seq)
	random.shuffle(seq)
	print('>negative_%s' % i)
	print(''.join(seq))
	i += 1
