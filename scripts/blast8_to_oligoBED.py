#!/usr/bin/env python

###############################################################
# Reformat BLAST8 alignments into BED coordinates of fixed size
# Author: Colin Shew
# Updated: 14 Oct 2020
###############################################################

import sys
import argparse
import pandas as pd

# parse arguments
parser = argparse.ArgumentParser(description='Reformat BLAST8 alignments into BED coordinates of fixed size')

parser.add_argument('blat', metavar='blat', type=str, help='BLAT output (BLAST8 format, single header line)')
parser.add_argument('minL', metavar='min', type=int, help='Min alignment length to keep')
parser.add_argument('maxL', metavar='max', type=int, help='Max alignment length to keep')
parser.add_argument('k', metavar='length', type=int, help='Output oligo length')
parser.add_argument('genome', metavar='genome', type=str, help='Genome name')

args = parser.parse_args()

def printBED(quer, subj, sStart, sEnd, genome, strand):
	print('\t'.join([subj, str(sStart), str(sEnd), '%s_%s.%s:%s-%s' % (quer, genome, subj, sStart, sEnd), '0', strand]))

# read BLAT output
for row in pd.read_table(args.blat, header=0).iterrows():
	quer, subj, percID, alignL, mismatch, gapOpen, qStart, qEnd, sStart, sEnd, eVal, bitScore = row[1]
	strand = '+'
	# flip back coords on minus strand
	if int(sStart) > int(sEnd):
		sStart, sEnd = sEnd, sStart
		strand = '-'
	# change to 0-based
	sStart = int(sStart-1)
	# skip if range too large
	if not args.minL <= int(sEnd)-int(sStart) <= args.maxL:
		continue
	# adjust sizes to k
	if int(sEnd)-int(sStart) == args.k:
		printBED(quer, subj, sStart, sEnd, args.genome, strand)
	else: # expand or trim
		midpoint = int((sStart+sEnd)/2) # for odd sizes int will always round down, but just as arbitrary aa rounding up
		sStart = midpoint - 100
		sEnd = midpoint + 100
		printBED(quer, subj, sStart, sEnd, args.genome, strand)
