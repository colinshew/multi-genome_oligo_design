#!/usr/bin/env python

#########################################################################
# Filter BLAT8 output in two steps by percent identity + alignment length
# Author: Colin Shew
# Created: 24 Mar 2020
#########################################################################

import sys
import argparse
import pandas as pd

# parse arguments
parser = argparse.ArgumentParser(description='Filter BLAT8 output in two steps by percent identity + alignment length')

parser.add_argument('blat', metavar='blat', type=str, help='BLAT output (BLAST8 format, no header)')
parser.add_argument('id1', metavar='id1', type=float, help='First filter minimum percent identity')
parser.add_argument('L1', metavar='length1', type=int, help='First filter minimum length')
parser.add_argument('n1', metavar='n1', type=int, help='First filter maximum number of hits')
parser.add_argument('id2', metavar='id2', type=float, help='Second filter minimum percent identity')
parser.add_argument('L2', metavar='length2', type=int, help='Second filter minimum length')
parser.add_argument('n2', metavar='n2', type=int, help='Second filter maximum number of hits')
parser.add_argument('logfile', metavar='logfile', type=str, help='Path to log file (contains number of hits for each query for filters 1 and 2')

args = parser.parse_args()

def printDF(df):
	"""print dataframe to stdout (.to_csv() function adds double line break)"""
	dfLOL = df.values.astype(str).tolist()
	for list in dfLOL:
		print('\t'.join(list))

def writeLog(q, nRow1, nRow2, string):
	"""log number of hits after each filter"""
	f = open(args.logfile, 'a')
	f.write('\t'.join([q, str(nRow1), str(nRow2), string]))
	f.write('\n')
	f.close()

# initialize log file
f = open(args.logfile, 'w') # overwrite if existing
f.close()

# read and format BLAT output
blat = pd.read_table(args.blat, header=None)
blat.columns = ['quer', 'subj', 'percID', 'alignL', 'mismatch', 'gapOpen', 'qStart', 'qEnd', 'sStart', 'sEnd', 'eVal', 'bitScore']
print('\t'.join(blat.columns)) # print header

# process hits one query at a time
queries = blat.quer.unique()
for q in queries:
	subset = blat[blat.quer.eq(q)]

	# filter 1
	filt = subset[(subset.percID >= args.id1) & (subset.alignL >= args.L1)]
	nRow1 = filt.shape[0]
	if 0 < nRow1 <= args.n1:
		printDF(filt)
		writeLog(q, nRow1, '.', 'pass@filt1')
		continue

	# filter 2
	filt = subset[(subset.percID >= args.id2) & (subset.alignL >= args.L2)]
	nRow2 = filt.shape[0]
	if 0 < nRow2 <= args.n2:
		printDF(filt)
		writeLog(q, nRow1, nRow2, 'pass@filt2')
		continue

	# fails both filters
	writeLog(q, nRow1, nRow2, 'fail')
