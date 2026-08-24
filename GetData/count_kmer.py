#!/usr/bin/env python

# amnonscript
# count k-mer header freq in a directory of fasta files

import argparse
import sys
import os
from collections import defaultdict

__version__ = "0.1"


def iterfastaseqs(filename):
	"""
	iterate a fasta file and return header,sequence
	input:
	filename - the fasta file name

	output:
	seq - the sequence
	header - the header
	"""

	fl = open(filename, "r")
	cseq = ''
	chead = ''
	for cline in fl:
		if cline[0] == '>':
			if chead:
				yield(cseq, chead)
			cseq = ''
			chead = cline[1:].rstrip()
		else:
			cseq += cline.strip().replace('U', 'T')
	if cseq:
		yield(cseq, chead)
	fl.close()

def count_kmer_head(files, kmer_len=4, max_seqs=0):
	kmers = defaultdict(float)
	tot_seqs = 0
	for cfile in files:
		file_seqs = 0
		for cseq, chead in iterfastaseqs(cfile):
			chead = cseq[:kmer_len]
			kmers[chead] += 1
			tot_seqs += 1
			file_seqs += 1
			if max_seqs > 0:
				if file_seqs > max_seqs:
					break
	for k,v in kmers.items():
		kmers[k] = kmers[k] * 100 / tot_seqs
	return sorted(kmers.items(), key=lambda x: x[1], reverse=True)

def main(argv):
	parser = argparse.ArgumentParser(description='count_kmer_head version %s.' % __version__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i', '--input', help='name of input fasta file')
	parser.add_argument('-d','--dir', help='name of input directory (instead of single file')
	parser.add_argument('-l', '--len', help='length of k-mer to count', default=4)
	parser.add_argument('--max-seqs', help='number of seqs to check per file (0 for all)', default=0)

	args = parser.parse_args(argv)
	if args.dir is None:
		files = [args.input]
	else:
		files = [os.path.join(args.dir, f) for f in os.listdir(args.dir) if f.endswith('.fasta')]
	kmers = count_kmer_head(files, kmers=args.len, max_seqs=args.max_seqs)
	for k,v in kmers:
		print("%s: %s" % (k,v))


if __name__ == "__main__":
	main(sys.argv[1:])
