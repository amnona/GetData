#!/usr/bin/env python

"""
#amnonscript.py
get all the fasta files from the sra runinfo file
input:
a tab or comma delimited file (runinfo table - SraRunTable.txt - from the study runs links - i.e. http://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP056779 -> http://trace.ncbi.nlm.nih.gov/Traces/study/?acc=SRP056779)

"""

import argparse
import sys
import os.path

import csv
import subprocess

__version__ = "1.1"


def GetSRA(inputname, path, skipifthere=False, fastq=False, delimiter=None, outdir='fasta'):
	if delimiter is None:
		with open(inputname) as csvfile:
			xx = csv.Sniffer()
			res = xx.sniff(csvfile.readline(), delimiters=',\t')
			delimiter = res.delimiter
			print('Detected delimiter %s' % delimiter)

	ifile = csv.DictReader(open(inputname, 'rU'), delimiter=delimiter)
	num_files = 0
	for cline in ifile:
		if 'Run_s' in cline:
			csamp = cline['Run_s']
		elif 'Run' in cline:
			csamp = cline['Run']
		elif 'acc' in cline:
			csamp = cline['acc']
		num_files += 1

		if skipifthere:
			if os.path.isfile(os.path.join(outdir, csamp) + '.fasta'):
				print("skipping file %s. file exists" % csamp)
				continue

		print("getting file %s" % csamp)
		params = [path + 'fastq-dump', '--disable-multithreading']
		params += ['--outdir', outdir]
		if not fastq:
			params += ['--fasta']
		params += [csamp]
		print(params)
		subprocess.call(params)
		print("got file %s" % csamp)
	print('got %d files.' % num_files)
	return num_files


def main(argv):
	parser = argparse.ArgumentParser(description='Get all samples of a study from the SRA version ' + __version__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i', '--input', help='sra runinfo table for the study')
	parser.add_argument('-o', '--outdir', help='directory to place the fasta/fastq files', default='fasta')
	parser.add_argument('-p', '--path', help='path to the sratoolkit binary', default='/home/amam7564/bin/sratoolkit.2.8.0-centos_linux64/bin/')
	parser.add_argument('-s', '--skipifhere', help='if set, dont reload files already in the dir', action='store_true')
	parser.add_argument('-q', '--fastq', help='if set, output fastq instead of fasta', action='store_true')
	args = parser.parse_args(argv)
	GetSRA(args.input, args.path, args.skipifhere, fastq=args.fastq)


if __name__ == "__main__":
	main(sys.argv[1:])
