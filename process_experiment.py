#!/usr/bin/env python

# amnonscript

import argparse
import sys
import os
import re
from collections import defaultdict
import subprocess

import numpy as np

import get_sra
import get_region

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


def test_fasta_file(files, primers={'AGAGTTTGATC[AC]TGG[CT]TCAG': 'v1', 'CCTACGGG[ACGT][CGT]GC[AT][CG]CAG': 'v3', 'GTGCCAGC[AC]GCCGCGGTAA': 'v4'}, max_start=10, min_primer_len=10, num_reads=1000, min_fraction=0.25, min_files_fraction=0.5):
	'''Check if the fasta file starts with one of a given set of primers.

	Parameters
	----------
	filename: str
		the fasta file name to test
	primers: dict of {primer(str): region_name(str)}, optional
		the primers to test for
	max_start: int, optional
		maximal start position for the primer (i.e. do not return if primer starts after position max_start)
	min_primer_len: int, optional
		trim primers to keep only min_primer_len last chars
	num_reads: int, optional
		the number of reads in the fasta file to process
	min_fraction: float, optional
		need at least min_fraction sequence matches in order to return the primer

	Returns
	-------
	primer: str or None
		the primer identified as maximal (if >min_fraction matches in tested reads)
	primer_name: str or None
		the name of the primer region identified
	'''
	# trim the primers if needed
	if min_primer_len is not None:
		new_primers = {}
		for k, v in primers.items():
			pos = len(k)
			numchars = 0
			newp = ''
			while True:
				if numchars >= min_primer_len:
					break
				pos = pos - 1
				if pos < 0:
					break
				if k[pos] != ']':
					newp = k[pos] + newp
					numchars += 1
					continue
				while k[pos] != '[':
					newp = k[pos] + newp
					pos = pos - 1
				newp = k[pos] + newp
				numchars += 1
			new_primers[newp] = v
		primers = new_primers

	# scan the files
	all_matches = defaultdict(float)
	for cfile in files:
		matches = defaultdict(int)
		num_tested = 0
		for cseq, chead in iterfastaseqs(cfile):
			for cprimer in primers.keys():
				ccseq = cseq[:max_start + len(cprimer)]
				match = re.search(cprimer, ccseq)
				if match is not None:
					matches[cprimer] += 1
			num_tested += 1
			if num_tested > num_reads:
				break
		if len(matches) > 0:
			max_primer = max(matches, key=matches.get)
			match_fraction = matches[max_primer] / num_tested
			if match_fraction > min_fraction:
				all_matches[max_primer] += 1
	if len(all_matches) > 0:
		maxregion = max(all_matches, key=all_matches.get)
		if all_matches[maxregion] / len(files) > min_files_fraction:
			return maxregion, primers[maxregion]
	return None, None


def test_read_length(files, num_reads=1000, prctile=75):
	'''get the typical read length for files
	Parameters
	----------
	files: list of str
		the fasta files to test
	num_reads: int, optional
		the number of reads in the fasta file to process
	prctile: float, optional
		the percentile to use as the typical read length

	Returns
	-------
	read_length: int
		the read length
	'''
	all_reads = []
	for cfile in files:
		num_tested = 0
		for cseq, chead in iterfastaseqs(cfile):
			all_reads.append(len(cseq))
			num_tested += 1
			if num_tested > num_reads:
				break
	return int(np.percentile(all_reads, prctile))


def test_kmer_head_region(files, kmers={'v4': ['TACG'], 'v3': ['TGGG', 'TGAG'], 'v1': ['GACG', 'GATG', 'ATTG']}, num_reads=1000, min_fraction=0.5, min_files_fraction=0.5):
	'''Test if a fasta file starts with known region k-mers

	Parameters
	----------
	infile: str
		name of fasta file
	kmers: dict of {region_name(str): [kmers(str)]}
		the kmers expected to begin each primer region
		get the values using ~/scripts/count_kmer.py -d v1
	num_reads: int
		the maximal number of reads per fileto test
	min_fraction: float
		the minimal expected fraction of reads in the region starting with the kmers (summed over all kmers of the region)
	min_files_fraction: float
		the minimal fraction of files positive for the primer in order to identify the experiment

	Returns
	-------
	region: str or None
		the matching region or None if no region
	'''
	file_primers = defaultdict(float)
	kmer_len = len(kmers['v4'][0])
	for cfile in files:
		print(cfile)
		num_tested = 0
		kmer_dist = defaultdict(float)
		for cseq, chead in iterfastaseqs(cfile):
			cseq = cseq[:kmer_len]
			for cregion, ckmers in kmers.items():
				if cseq in ckmers:
					kmer_dist[cregion] += 1
			num_tested += 1
			if num_tested > num_reads:
				break

		if len(kmer_dist) > 0:
			for k, v in kmers.items():
				kmer_dist[k] = kmer_dist[k] / num_tested
			maxregion = max(kmer_dist, key=kmer_dist.get)
			if kmer_dist[maxregion] > min_fraction:
				file_primers[maxregion] += 1
	if len(file_primers) == 0:
		print('no exact regions identified in files')
		return None
	maxregion = max(file_primers, key=file_primers.get)
	if file_primers[maxregion] / len(files) > min_files_fraction:
		return maxregion
	return None


def process_experiment(infile, sra_path, fasta_dir='fasta', max_test=10, skip_get=False, seq_len=150):
	'''download the Sra table, convert to known region, and deblur

	Parameters
	----------
	infile: str
		name of the input SraRunTable (tab or comma delimited)
	fasta_dir: str, optional
		name of the output fasta directory for the sra download
	max_test: int, optional
		maximal number of files to check for primer/region
	'''
	# get all the fasta files
	if not skip_get:
		print('processing sratable %s' % infile)
		num_files = get_sra.GetSRA(infile, sra_path, skipifthere=True, outdir=fasta_dir)
		print('downloaded %d files' % num_files)
	else:
		print('skipping getting files from sra')

	# check if known region / if we need to trim primer
	files = [os.path.join(fasta_dir, f) for f in os.listdir(fasta_dir) if f.endswith('.fasta')]
	if len(files) == 0:
		raise ValueError('no fasta files found in %d' % fasta_dir)
	if len(files) > max_test:
		files = [files[x] for x in np.random.permutation(len(files))[:max_test]]
	region = test_kmer_head_region(files)
	if region is not None:
		print('region is %s and no primer trimming needed' % region)
	else:
		match_primer, match_primer_name = test_fasta_file(files)
		if match_primer is not None:
			print('trimming with primer %s for region %s' % (match_primer, match_primer_name))
			trim_dir = 'trim'
			get_region.get_region(fasta_dir, outputname=trim_dir, fprimer=match_primer, skip_reverse=True)
			fasta_dir = trim_dir
			print('finished trimming')
		else:
			raise ValueError('No matching regions/primers. please check manually!')

	# check the length of typical reads
	read_len = test_read_length(files)
	print('typical read length = %d' % read_len)
	if read_len < 100:
		raise ValueError('Read length %d too short' % read_len)
	read_len = np.min([read_len, seq_len])
	print('deblurring')
	# deblur workflow --seqs-fp fasta --output-dir deblur -w -t 150 -O 32 --min-reads 10 --pos-ref-db-fp /home/amam7564/data/icu/deblur/deblur_working_dir/88_otus --neg-ref-db-fp /home/amam7564/data/icu/deblur/deblur_working_dir/artifacts
	params = []
	# params += ['qsub', '-d', '$PWD', '-V', '-m', 'abe', '-M', 'amnonimjobs@gmail.com', '-j', 'eo', '-e', 'process.err', '-l', 'walltime=48:00:00,nodes=1:ppn=32,mem=250gb', '-N', 'process']
	params += ['deblur', 'workflow']
	params += ['--seqs-fp', fasta_dir]
	params += ['--output-dir', 'deblur']
	params += ['-w', '-t', str(read_len)]
	params += ['-O', '32', '--min-reads', '10']
	params += ['--pos-ref-db-fp', '/home/amam7564/data/icu/deblur/deblur_working_dir/88_otus']
	params += ['--neg-ref-db-fp', '/home/amam7564/data/icu/deblur/deblur_working_dir/artifacts']
	subprocess.call(params)
	print('done')


def main(argv):
	parser = argparse.ArgumentParser(description='Process experiment version ' + __version__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i', '--input', help='name of input SraRunTable file (with the sample accessions')
	parser.add_argument('-p', '--sra-path', help='path to the sratoolkit binary', default='/home/amam7564/bin/sratoolkit.2.8.0-centos_linux64/bin/')
	parser.add_argument('--skip-get', help='if set, skip getting the fasta files from SRA', action='store_true')

	args = parser.parse_args(argv)
	process_experiment(infile=args.input, sra_path=args.sra_path, skip_get=args.skip_get)


if __name__ == "__main__":
	main(sys.argv[1:])
