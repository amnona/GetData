#!/usr/bin/env python

# amnonscript

import argparse
import sys
import os
import re
from collections import defaultdict
import subprocess
import logging

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


def rev_comp_fasta(infile, outdir, reverse=True, complement=True):
	if not os.path.exists(outdir):
		os.makedirs(outdir)
	outfile = os.path.join(outdir, os.path.basename(infile))
	logging.debug('revcomp file %s into %s' % (infile, outfile))
	with open(outfile, 'w') as ofl:
		for cseq, chead in iterfastaseqs(infile):
			if complement:
				cseq = cseq.lower()
				cseq = cseq.replace('a', 'T')
				cseq = cseq.replace('c', 'G')
				cseq = cseq.replace('g', 'C')
				cseq = cseq.replace('t', 'A')
				cseq = cseq.upper()
			if reverse:
				cseq = cseq[::-1]
			ofl.write('>' + chead + '\n')
			ofl.write(cseq + '\n')


def test_fasta_file(files, base_dir=None, primers={'AGAGTTTGATC[AC]TGG[CT]TCAG': 'v1', 'CCTACGGG[ACGT][CGT]GC[AT][CG]CAG': 'v3', 'GTGCCAGC[AC]GCCGCGGTAA': 'v4'}, max_start=25, min_primer_len=10, num_reads=1000, min_fraction=0.25, min_files_fraction=0.2):
	'''Check if the fasta file starts with one of a given set of primers.

	Parameters
	----------
	filename: str
		the fasta file name to test
	base_dir: str or None, optional
		the directory where the files reside, or None to assume exact (full path) file names in files
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
	# attach the base_dir if needed
	logging.debug('Testing %d files for %d primers' % (len(files), len(primers)))

	if base_dir is not None:
		files = [os.path.join(base_dir, x) for x in files]

	# trim the primers if needed
	if min_primer_len is not None:
		logging.debug('Trimming primers before test to length %d' % min_primer_len)
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
		logging.debug('Trimmed primers are: %s' % primers)

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
		logging.debug('best matching region is %s' % maxregion)
		if all_matches[maxregion] / len(files) >= min_files_fraction:
			logging.debug('enough matches found for primer %s: %s' % (maxregion, all_matches[maxregion]))
			return maxregion, primers[maxregion]
		else:
			logging.debug('not enough matches per primer. details: %s' % primers)
	else:
		logging.debug('no matches found for any primer')
	logging.info('No match for any of %d primers found' % len(primers))
	return None, None


def test_read_length(files, base_dir=None, num_reads=1000, prctile=75):
	'''get the typical read length for files
	Parameters
	----------
	base_dir: str or None, optional
		the directory where the files reside, or None to assume exact (full path) file names in files
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
	if base_dir is not None:
		files = [os.path.join(base_dir, x) for x in files]
	all_reads = []
	for cfile in files:
		num_tested = 0
		for cseq, chead in iterfastaseqs(cfile):
			all_reads.append(len(cseq))
			num_tested += 1
			if num_tested > num_reads:
				break
	return int(np.percentile(all_reads, prctile))


def test_kmer_head_region(files, base_dir=None, kmers={'v4': ['TACG'], 'v3': ['TGGG', 'TGAG'], 'v1': ['GACG', 'GATG', 'ATTG']}, num_reads=1000, min_fraction=0.5, min_files_fraction=0.5):
	'''Test if a fasta file starts with known region k-mers

	Parameters
	----------
	infile: str
		name of fasta file
	base_dir: str or None, optional
		the directory where the files reside, or None to assume exact (full path) file names in files
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
	if base_dir is not None:
		files = [os.path.join(base_dir, x) for x in files]
	file_primers = defaultdict(float)
	kmer_len = len(kmers['v4'][0])
	logging.debug('testing kmer head for region on %d files:' % len(files))
	for cfile in files:
		logging.debug(cfile)
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
		logging.info('no exact regions identified in files')
		return None
	maxregion = max(file_primers, key=file_primers.get)
	if file_primers[maxregion] / len(files) > min_files_fraction:
		return maxregion
	return None


def process_experiment(infile, sra_path, fasta_dir='fasta', max_test=10, skip_get=False, seq_len=150, skip_16s_check=False, skip_region=False, deblur_path=None, num_threads=1, max_primer_start=25, skip_exact=False):
	'''download the Sra table, convert to known region, and deblur

	Parameters
	----------
	infile: str
		name of the input SraRunTable (tab or comma delimited)
	fasta_dir: str, optional
		name of the output fasta directory for the sra download
	max_test: int, optional
		maximal number of files to check for primer/region
	skip_get: bool, optional
		True to skip the SRA downloading step (assumes fasta files are in the fasta/ dir)
	seq_len: int, optional
		the length to trim each sequence after primer removal (actual length is min(seq_len, actual length))
	skip_16s_check: bool, optional
		True to skip the validation step that the sample is not WGS before downloading
	skip_region:bool, optional
		True to skip the region validation step (always process initial fasta reads)
	deblur_path: str or None, optional
		if not None, path to the directory containined the preprocessed artifacts/repseqs files for deblur
		(88_otus.bursttrie_0.dat, 88_otus.kmer_0.dat, 88_otus.pos_0.dat, 88_otus.stats, artifacts.bursttrie_0.dat, artifacts.kmer_0.dat, artifacts.pos_0.dat, artifacts.stats)
		if None, deblur will preprocess
	num_threads: int, optional
		the number of threads to run deblur with
	max_primer_start: int, optional
		the maximal allowed offset for the primer within the reads (so primer does not appear in the middle of the sequence - i.e. v4 in v34)
	skip_exact: bool, optional
		if True, skip the exact region match test (no trimming) - assume it is not exact
	'''
	# get all the fasta files
	if not skip_get:
		logging.info('processing sratable %s' % infile)
		num_files = get_sra.GetSRA(infile, sra_path, skipifthere=True, outdir=fasta_dir, skip_16s_check=skip_16s_check)
		logging.info('downloaded %d files' % num_files)
	else:
		logging.info('skipping getting files from sra')

	# check if known region / if we need to trim primer
	files = [f for f in os.listdir(fasta_dir) if f.endswith('.fasta')]
	found_it = False
	if not skip_region:
		if len(files) == 0:
			raise ValueError('no fasta files found in %d' % fasta_dir)
		if len(files) > max_test:
			test_files = [files[x] for x in np.random.permutation(len(files))[:max_test]]
		else:
			test_files = files

		if not skip_exact:
			# test if the sequences are of some known region
			region = test_kmer_head_region(test_files, fasta_dir)
		else:
			region = None
		if region is not None:
			logging.info('region is %s and no primer trimming needed' % region)
			found_it = True
		else:
			# test if sequences contain known primer
			match_primer, match_primer_name = test_fasta_file(test_files, fasta_dir, max_start=max_primer_start)

			# no match for primer - let's try reverse-complement
			if match_primer is None:
				rc_dir = 'revcomp'
				for cfile in files:
					rev_comp_fasta(os.path.join(fasta_dir, cfile), rc_dir)
				fasta_dir = rc_dir
				region = test_kmer_head_region(test_files, fasta_dir)
				if region is not None:
					logging.info('Found exact region %s after reverse complement')
					found_it = True
				else:
					# test if sequences contain known primer
					match_primer, match_primer_name = test_fasta_file(test_files, fasta_dir, max_start=max_primer_start)

			# if found matching primer in sequences, trim it
			if match_primer is not None:
				logging.info('trimming with primer %s for region %s' % (match_primer, match_primer_name))
				trim_dir = 'trim'
				get_region.get_region(fasta_dir, outputname=trim_dir, fprimer=match_primer, skip_reverse=True)
				fasta_dir = trim_dir
				logging.info('finished trimming')
				found_it = True
			# after all these tries didn't identify reads as coming from any known region
			if not found_it:
				raise ValueError('No matching regions/primers. please check manually!')

	# check the length of typical reads
	read_len = test_read_length(files, fasta_dir)
	logging.info('typical read length = %d' % read_len)
	if read_len < 100:
		raise ValueError('Read length %d too short' % read_len)
	read_len = np.min([read_len, seq_len])
	logging.info('deblurring')
	# deblur workflow --seqs-fp fasta --output-dir deblur -w -t 150 -O 32 --min-reads 10 --pos-ref-db-fp /home/amam7564/data/icu/deblur/deblur_working_dir/88_otus --neg-ref-db-fp /home/amam7564/data/icu/deblur/deblur_working_dir/artifacts
	params = []
	# params += ['qsub', '-d', '$PWD', '-V', '-m', 'abe', '-M', 'amnonimjobs@gmail.com', '-j', 'eo', '-e', 'process.err', '-l', 'walltime=48:00:00,nodes=1:ppn=32,mem=250gb', '-N', 'process']
	params += ['deblur', 'workflow']
	params += ['--seqs-fp', fasta_dir]
	params += ['--output-dir', 'deblur']
	params += ['-w', '-t', str(read_len)]
	params += ['-O', str(num_threads), '--min-reads', '10']
	if deblur_path is not None:
		params += ['--pos-ref-db-fp', os.path.join(deblur_path, '88_otus')]
		params += ['--pos-ref-db-fp', os.path.join(deblur_path, 'artifacts')]
	subprocess.call(params)
	logging.info('done')


def main(argv):
	parser = argparse.ArgumentParser(description='Process experiment version ' + __version__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i', '--input', help='name of input SraRunTable file (with the sample accessions')
	parser.add_argument('-p', '--sra-path', help='path to the sratoolkit binary', default='/home/amam7564/bin/sratoolkit.2.8.0-centos_linux64/bin/')
	parser.add_argument('-t', '--trim-length', help='length to trim seqs after primer removal', default=150, type=int)
	parser.add_argument('--skip-16s-check', help='download also samples that seem to be non-16s', action='store_true')
	parser.add_argument('--skip-get', help='if set, skip getting the fasta files from SRA', action='store_true')
	parser.add_argument('--skip-region', help='if set, skip validating/trimming region for primers (just process fasta)', action='store_true')
	parser.add_argument('--skip-exact', help='if set, assume the sequence is not an exact region match (force primer checking)', action='store_true')
	parser.add_argument('--max-primer-start', help='the maximal offset (from read start) for primer end', default=25, type=int)
	parser.add_argument('--log-file', help='log file for the run', default='process_experiment.log')
	parser.add_argument('--log-level', help='level of log file msgs (10=debug, 20=info ... 50=critical', type=int, default=20)
	parser.add_argument('--deblur-path', help='location of deblur pre-compiled artifacts/rep seqs')
	parser.add_argument('--num-threads', help='number of threads to run for deblur', default=1)

	args = parser.parse_args(argv)

	logging.basicConfig(filename=args.log_file, filemode='w', format='%(levelname)s:%(message)s', level=args.log_level)
	logging.debug('process_experiment started')
	process_experiment(infile=args.input, sra_path=args.sra_path, skip_get=args.skip_get, seq_len=args.trim_length, skip_16s_check=args.skip_16s_check, skip_region=args.skip_region, deblur_path=args.deblur_path, num_threads=args.num_threads, max_primer_start=args.max_primer_start, skip_exact=args.skip_exact)


if __name__ == "__main__":
	main(sys.argv[1:])
