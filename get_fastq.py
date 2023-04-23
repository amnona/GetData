#!/usr/bin/env python

# amnonscript

import argparse
import sys
from collections import defaultdict
import logging

import get_sra


def main(argv):
	parser = argparse.ArgumentParser(description='Get fastq/fasta from SRA metadata table. version ' + __version__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i', '--input', help='name of input SraRunTable file (with the sample accessions', default='SraRunTable.txt')
	parser.add_argument('-p', '--sra-path', help='path to the sratoolkit binary', default='/home/amam7564/bin/sratoolkit.3.0.0-centos_linux64/bin/')
	parser.add_argument('--log-file', help='log file for the run', default='process_experiment.log')
	parser.add_argument('--log-level', help='level of log file msgs (10=debug, 20=info ... 50=critical', type=int, default=20)

	args = parser.parse_args(argv)

	logging.basicConfig(filename=args.log_file, filemode='w', format='%(levelname)s:%(message)s', level=args.log_level)
	logging.debug('get_fastq started')
	num_files = get_sra.GetSRA(infile, sra_path, skipifthere=True, outdir=fasta_dir, skip_16s_check=skip_16s_check,fastq=get_fastq)
	logging.info('downloaded %d files' % num_files)


if __name__ == "__main__":
	main(sys.argv[1:])
