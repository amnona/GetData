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

__version__ = "1.2"


def GetSRA(inputname, path, skipifthere=False, fastq=False, delimiter=None, outdir='fasta', skip_16s_check=True, split_files=False):
        '''Get all the samples from the SRA. Using the Metadata (runinfo metadata) table SraRunTable.txt from the run browser

        Parameters
        ----------
        inputname: str
                the SraRunInfo.txt file. A table containing a column with Run_s/Run/acc column that contains the SRR accession numbers.
        path: str
                path to the SraToolKit binary directory
        skipifthere: bool, optional
                if true, do not download files that already exist
        fastq: bool, optional
                if true, download fastq instead of fasta
        delimiter: str or None, optional
                delimiter for the table. If none, autodetect
        outdir: str, optional
                name of the output directory for the downloads
        skip_16s_check: bool, optional
                if True, try to identify which samples are WGS and not 16s (>500M reads, not PCR/AMPLICON) and ignore them
        split_files: bool, optional
                if True, split the samples into forward and reverse reads

        Returns
        -------
        num_files: int
                number of files downloaded
        '''
        if delimiter is None:
                with open(inputname) as csvfile:
                        xx = csv.Sniffer()
                        res = xx.sniff(csvfile.readline(), delimiters=',\t')
                        delimiter = res.delimiter
                        print('Detected delimiter %s' % delimiter)

        ifile = csv.DictReader(open(inputname, 'rU'), delimiter=delimiter)
        num_files = 0
        num_skipped = 0
        for cline in ifile:
                if 'Run_s' in cline:
                        csamp = cline['Run_s']
                elif 'Run' in cline:
                        csamp = cline['Run']
                elif 'acc' in cline:
                        csamp = cline['acc']
                num_files += 1

                # test if the sample is 16s or shotgun
                # look for some clues and also only if it is big (>500Mb)
                suspicious = False
                if 'LibrarySelection' in cline:
                        if cline['LibrarySelection'] != 'PCR':
                                suspicious = True

                if 'Assay_Type' in cline:
                        if cline['Assay_Type'] != 'AMPLICON':
                                suspicious = True

                if not skip_16s_check:
                        if suspicious:
                                if 'MBases' in cline:
                                        if int(cline['MBases']) > 500:
                                                print("skipping sample %s since it seems not 16S")
                                                num_skipped += 1
                                                continue

                if skipifthere:
                        if os.path.isfile(os.path.join(outdir, csamp) + '.fasta'):
                                print("skipping sample %s. file exists" % csamp)
                                continue

                print("getting file %s" % csamp)
                params = [os.path.join(path, 'fastq-dump'), '--disable-multithreading']
                params += ['--outdir', outdir]
                if split_files:
                        params += ['--split-files']
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
        parser.add_argument('-r', '--split-files', help='if set, split forward and reverse reads', action='store_true')
        args = parser.parse_args(argv)
        GetSRA(args.input, args.path, args.skipifhere, fastq=args.fastq, split_files=args.split_files)


if __name__ == "__main__":
        main(sys.argv[1:])
