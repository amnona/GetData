#!/usr/bin/env python

"""
amnonscript
Trim primer region from sequences
original script by Jenya
modified by Amnon for command line interface and trimming the primer sequences
"""

__version__ = "1.6"

import argparse

# from cogent.parse.fasta import MinimalFastaParser
# import skbio

import sys
import re
import os


def get_region_single(inputname, fprimer, rprimer, length, remove_ambig, keep_primers, skip_reverse, output_mismatch=False, output_file=None):
    '''
    output_file: str or None, optional
        if not none, save to output_file. otherwise, print to stdout
    '''
    if output_file is None:
        outf = None
    else:
        outf = open(output_file, 'w')

    fplen = 0
    for cseq, seqid in iterfastaseqs(inputname):
        cseq = cseq.upper()

        reverse_primer_seq = ''
        # find the start of the primer and output all following sequence
        try:
            match = re.search(fprimer, cseq)
            if keep_primers:
                forward_primer_seq = cseq[match.start():]
                fplen = match.end() - match.start()
            else:
                forward_primer_seq = cseq[match.end():]
        except AttributeError:
            forward_primer_seq = ''
            if output_mismatch:
                    print(">%s\n%s" % (seqid, cseq))

        # if sequence can be amplified, search for reverse primer
        if forward_primer_seq != '':
            if not skip_reverse:
                # find the reverse primer match and use only up to it
                try:
                    match = re.search(rprimer, forward_primer_seq)
                    if keep_primers:
                        reverse_primer_seq = forward_primer_seq[:match.end()]
                    else:
                        reverse_primer_seq = forward_primer_seq[:match.start()]
                except AttributeError:
                    reverse_primer_seq = ''
            else:
                # reverse_primer_seq = forward_primer_seq[:length + fplen]
                reverse_primer_seq = forward_primer_seq[:]

            if reverse_primer_seq != '':
                # trim length if needed
                if length > 0:
                    if keep_primers:
                        reverse_primer_seq = reverse_primer_seq[:length + fplen]
                    else:
                        reverse_primer_seq = reverse_primer_seq[:length]
                printit = True
                if output_mismatch:
                    printit = not printit
                if remove_ambig:
                    if 'N' in reverse_primer_seq:
                        printit = False
                if printit:
                    outstr = ">%s\n%s" % (seqid, reverse_primer_seq)
                    if outf is None:
                        print(outstr)
                    else:
                        outf.write(outstr + '\n')


def get_region(inputname, outputname, fprimer, rprimer=None, length=0, remove_ambig=False, keep_primers=False, skip_reverse=True, output_mismatch=False, output_file=None):
    if os.path.isfile(inputname):
        filelist = [inputname]
        outfilelist = [outputname]
    else:
        if not os.path.exists(outputname):
            os.makedirs(outputname)
        filelist = [os.path.join(inputname, f) for f in os.listdir(inputname) if os.path.isfile(os.path.join(inputname, f))]
        outfilelist = [os.path.join(outputname, f) for f in os.listdir(inputname) if os.path.isfile(os.path.join(inputname, f))]

    for cfile, coutfile in zip(filelist, outfilelist):
        # add only .fasta/.fa files from the dir
        if cfile[-6:] == '.fasta' or cfile[-3:] == '.fa':
            get_region_single(cfile, fprimer, rprimer, length, remove_ambig, keep_primers, skip_reverse, output_mismatch, output_file=coutfile)


def iterfastaseqs(filename):
    """
    iterate a fasta file and return header,sequence
    input:
    filename - the fasta file name

    output:
    seq - the sequence
    header - the header
    """

    fl = open(filename, "rU")
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


def main(argv):
    parser = argparse.ArgumentParser(description='Extract an amplified region from a set of primers version ' + __version__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', help='name of input fasta file or directory')
    parser.add_argument('-o', '--output', help='name of output fasta file or directory (created if not exist)', default='region')
    parser.add_argument('-f', '--fprimer', help='forward primer sequence (default is V4f)', default='GTGCCAGC[AC]GCCGCGGTAA')
    parser.add_argument('-r', '--rprimer', help='reverse primer sequence (in reverse complement) (default is V4r)', default='ATTAGA[AT]ACCC[CGT][AGT]GTAGTCC')
    parser.add_argument('-l', '--length', help='Trim all sequences to length (0 for full length)', default='0')
    parser.add_argument('-n', '--remove_ambig', help='Remove seqs containing ambiguous characters (N)', action='store_true')
    parser.add_argument('-s', '--skip_reverse', help='Dont look for reverse primer', action='store_true')
    parser.add_argument('-k', '--keep_primers', help="Don't remove the primer sequences", action='store_true')
    parser.add_argument('-m', '--output_mismatch', help="show sequences with no primer instead", action='store_true')

    args = parser.parse_args(argv)

    get_region(args.input, args.output, args.fprimer, args.rprimer, int(args.length), args.remove_ambig, args.keep_primers, args.skip_reverse, args.output_mismatch)

    # if os.path.isfile(args.input):
    #     get_region(args.input, args.fprimer, args.rprimer, int(args.length), args.remove_ambig, args.keep_primers, args.skip_reverse, args.output_mismatch, output_file=args.output)
    # else:
    #     if not os.path.exists(args.output):
    #         os.makedirs(args.output)
    #     filelist = [f for f in os.listdir(args.input) if os.path.isfile(os.path.join(args.input, f))]
    #     for cfile in filelist:
    #         # add only .fasta/.fa files from the dir
    #         if cfile[-6:] == '.fasta' or cfile[-3:] == '.fa':
    #             get_region(os.path.join(args.input, cfile), args.fprimer, args.rprimer, int(args.length), args.remove_ambig, args.keep_primers, args.skip_reverse, args.output_mismatch, output_file=os.path.join(args.output, cfile))


if __name__ == "__main__":
    main(sys.argv[1:])
