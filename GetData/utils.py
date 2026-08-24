def iterfastaseqs(filename):
	"""
	iterate a fasta or fastq file and return header,sequence
	input:
	filename - the fasta/fastq file name

	output:
	seq - the sequence
	header - the header
	"""

	fl = open(filename, "r")
	cseq = ''
	chead = ''
	skip_next = False

	for cline in fl:
		if skip_next:
			skip_next = False
			continue
		if cline[0] == '>' or cline[0] == '@':
			if chead:
				yield(cseq, chead)
			cseq = ''
			chead = cline[1:].rstrip()
		elif cline[0] == '+':
			# need to skip the next line since it is the quality scores
			skip_next = True
		else:
			cseq += cline.strip().replace('U', 'T')
	if cseq:
		yield(cseq, chead)
	fl.close()
