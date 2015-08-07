#!/usr/bin/env python

"""
Needleman-Wunsch implementation in Python
"""

from sys import argv
from Bio import SeqIO


def align(seq1, seq2):

	# Initialization

	# Calculation

	# Termination

	# Traceback


def read_sequences(filename):
	"""
	Read FASTA file and obtain the two sequences
	"""

	fa = list(SeqIO.parse(filename, "fasta"))
	seq1 = str(fa[0].seq)
	seq2 = str(fa[1].seq)

	return seq1, seq2


def main():

	filename = argv[1]
	seq1, seq2 = read_sequences(filename)
	align(seq1, seq2)



if __name__ == "__main__":
    main()
