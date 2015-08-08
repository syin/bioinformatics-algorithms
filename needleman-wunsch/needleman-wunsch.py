#!/usr/bin/env python

"""
Needleman-Wunsch implementation in Python
"""
from sys import argv
from Bio import SeqIO
from Bio.SubsMat import MatrixInfo


def align(seq1, seq2, d):
    # d is gap penalty

    # Initialize dynamic programming matrix
    S = [[0 for i in range(len(seq1)+1)] for j in range(len(seq2)+1)]
    pointers = [["" for i in range(len(seq1)+1)] for j in range(len(seq2)+1)]


    # Initialization
    # set first row
    S[0] = [-d*i for i in range(len(seq1)+1)]
    # set first column
    for j,row in enumerate(S):
        row[0] = -d*float(j)

    # Calculation
    # i is row, j is column

    for i in range(len(seq2)+1)[1:]:
        for j in range(len(seq1)+1)[1:]:

            f = [S[i-1][j-1] + score((seq2[i-1],seq1[j-1])), # match
                 S[i-1][j] - d,                              # gap
                 S[i][j-1] - d]                              # gap

            S[i][j] = max(f)


            argmax = f.index(max(f))
            pointers[i][j] = argmax


    # Termination
    final_score = S[len(seq2)][len(seq1)]

    # Traceback
    align1 = "", align2 = ""
    i = len(seq2), j = len(seq1)

    while True:
        if i==0 or j==0:
            break

        # argmax == 0: diagonal
        # argmax == 1: up
        # argmax == 2: left

        if pointers[i][j] == 0:
            align1 = seq1[j-1] + align1
            align2 = seq2[i-1] + align2

            i = i-1
            j = j-1
        elif pointers[i][j] == 1:
            align1 = "-" + align1
            align2 = seq2[i-1] + align2

            i = i-1

        elif pointers[i][j] == 2:
            align1 = seq1[j-1] + align1
            align2 = "-" + align2

            j = j-1

    print align1
    print align2



def score(pair):

    blosum = MatrixInfo.blosum62

    if pair in blosum:
        return blosum[pair]
    else:
        return blosum[tuple(reversed(pair))]


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
    d = float(argv[2])
    seq1, seq2 = read_sequences(filename)
    align(seq1, seq2, d)



if __name__ == "__main__":
    main()
