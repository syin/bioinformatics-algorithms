#!/usr/bin/env python

"""
Needleman-Wunsch implementation in Python
"""
from sys import argv
from Bio import SeqIO
from Bio.SubsMat import MatrixInfo
import argparse

def align(seq1, seq2, d, outfile):
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
    align1 = ""
    align2 = ""
    i = len(seq2)
    j = len(seq1)

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


    write_output(outfile, align1, align2)


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
    if len(fa) == 2:
        seq1 = str(fa[0].seq)
        seq2 = str(fa[1].seq)
    else:
        print "Incorrect number of sequences in input file (should be 2)"
        exit()

    return seq1, seq2


def write_output(outfile, s1, s2):
    f = open(outfile, "w")
    f.write(s1+"\n"+s2)
    f.close()

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", nargs="+", help="Input FASTA file containing two sequences", required=True)
    parser.add_argument("-o", "--output", nargs="+", help="Output file", required=True)
    parser.add_argument("-d", "--gap-open", nargs="+", help="Gap open penalty", required=True, type=int)
    parser.add_argument("-e", "--gap-extend", nargs="+", help="Gap extend penalty. If not specified, assumed linear gap penalty.", type=int)

    args = parser.parse_args()

    try:
        filename = args.input[0]
        seq1, seq2 = read_sequences(filename)
        d = int(args.gap_open[0])
        outfile = args.output[0]
        align(seq1, seq2, d, outfile)
        
    except:
        pass


if __name__ == "__main__":
    main()
