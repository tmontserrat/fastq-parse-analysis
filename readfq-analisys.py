#!/usr/bin/env python
# nuccount.py -- tally nucleotides in a file
import sys
from collections import Counter
from readfq import readfq
# readfq Ref: https://github.com/lh3/readfq/blob/master/readfq.py

IUPAC_BASES = "ACGTRYSWKMBDHVN-."
PATH = "contam.fastq"

# intialize counter
counts = {}

# counting variable for count the number of sequences
s = 0

# list for soft-masked characters (lower case letters)
lower = []

# lists to append the sequences length and the gc content per sequence
seqs_length = []
seqs_gc_content = []

with open("contam.fastq") as file_object:
    contents = file_object.read()
    contents = str(contents.strip())

with open("contam.fastq") as file:

    for name, seq, qual in readfq(file):
        number_lower = 0
        number_gc = 0
        # each time there's a new name, add 1
        s += 1
        # each time there's a new sequence, append its length
        seqs_length.append(len(seq))
        # in eac sequence, compute the number of G's and C's
        for nucleotide in seq:
            # make nucleotide case insensitive
            if nucleotide.upper() == "G" or nucleotide.upper() == "C":
                number_gc += 1
            # Create the dictionary on the fly
            if nucleotide.upper() in counts:
                counts[str(nucleotide.upper())] += 1
            else:
                counts[str(nucleotide.upper())] = 1
            if nucleotide.islower():
                number_lower += 1
        lower.append(number_lower)
            
        # once compute the number of G's and C's, compute its proportion
        seqs_gc_content.append(number_gc/len(seq))

# print the results
for base in IUPAC_BASES:
    if base in counts:
        print(base + "\t" + str(counts[base]))
    else:
        counts[base] = 0
        print(base + "\t" + str(counts[base]))

print(f"\n---Number of sequences---\n{s}")

number_seq = 0

print("\n---Sequence lengths---")
for sequence in seqs_length:
    number_seq += 1
    print(f"Sequence number {number_seq}: {sequence}")

number_seq = 0

print("\n---GC content---")
for gc_content in seqs_gc_content:
    number_seq += 1
    print(f"Sequence number {number_seq}: {round(gc_content, 4)}")

number_seq = 0

print ("\n---Soft-masked characters---")
for sequence_number in range(s):
    print(f"Sequence number {sequence_number+1}: {lower[sequence_number]}")
    sequence_number += 1