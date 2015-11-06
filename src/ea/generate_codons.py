# -*- coding: utf-8 -*-
"""
Read and load a black and white image into a table associating
pairs of x and y coordinates (starting from top-leftmost pixel) to the
corresponding color/label (white=0, black=1)
"""

import sys
import numpy as np


if __name__ == "__main__":

    size_codons = int(sys.argv[1])
    number_codon_pairs = int(sys.argv[2])
    list_codons = []

    while len(list_codons) < number_codon_pairs:
        number_ones = 0
        number_zeros = 0
        start = ""
        stop = ""
        while (len(start) < size_codons):
            if (number_zeros < (size_codons/2)) and (number_ones <
                                                     (size_codons/2)):
                character = np.random.choice(['0', '1'])
            else:
                if number_zeros == (size_codons/2):
                    character = '1'
                else:
                    if number_ones == (size_codons/2):
                        character = '0'

            if character == '0':
                number_zeros = number_zeros + 1
            else:
                if character == '1':
                    number_ones = number_ones + 1
            start = start + character

        number_ones = 0
        number_zeros = 0

        while len(stop) < size_codons:
            if (number_zeros < (size_codons/2)) and (number_ones <
                                                     (size_codons/2)):
                character = np.random.choice(['0', '1'])
            else:
                if number_zeros == (size_codons/2):
                    character = '1'
                else:
                    if number_ones == (size_codons/2):
                        character = '0'
            if character == '0':
                number_zeros = number_zeros + 1
            else:
                if character == '1':
                    number_ones = number_ones + 1

            stop = stop + character
        if start != stop:
            list_codons.append((start, stop))

    for l in list_codons:
        print l[0], " ", l[1]
