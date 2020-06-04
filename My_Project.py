# CS_22B_Project
# Author: Punit Sundar
# Last Updated: May 25, 2020
# Purpose: Write an object-oriented program for analyzing DNA sequence,
# including CG content, k-mer counting, complemnentary sequence, and
# translation of all 6 reading frames into protein sequences.

# Given a DNA sequence, do the following:
# 1. CG content
    # count C's and G's to get CG content
# 2. k-mer count
    # User gives a k-mer count to calculate
# 3. Complementary sequence
    # Given the DNA sequence, give the complementary sequence
    # Ask user if 5' to 3' or 3' to 5'
# 4. Translate in all 6 reading frames
    # AGTCGTAGATGATCGTAGTAGATCGTAGTAG
    # 3 in forward strand and 3 in reverse strand
    # Make program to ask about start or end codon


## Warnings:
# 1. Make sure DNA sequence does not have any odd characters
# 2. Length is a multiple of 3 to ensure correct codon count
# 3. Make sure k-mer size does not exceed length of the DNA sequence
# 4. If user inputs are not the correct type or value causing exceptions
# 5. 

import re

class Analyze_DNA_Sequence:
    
    def __init__(self, sequence): # initializes all the attributes of the object
        self.sequence = sequence
        self.length = len(sequence)
        
    def GC_Content(self):
        """Calculates the GC% of a given DNA sequence"""
        GC_content = lambda dna: (dna.count('G')+dna.count('C'))\
                         /self.length
        return round(GC_content(self.sequence),4)

    def kmer_count(self,size):
        """Calculates all possible kmers given a size"""
        if size == 1:
            return ['A','T','C','G']
        else:
            result = []
            for seq in Analyze_DNA_Sequence.kmer_count(self,size-1):
                for base in ['A','T','C','G']:
                    result.append(seq+base)
            return result
    
    def Complementary(self,direction):
        """Gives the complementary sequence to the given DNA sequence"""
        complementary = self.sequence.upper().replace('A','t')\
                        .replace('T','a').replace('G','c')\
                        .replace('C','g').upper()
        if direction == '3-5':
            return complementary
        else:
            return complementary[::-1]

    def Translate(self):
        """Gives all 6 reading frame protein sequences for given DNA sequence"""
        pass

dna_seq = input("Enter a DNA sequence: ")
match = re.search(r"[^ATCGN]", dna_seq)
if match:
    print("Not a DNA sequence. Try again with a DNA sequence")
else:
    dna_seq = dna_seq
    Seq = Analyze_DNA_Sequence(dna_seq)
    print(Seq.GC_Content())
    print(Seq.Complementary('5-3'))
    print(Seq.kmer_count(4))









