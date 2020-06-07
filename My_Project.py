# CS_22B_Project
# Author: Punit Sundar
# Last Updated: June 6, 2020
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


## Warnings:
# 1. Make sure DNA sequence does not have any odd characters
# 2. Make sure k-mer size does not exceed length of the DNA sequence
# 3. If user inputs are not the correct type or value causing exceptions


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
        dna_to_protein = {
            'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
            'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
            'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
            'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
            'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
            'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
            'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
            'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
            'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
            'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
            'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
            'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
            'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
            'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
            'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
            'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
            }
        length = self.length
        reading = {}
        for i in range(3):
            reading['frame_'+str(i+1)] = tuple([dna_to_protein[self.sequence[index:index+3]] for index in range(i,length-2,3)])
        reverse_strand = Analyze_DNA_Sequence.Complementary(self,'5-3')
        for i in range(3):
            reading['frame_'+str(i+4)] = tuple([dna_to_protein[reverse_strand[index:index+3]] for index in range(i,length-2,3)])

        return reading
        #ATGTTTGGATAG = 1. MFG*
        #               2. CLD
        #               3. VWI
        #CTATCCAAACAT = 4. LSKH
        #               5. YPN
        #               6. IQT



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
    print(Seq.Translate())









