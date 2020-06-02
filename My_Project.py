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
 
        match = re.search(r"[^ATCGN]", self.sequence)
        self.my_seq = []
        if match:
            protein_seq = self.sequence
            self.my_seq.append("Protein")
        else:
            dna_seq = self.sequence
            self.my_seq.append("DNA")
      
    def type_seq(self):
        
        if self.my_seq[0] == "DNA":
            return dna_seq
        else:
            return "Not a DNA sequence. Run program again with DNA sequence"
        
    def GC_Content(dna_seq,self):
        c_count = dna.count('C')
        g_count = dna.count('G')
        GC_content = (c_count + g_count/self.length)*100
        return GC_content
    
    def kmer_count(self,size):
        pass

    def Complementary(self):
        pass

    def Translate(self):
        pass

Sequence = input("Enter a DNA sequence: ")

Seq = Analyze_DNA_Sequence("ATGCT")
print(Seq.GC_Content(Seq))

# Sequence = "ATGCTCGTAGAT" # Start-Leucine-Valine-Asparagine
        
