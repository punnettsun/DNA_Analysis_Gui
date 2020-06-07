# CS_22B_Project
# Author: Punit Sundar
# Last Updated: June 6, 2020
# Purpose: Write an object-oriented program for analyzing DNA sequence,
# including CG content, k-mer counting, complemnentary sequence, and
# translation of all 6 reading frames into protein sequences.

from tkinter import *
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


def get_sequence():
    sequence = entry.get()
    match = re.search(r"[^ATCG]", sequence)
    if match:
        print("Letters other than A,T,C, or G are not accepted. Try again.")
    else:
        Seq = Analyze_DNA_Sequence(sequence)
        gc_content = Seq.GC_Content()
        complement = Seq.Complementary('5-3')
        kmer = Seq.kmer_count(3)
        reading_frames = Seq.Translate()
        print("-----------GC%------------")
        print(gc_content)
        print("--------Complement--------")
        print(complement)
        print("-----------kmer-----------")
        print(kmer)
        print("----All Reading Frames----")
        for keys,values in reading_frames.items():
            print(keys,values)


root = Tk()
dna_label = Label(text = 'DNA Sequence', font=('Verdana',12,'bold'),
                  bg='blue',fg='white')
dna_label.grid(row=0, column=0)

entry = Entry()
entry.grid(row=0,column=1)

b = Button(root, text='Enter', command = get_sequence)
b.grid(row=0,column=2)

mainloop()

