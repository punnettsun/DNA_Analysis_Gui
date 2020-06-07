# CS_22B_Project
# Author: Punit Sundar
# Last Updated: June 6, 2020
# Purpose: Write an object-oriented program for analyzing DNA sequence,
# including CG content, k-mer counting, complemnentary sequence, and
# translation of all 6 reading frames into protein sequences.

from tkinter import *
from tkinter.scrolledtext import ScrolledText
import re

#tkinter add
#1. kmer size option

#Warnings
#1. Incorrect kmer size value

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
    sequence = dna_entry.get()
    kmer_size = kmer_size_entry.get()
    
    match = re.search(r"[^ATCG]", sequence)
    if match:
        print("Letters other than A,T,C, or G are not accepted. Try again.")
    else:
        Seq = Analyze_DNA_Sequence(sequence)
        gc_content = Seq.GC_Content()
        complement = Seq.Complementary('5-3')
        kmer = Seq.kmer_count(int(kmer_size))
        reading_frames = Seq.Translate()
        gc_content_text.insert(1.0,gc_content)
        complement_text.insert(END,complement)
        reading_frames_text.insert(END,reading_frames)
        textbox.insert(END,kmer)

root = Tk()
root.title("DNA GUI")

Directions = Label(text = "Directions: Paste a DNA sequence and kmer size below. Click 'Enter'",
                   font=('Verdana',12,'bold'),
                   bg='black',fg='white',width = 55, height=5)
Directions.grid(row=0,column=0,columnspan=2)

my_list = ['DNA Sequence','GC%','Type Kmer Size','Complement','Reading Frames']
for i in range(1,6):
    label = Label(text = my_list[i-1], font=('Verdana',12,'bold'),
                  fg='black',width = 25, height=5)
    label.grid(row= i, column=0)

dna_entry = Entry()
dna_entry.grid(row=1,column=1)

gc_content_text = Text(font = ('Verdana',12), height = 1, width = 20)
gc_content_text.grid(row=2,column=1)

kmer_size_entry = Entry()
kmer_size_entry.grid(row=3,column=1)

complement_text = Text(font = ('Verdana',12), height = 1, width = 20)
complement_text.grid(row=4,column=1)

reading_frames_text = Text(font = ('Verdana',12), height = 6, width = 40)
reading_frames_text.grid(row=5,column=1)

textbox = Text(font = ('Verdana',12), height = 6, width = 40)
textbox.grid(row=7,column=1)


dna_button = Button(root, text='Enter', command = get_sequence)
dna_button.grid(row=6,column=0)


mainloop()

