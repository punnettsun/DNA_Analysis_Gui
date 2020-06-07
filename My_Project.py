# CS_22B_Project
# Author: Punit Sundar
# Last Updated: June 7, 2020
# Purpose: Write an object-oriented program for analyzing DNA sequence,
# including CG content, k-mer counting, complemnentary sequence, and
# translation of all 6 reading frames into protein sequences.

import sys
from tkinter import *
import re

# Fixed
# 1. Incorrect kmer size value
# Still needs fixing
# 1. Delete textbox after consecutive Enters

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

class TooLow(Exception):
    pass
class TooHigh(Exception):
    pass
class Not_DNA_Sequence(Exception):
    pass

def get_sequence():
    sequence = dna_entry.get()
    try:
        match = re.search(r"[^ATCG]", sequence)
        if match:
            raise Not_DNA_Sequence
    except Not_DNA_Sequence:
        textbox.insert(END,'Sequence entered in not a DNA sequence')
        sys.exit()
    try:
        kmer_size = kmer_size_entry.get()
        kmer_size = int(kmer_size)
        if kmer_size < 1:
            raise TooLow
        if kmer_size > 10:
            raise TooHigh
    except ValueError:
        textbox.insert(END,' ValueError. Please type a number only ')
        sys.exit()
    except TooLow:
        textbox.insert(END,' Kmer number too low ')
        sys.exit()
    except TooHigh:
        textbox.insert(END,' High Kmer values might cause slow processing ')
        sys.exit()

    else:
        Seq = Analyze_DNA_Sequence(sequence)
        gc_content = Seq.GC_Content()
        complement = Seq.Complementary('5-3')
        kmer = Seq.kmer_count(int(kmer_size))
        reading_frames = Seq.Translate()
        gc_content_text.insert(1.0,gc_content)
        complement_text.insert(END,complement)
        for keys in reading_frames:
            values = reading_frames[keys]
            final = keys + ':' + str(values) + '\n'
            reading_frames_text.insert(END,final)
        textbox.insert(END,' '+ str(kmer))

def tkinter_open_window():
    root = Tk()
    root.title("DNA to Protein")
    root.configure(background = 'SteelBlue4')

    Directions = Label(text = "Directions: Paste ONLY a DNA sequence and kmer size below. Then click 'Enter'",
                       font = ('Verdana',12,'bold'), bg = 'Salmon',fg = 'white', width = 60, height=4)
    Directions.grid(row = 0, column = 0, columnspan = 2)

    my_list = ['DNA Sequence','GC%','Type Kmer Size',"5'-3' Complement",'Reading Frames']
    for i in range(1,6):
        label = Label(text = my_list[i-1], font = ('Verdana',12,'bold'),bg = 'SteelBlue4', fg = 'white', width = 25, height = 6)
        label.grid(row= i, column=0)
    global dna_entry
    dna_entry = Entry(bg = 'SteelBlue4', fg = 'white', width = 25)
    dna_entry.grid(row = 1, column = 1)

    global gc_content_text 
    gc_content_text = Text(font = ('Verdana',12), height = 1, width = 25, bg = 'SteelBlue4', fg = 'white')
    gc_content_text.grid(row = 2, column = 1)

    global kmer_size_entry
    kmer_size_entry = Entry(bg = 'SteelBlue4', fg = 'white')
    kmer_size_entry.grid(row = 3, column = 1)

    global complement_text
    complement_text = Text(font = ('Verdana',12), height = 5, width = 25, bg = 'SteelBlue4', fg = 'white')
    complement_text.grid(row = 4, column = 1)

    global reading_frames_text
    reading_frames_text = Text(font = ('Verdana',12), height = 6, width = 25, bg = 'SteelBlue4', fg = 'white')
    reading_frames_text.grid(row = 5, column = 1)

    global textbox
    textbox = Text(font = ('Verdana',12), height = 6, width = 40, fg = 'white', bg = 'salmon')
    textbox.grid(row = 7, column = 1)
    textbox.insert(1.0,'')


    dna_button = Button(root, text='Enter', command = get_sequence)
    dna_button.grid(row = 6, column = 0)


    mainloop()

tkinter_open_window()
