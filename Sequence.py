#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2022 by Group 7 (MSc Bioinformatics - University of Minho).  All rights reserved.

"""
This module provides the :class:`Sequence` class, for verifying sequences and determing the type (DNA, RNA or PROTEIN).
This class includes diverse methods, including: transcription, comp_inverse, transcript, get_orfs, translation, get_all_prots, and others.
"""

import re

class Sequence:
    def __init__(self, seq: str) -> None:
        '''
        Verifies and capitalizes the sequence inputted by the user implementing an error in case the sequence submitted is invalid
        
        Parameters
        ------------
        seq: str
            Sequence inputted by the user
        '''
        self.seq = seq.upper()
        dna = "ACGT-"
        rna = "ACGU-"
        amino = "FLIMVSPTAY_HQNKDECWRG-"
        
        if all (i in dna for i in self.seq) is True:
            self.check = "DNA"
        elif all (i in rna for i in self.seq) is True:
            self.check = "RNA"
        elif all (i in amino for i in self.seq) is True:
            self.check = "AMINO"
        else:
            raise TypeError("Invalid input string!")
            
    def __str__(self) -> str:
        '''
        Writing the string with the sequence and its type
        
        Returns
        ------------
        str
            String with the information of the object
        '''
        return self.check + ":" + self.seq
        
    def print_seq(self) -> str:
        '''
        Printing the complete sequence
        
        Returns
        ------------
        str
            String of the complete sequence
        '''
        print(self.seq)

    def alphabet(self) -> str:
        if self.check == "DNA":
            return 'ACGT'
        if self.check == "RNA":
            return 'ACGU'
        if self.check == "amino":
            return 'FLIMVSPTAY_HQNKDECWRG'

    
    def get_slice(self, i: str, j: str) -> str:
        '''
        Slice of sthe sequence between i and j
        
        Returns
        ------------
        str
            String with the slice of the sequence
        '''
        return self.seq[i:j]
    	
    def percentage(self) -> dict:
        '''
        Construction of a dictionary with the elements of the sequence and corresponding percentage
        
        Returns
        ------------
        dict
            Dictionary with the characters and percentage
        '''
        count = {c: 0 for c in self.seq}
        for i in self.seq:
            count[i] += 1
        
        for c, v in count.items(): 
            count[c] = round(v/len(self.seq), 3)*100
            
        return count
    
    def ind_percentage(self, n: str) -> str:
        '''
        Iterates over the dictionary with the characters and percentages
        
        Returns
        ------------
        str
            Percentage of the character introduced in the object sequence
        '''
        count = Sequence.percentage(self)
        return f"{count[n.upper()]}%"
        
    
    def transcription(self) -> str:
        '''
        Replacing of the thymine base 'T' for the uracil base 'U'
        Only executes in case the object is a DNA sequence
        
        Parameters
        ----------
        str
            DNA sequence
        
        Returns
        ------------
        str
            RNA sequence
        '''
        
        assert self.check == "DNA", "Introduced sequence must be DNA!"
        if self.seq:
        
            return self.seq.upper().replace("T", "U")
        
        else:
             raise Exception()
    
    
    def comp_inverse(self) -> str:
        '''
        Inversion and complementarity of the sequence
        Only executes in case the object is a DNA sequence
        
        Parameters
        ----------
        seq: str
            DNA sequence
        
        Returns
        ------------
        str
            String corresponding to the inverse and complement of the original object sequence
        '''
        assert self.check == "DNA", "Introduced sequence must be DNA!"
        
        if self.seq:
            change = {"A": "T", "T": "A", "G": "C", "C": "G"}
            new = "".join(change[i] for i in self.seq)
        
        else:
            raise Exception() #aqui devia ter mensagem do erro do exception, não devia estar vazio
        
        return new[::-1]
    
    def cDNA(self) -> str:
        '''
        Construction of a dictionary with all the characters and corresponding replaces and then implementing on a new string
        Only executes in case the object is a RNA sequence
        
        Parameters
        ----------
        seq: str
            RNA sequence  
        
        Returns
        ------------
        str
            Returns the string that represents the complementary DNA of the original object sequence
        '''
        assert self.check == "RNA", "Introduced sequence must be RNA!"
        
        change = {"A": "T", "U": "A", "G": "C", "C": "G"}
        new = "".join(change[i] for i in self.seq)
        
        return new
    
    def get_orfs(self) -> dict:
        '''
        Construction of a dictionary with all the 6 orfs (3 of each chain)
        Only executes in case the object is a DNA sequence
        
        Parameters
        ----------
        seq: str
            DNA sequence of which we want to know the reading frames (orfs)
               
        Returns
        ------------
        dict
            Returns the dictionary with each orf and the corresponding list of codons
        '''
        assert self.check == "DNA", "Introduced sequence must be DNA!"
        seq_rev = Sequence.comp_inverse(self)
        
        orfs = {"ORF +1": None, "ORF +2": None, "ORF +3": None, "ORF -1": None, "ORF -2": None, "ORF -3": None}
    
        orfs["ORF +1"] = re.findall("(...)", self.seq)
        orfs["ORF +2"] = re.findall("(...)", self.seq[1:])
        orfs["ORF +3"] = re.findall("(...)", self.seq[2:])
        orfs["ORF -1"] = re.findall("(...)", seq_rev)
        orfs["ORF -2"] = re.findall("(...)", seq_rev[1:])
        orfs["ORF -3"] = re.findall("(...)", seq_rev[2:])
            
        return orfs
    
    def get_codons(self) -> list:
        '''
        Obtaining of all the codons for the ORF +1
        Only executes in case the object is a DNA sequence
        
        Parameters
        ----------
        seq: str
            DNA sequence
        
        Returns
        ------------
        list
            List with all the codons of the sequence
        '''
        assert self.check == "DNA", "Introduced sequence must be DNA!"
        return re.findall("(...)", self.seq) #aqui devia dar para fazer codoes das outras sequencias (invertidas e assim)
    
    def translation(self) -> str:
        '''
        Obtaining the sequence codons and building an amino acid chain across the codons
        Transforms codons into amino acids through a dictionary with the relevant information
        Only executes in case the object is a DNA sequence
        
        Parameters
        ----------
        seq: str
            DNA sequence
        
        Returns
        ------------
        str
            Amino acid chain string resulting from the translation of the DNA sequence
        '''
        assert self.check == "DNA", "Introduced sequence must be DNA!"
        amino_codons = {"F":["TTT", "TTC"], "L":["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"], "I":["ATT", "ATC", "ATA"], "M":"ATG", "V":["GTT", "GTC", "GTA", "GTG"], "S":["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"], "P":["CCT", "CCC", "CCA", "CCG"], "T":["ACT", "ACC", "ACA", "ACG"], "A":["GCT", "GCC", "GCA", "GCG"], "Y":["TAT", "TAC"], "_":["TAA", "TAG", "TGA"], "H":["CAT", "CAC"], "Q":["CAA", "CAG"], "N":["AAT", "AAC"], "K":["AAA", "AAG"], "D":["GAT", "GAC"], "E":["GAA", "GAG"], "C":["TGT", "TGC"], "W":"TGG", "R":["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"], "G":["GGT", "GGC", "GGA", "GGG"]}
        aa = ""
        codons = Sequence.get_codons(self) 
        
        if self.seq:
            for codon in codons:
                for amino, codoes in amino_codons.items():
                    if codon in codoes:
                        aa += amino
        
            return aa
        
        else:
             raise Exception()
    
    
    def get_all_prots(self) -> list:
        '''
        Construction of a list with all the existing and possible proteins in the sequence. The sequence can be proteic or DNA
        Only executes in case the object is a DNA or AMINO sequence
        
        Parameters
        ----------
        seq : str
            DNA or amino sequence

        Returns
        ------------
        list
            Returns a list with all the proteins in the sequence
        '''
        assert self.check == "DNA" or self.check == "AMINO", "Introduced sequence must be DNA or AMINO!"
        
        list_prots = []
        
        if self.check == "DNA":
            aa = Sequence.translation(self)
            boolean = False
            proteina = ""
            for i in aa:
                if i == "M": #Aqui como repete em baixo, devia ser uma função auxiliar que é chamada nos dois sitios
                    proteina += i
                    boolean = True
                elif i == "_" and boolean is True:
                    proteina += i
                    boolean = False
                    list_prots.append(proteina)
                    proteina = ""
                elif boolean is True:
                    proteina += i
        
        elif self.check == "AMINO":
            boolean = False
            proteina = ""
            for i in self.seq:
                if i == "M":
                    proteina += i
                    boolean = True
                elif i == "_" and boolean is True:
                    proteina += i
                    boolean = False
                    list_prots.append(proteina)
                    proteina = ""
                elif boolean is True:
                    proteina += i
                
        return list_prots
        
    def get_bigger_prot(self) -> str:
        '''
        Getting the bigger protein in the sequence
        
        Returns
        ------------
        str
            String with the bigger protein
        '''
        list = self.get_all_prots()
        bigger = ""
        for p in list:
            if len(p)>len(bigger):
                bigger = p
        return bigger
    
    def get_aa_orfs(self) -> dict:
        '''
        Construction of a dictionary with all the ORFs and the corresponding chain of amino acids
        Transforms codons into amino acids through a dictionary with the relevant information
        Only executes in case the object is a DNA sequence

        Parameters
        ----------
        seq: str
            DNA sequence     

        Returns
        ------------
        dict
            Returns a dictionary with the ORFs and amino acid chain
        '''
        assert self.check == "DNA", "Introduced sequence must be DNA!"
        orfs = Sequence.get_orfs(self)
        amino_codons = {"F":["TTT", "TTC"], "L":["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"], "I":["ATT", "ATC", "ATA"], "M":"ATG", "V":["GTT", "GTC", "GTA", "GTG"], "S":["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"], "P":["CCT", "CCC", "CCA", "CCG"], "T":["ACT", "ACC", "ACA", "ACG"], "A":["GCT", "GCC", "GCA", "GCG"], "Y":["TAT", "TAC"], "_":["TAA", "TAG", "TGA"], "H":["CAT", "CAC"], "Q":["CAA", "CAG"], "N":["AAT", "AAC"], "K":["AAA", "AAG"], "D":["GAT", "GAC"], "E":["GAA", "GAG"], "C":["TGT", "TGC"], "W":"TGG", "R":["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"], "G":["GGT", "GGC", "GGA", "GGG"]}
        
        orfs_aa = {"ORF +1": None, "ORF +2": None, "ORF +3": None, "ORF -1": None, "ORF -2": None, "ORF -3": None}
        for c, v in orfs.items():
            aa = ""
            for codon in v:
                for amino, codons in amino_codons.items():
                    if codon in codons: aa += amino
            orfs_aa[c] = aa
        
        return orfs_aa
   
    
    def get_all_prots_orfs(self) -> dict:
        '''
        Construction of a dictionary with all the ORFs and the corresponding list with all the existing and possible proteins in the sequence ORF
        Only executes in case the object is a DNA sequence
        
        Parameters
        ----------
        seq: str
            DNA sequence  
       
        Returns
        ------------
        dict
            Returns a dictionary with the ORFs and possible proteins
        '''
        assert self.check == "DNA", "Introduced sequence must be DNA!"
        orfs_aa = Sequence.get_aa_orfs(self)
        
        orfs_prot = {"ORF +1": None, "ORF +2": None, "ORF +3": None, "ORF -1": None, "ORF -2": None, "ORF -3": None}
        for c, aa in orfs_aa.items():
            boolean = False
            proteina = ""
            list_prots = []
            for i in aa:
                if i == "M": #outra vez o mesmo codigo que está acima, devia ter a tal 
                    proteina += i
                    boolean = True
                elif i == "_" and boolean is True:
                    proteina += i
                    boolean = False
                    list_prots.append(proteina)
                    proteina = ""
                elif boolean is True:
                    proteina += i
                    
            orfs_prot[c] = list_prots
            
        return orfs_prot



# teste
def teste():
    # seq_dna = input("Sequencia:")
    seq_dna = 'ATTTAATTACAAGTCTTCAGAATGCCAGAGATATACAGGATCTAACCAATGTTAGAGTTATTAAAAAGTCTGGTATTCGCCGTAATCATGGTACCTGTCGTGATGGCCATCATCCTGGGTCTGATTTACGGTCTTGGTGAAGTATTCAACATCTTTTCTGGTGTTGGTAAAAAAGACCAGCCCGGACAAAATCATTGA'
    s1 = Sequence(seq_dna)
    s1.print_seq()

    print("Transcricao: ")
    print(s1.transcription())
    print("Complemento inverso:") 
    print(s1.comp_inverse())
    print("Traducao: ") 
    print(s1.translation())
    print("ORFs:")
    print(s1.get_aa_orfs())
    print("Proteinas:")
    print(s1.get_all_prots())
    print("Maior proteina:")
    print(s1.get_bigger_prot()) 

if __name__ == "__main__": 
    teste()
