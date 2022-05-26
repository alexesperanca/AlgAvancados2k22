# -*- coding: utf-8 -*-
# Copyright 2022 by Group 7 (MSc Bioinformatics - University of Minho).  All rights reserved.

"""
This module provides the :class:`BoyerMoore` class.
This class includes diverse strategies to apply different methods to facilitate pattern matching in a text or sequence, such as:
    - Bad Character Rule and Good Suffix Rule process of a given pattern
    - Search matches of pattern in the provided text/sequence applying the BCR and GSR methods 
"""

class BoyerMoore:
    def __init__(self, alphabet: str, pattern: str):
        '''Class initialization with pre-process started

        Parameters
        ----------
        alphabet : str
            Possible characters present in pattern and future provided texts/sequences
        pattern : str
            String to search for matches in further text/sequences provided in functions 
        '''
        self.alphabet = alphabet.upper()
        self.pattern = pattern.upper()
        self.preprocess()

    def preprocess(self):
        '''Bad Character Rule (BCR) and Good Suffix Rule (GSR) process initiated
        '''
        self.process_bcr()
        self.process_gsr()
        
    def process_bcr(self) -> dict:
        '''Method Bad Character Rule that constructs a dictionary with all the elements of the alphabet and respective steps to advance if pattern doed not match the text/sequence

        Returns
        -------
        dict
            Dictionary with all the elements of the alphabet and respective steps to advance
        '''
        self.occ = {self.alphabet[s]: -1 for s in range(len(self.alphabet))}
        for j in range(0, len(self.pattern)):
            for c in self.occ:
                if self.pattern[j] in self.occ.keys():
                    self.occ[c] = j
        return self.occ

    def process_gsr(self) -> list:
        '''Method Good Suffix Rule that advances according to the suffix of the match obtained. It moves to the next occurrence of the pattern in the match before failure.

        Returns
        -------
        list
            List of movements to advance according to the size of match obtained
        '''
        self.f = [0 for value in range(len(self.pattern) +1)]
        self.s = [0 for value in range(len(self.pattern) +1)]
        i = len(self.pattern)
        j = len(self.pattern) + 1
        self.f[i] = j
        while i > 0:
            while j <= len(self.pattern) and self.pattern[i-1] != self.pattern[j-1]:
                if self.s[j] == 0:
                    self.s[j] = j - i
                j = self.f[j]
            i = i - 1
            j = j - 1
            self.f[i] = j
        j = self.f[0]
        for i in range(0, len(self.pattern)):
            if self.s[i] == 0:
                self.s[i] = j
            if i == j:
                j = self.f[j]

        return self.s
        
    def search_pattern(self, text: str) -> list:
        '''Method that searches the perfect match of the pattern in the given text by executing the Bad Character Rule and Good Suffix Rule 

        Parameters
        ----------
        text : str
            Text to find the place where the pattern matches

        Returns
        -------
        list
            List of positions that the pattern matched in the text
        '''
        res = []
        i = 0
        text=text.upper()
        while i <= (len(text)-len(self.pattern)):
            j = len(self.pattern) - 1
            while j >= 0 and self.pattern[j] == text[j+i]:
                j = j - 1
            if j < 0:
                res.append(i)
                i = i + self.s[0]
            else:
                c = text[j+i]
                i += max(self.s[j+1], j-self.occ[c])

        if len(res) == 0:
            return 'No match!'
        else: 
            return res

def test():
    padrao = "ACCA"
    texto = "ATAGAACCAATGAACCATGATGAACCATGGATACCCAACCACC"
    bm = BoyerMoore('ACTG', padrao)
    print (bm.process_bcr())
    print (bm.process_gsr())
    print (bm.search_pattern(texto))

if __name__ == '__main__':
    test()

# result: [5, 13, 23, 37]