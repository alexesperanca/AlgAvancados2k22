# -*- coding: utf-8*-
# Copyright 2022 by Group 7 (MSc Bioinformatics - University of Minho).  All rights reserved.

"""
This module provides the :class:`Automata` class that facilitates the encounter of patterns in sequences.
This class includes diverse strategies to build transition tables and find patterns, such as:
    - Transition table construction that indicates the states transition along the pattern
    - Obtain the next state for a given symbol and current state
    - Find occurrences of the pattern in a given sequence
    - Find the size of the overlap between 2 different sequences
"""

class Automata:
    '''Class that builds a transition table to iterate over for better performance and faster analysis
    '''
    def __init__(self, alphabet: str, pattern: str):
        '''Class initialization to construct the initial transition table for an "alphabet" (chatacters present) and "pattern"

        Parameters
        ----------
        alphabet : str
            Characters present in the pattern
        pattern : str
            Pattern for the transition table construction
        '''
        self.numstates = len(pattern) + 1
        self.alphabet = alphabet
        self.transitionTable = {}
        self.pattern = pattern
        self.buildTransitionTable()   
    
    def buildTransitionTable(self):
        '''Method that constructs the transition table for the given pattern and alphabet previoulsy provided
        '''
        for q in range(self.numstates):
            for a in self.alphabet:
                prefix = self.pattern[0:q] + a
                self.transitionTable[(q, a)] = overlap(prefix, self.pattern)

    def printAutomata(self):
        '''Pretty print the automata data. States, alphabet, and transition table information printed
        '''
        print ("States: " , self.numstates)
        print ("Alphabet: " , self.alphabet)
        print ("Transition table:")
        for k in self.transitionTable.keys():
            print (k[0], ",", k[1], " -> ", self.transitionTable[k])

    def nextState(self, current: int, symbol: str) -> int:
        '''Method that obtains the next state for a give current state "current" and character "symbol"

        Parameters
        ----------
        current : int
            Current state
        symbol : str
            Character of the alphabet

        Returns
        -------
        int
            Next state
        '''
        assert symbol in self.alphabet, "Symbol not found in alphabet"
        return self.transitionTable[(current, symbol)]
        
    def applySeq(self, seq: str) -> list:
        '''Method that returns the states along the given sequence

        Parameters
        ----------
        seq : str
            Sequence to apply and verify the state alteration

        Returns
        -------
        list
            List of alteration of states
        '''
        q = 0
        res = [q]
        for i in seq:
            res.append(self.transitionTable[(q, i)])
            q = self.transitionTable[(q, i)]
        return res
        
    def occurencesPattern(self, text: str) -> tuple:
        '''Method to get the occurences of the pattern in the "text" provided

        Parameters
        ----------
        text : str
            Model to cross the pattern and identify the occurences

        Returns
        -------
        tuple
            Tuple with the occurences of the pattern in the "text"
        '''
        q = 0
        res = []
        for i in range(len(text)):
            q = self.nextState(q, text[i])
            if q == self.numstates - 1:
                res.append(i - self.numstates + 2)
        return res

def overlap(s1: str, s2: str) -> int:
    '''Mathod to get the size of the overlap between two sequences. The biggest sufix of sequence "s1" that is prefix in "s2"

    Parameters
    ----------
    s1 : str
        First sequence
    s2 : str
        Second sequence

    Returns
    -------
    int
        Overlap size between the two sequences
    '''
    maxov = min(len(s1), len(s2))
    for i in range(maxov,0,-1):
        if s1[-i:] == s2[:i]: return i
    return 0
