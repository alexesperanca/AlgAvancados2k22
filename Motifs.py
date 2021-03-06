#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2022 by Group 7 (MSc Bioinformatics - University of Minho).  All rights reserved.

"""
This module provides the :class:`Motifs` class, for creating probabilistic profiles given a list of sequences.
This class includes diverse methods that allows exploring these, such as:
    - prob_seq(): returns the probability of a given sequence based on the profile originally created.
    - seq_most_probable(): returns the subsequence with the highest probability according to the associated profile.
    - consensus(): returns the consensus between two different sequences using the profile created.
"""

from Sequence import Sequence
import math


class Motifs():
    def __init__(self, list_seq: list, **kwargs) -> None:
        '''Class that allows creating and displaying probabilistic profiles (PWM or PSSM) and calculate the probability of a sequence and calculate the most probable sequence of the profile created.

        Parameters
        ----------
        list_seq : list
            Sequences of DNA, RNA or AMINO ACIS
        
        kwargs: dictionary
            May take as keys:
                profile -- pseudo
        '''
        self.list_seq = list_seq
        self.abc = self.validate_abc()
        self.profile_type = 'pwm'
        self.pseudo = 0
        self.n = 4
        
        if self.abc == 'DNA':
            self.abc_letters = "ACTG"
        elif self.abc == 'AMINO':
            self.abc_letters = "FLIMVSPTAYHQNKDECWRG_"
            self.n = 20
        elif self.abc == 'RNA':
            self.abc_letters = "ACGU"
        
        if 'profile_type' in kwargs:
            self.profile_type = kwargs['profile_type']
            
        if 'pseudo' in kwargs:
            self.pseudo = kwargs['pseudo']

        self.profile = self.create_profile()
        
    def validate_seq(self, seq: str) -> str:
        '''This method determines the type of the sequence: DNA, RNA or AMINO ACIDS
        It gives an error if the sequence is not recognized.

        Parameters
        ----------
        seq : str
            Sequence of DNA, RNA or AMINO ACIDS

        Returns
        -------
        str
            Type of the sequence ('DNA', 'RNA' or 'AMINO')
        '''
        test = Sequence(seq)
        return test.check
        
    
    def validate_abc(self) -> str:
        '''Through the Sequences Class, it validates if all the sequences of a given list are of the same type, be it DNA, RNA or PROTEIN (AMINO).
        If any sequence is not recognized or of the same type, an Error is displayed

        Returns
        -------
        str
            Sequence of DNA, RNA or AMINO ACIDS

        Raises
        ------
        TypeError
            None of the sequences are of the same type
        '''
        class_abc = Sequence(self.list_seq[0])
        abc = class_abc.check
        for i in self.list_seq[1:]:           
            valid = self.validate_seq(i)
            if valid != abc:
                raise TypeError("Sequences introduced are not of the same type!") 
        return abc
    
    
    def print_profile(self) -> repr:
        '''Prints a probabilistic profile legibly            

        Returns
        -------
        repr
            A pretty way of showing the probabilistic profile
        '''
        bases = sorted(self.profile[0].keys())
        tab = [[f"{p[b]:-5.2f}" for b in bases] for p in self.profile]
        for p in zip(*([bases] + tab)):
            print(*p)
            
    
    def __create_prof__(self) -> list:
        '''Auxiliary function of create_profile() that calculates the frequency of each position

        Returns
        -------
        list
            It returns list with dictionary of the profile
        '''
        total = len(self.list_seq) + self.n*self.pseudo
        matrix = []
        for line in zip(*self.list_seq):
            matrix.append({k: (line.count(k)+self.pseudo) for k in self.abc_letters})
        return matrix
    
    def create_profile(self) -> list:
        '''Calculates the probabilistic PWM (Position Weighted Matrix) or PSSM (Position Specific Scoring Matrix) profile of a given list of **DNA** sequences.
        Default pseudocount of 0

        Returns
        -------
        list
            It returns list with a dictionary of the probabilistic profile (PWM or PSSM)

        Raises
        ------
        TypeError
            If the profile selected isn't pwm or pssm
        '''
        total = len(self.list_seq) + self.n*self.pseudo
        matrix = self.__create_prof__()
        if self.profile_type == 'pwm':
            for dic in matrix:
                for i in dic:
                    dic[i]= dic[i]/total
        elif self.profile_type == 'pssm':
            for dic in matrix:
                for i in dic:
                    dic[i] = math.log2((dic[i]/total)/(1/self.n))
        else:
            raise TypeError("Type of profile invalid")   
        return matrix

    def prob_seq(self, seq: str) -> float:
        '''Calculates and returns the probability of a given sequence by the associated profile
        
        Parameters
        ----------
        seq : str
            Sequence of DNA, RNA or AMINO ACIDS
        profile : list
            List with a dictionary of the probabilistic profile (PWM or PSSM)

        Returns
        -------
        float
            A float number representative of the probability of the given sequence

        Raises
        ------
        TypeError
            If the profile selected isn't pwm or pssm
        '''
        assert len(seq) == len(self.profile), 'Sequence size does not match associated profile size'
        if self.profile_type == 'pwm':
            p = 1
            for b, c in zip(seq, self.profile):
                p *= c[b]
        elif self.profile_type == 'pssm':
            p = 0
            for b, c in zip(seq, self.profile):
                p += c[b]
        else:
            raise TypeError("Type of profile invalid")
        return round(p, 5)
    
    def seq_most_probable(self, seq: str) -> str:
        '''Calculates the probability of each subsequence of a given sequence, and returns the subsequence with the highest probability according to the associated profile

        Parameters
        ----------
        seq : str
            Sequence of DNA, RNA or AMINO ACIDS
        profile : list
            List with a dictionary of the probabilistic profile (PWM or PSSM)

        Returns
        -------
        str
            Subsequence with the highest probability
        '''
        list_seq = [seq[I:I+len(self.profile)] for I in range(len(seq) - len(self.profile) + 1)]
        probs = [self.prob_seq(s) for s in list_seq]
        score_max = max(probs)
        ind = probs.index(score_max)
        seq = [list_seq[I] for I,p in enumerate(probs) if p == score_max][0]
    
        return seq, ind


    def consensus(self) -> str:
        '''Creates a sequence consensus between two different sequences using the profile created.

        Parameters
        ----------
        s1 : str
            String sequence of DNA or Amino Acids
        s2 : str
            String sequence of DNA or Amino Acids

        Returns
        -------
        str
            Sequence consensus
        '''
        p_max = 0 
        cons = ''
        total = len(self.list_seq) + self.n*self.pseudo #s?? ?? aplic??vel ao pwm
        for dic in self.profile:
            m = max(*dic.values())
            p_max += int(m*total)
            key = [k for k, v in dic.items() if v == m][0]
            cons += key
        return cons #, p_max
