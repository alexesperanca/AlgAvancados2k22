#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2022 by Group 7 (MSc Bioinformatics - University of Minho).  All rights reserved.

"""
This module provides the `EAMotifsInt` and `EAMotifsReal` subclasses of `EvolAlgorithm` class, a part of a Genetic and Evolutionary Algorithm. 

The `EAMotifsInt` subclass is a binary representation, while the `EAMotifsReal` subclass is a real representation of the Evolutionary Algorithm for founding motifs.

"""

from EvolAlgorithm import EvolAlgorithm
from Popul import PopulInt, PopulReal
from MotifFinding import MotifFinding


def createMatZeros(nl: int, nc: int) -> list:
    '''Method that creates a matrix of zeros given the number of lines and columns.

    Parameters
    ----------
    nl : int
        number of lines
    nc : int
        number of columns

    Returns
    -------
    list
        matrix with 'nl' lines and 'nc' columns filled with zeros
    '''
    res = []
    for _ in range(0, nl):
        res.append([0]*nc)
    return res


def printMat(mat: list):
    '''Method that prints a matrix with better formatting

    Parameters
    ----------
    mat : list
        matrix
    '''
    for i in range(0, len(mat)):
        for j in range(len(mat[i])):
            print(f"{mat[i][j]:.3f}", end=' ')
        print()


class EAMotifsInt (EvolAlgorithm):
    def __init__(self, popsize: int, numits: int, noffspring: int, filename: str) -> None:
        '''SubClass EAMotifsInt, a binary representation of EA algorithm for motif finding. 

        Parameters
        ----------
        popsize : int
            size of population 
        numits : int
            number of iterations to perform
        noffspring : int
            number of new individuals (descendants)
        filename : str
            name of the file to read the sequences, example: 'file.txt'
        '''
        self.motifs = MotifFinding()
        self.motifs.readFile(filename)
        indsize = len(self.motifs)
        EvolAlgorithm.__init__(self, popsize, numits, noffspring, indsize)

    def initPopul(self, indsize: int) -> None:
        '''Method that initializes the population with a given size for individuals (binary representation).

        Parameters
        ----------
        indsize : int
            size of individuals (list of genes)
        '''
        maxvalue = self.motifs.seqSize(0) - self.motifs.motifSize
        self.popul = PopulInt(self.popsize, indsize, maxvalue, [])

    def evaluate(self, indivs: list) -> None:
        '''Method that calculates the score for each individual, setting its fitness.

        Parameters
        ----------
        indivs : list
            list that represents the individuals solution (list of lists with int numbers)
        '''
        for i in range(len(indivs)):
            ind = indivs[i]
            sol = ind.getGenes()
            fit = self.motifs.score(sol)
            ind.setFitness(fit)


class EAMotifsReal (EvolAlgorithm):
    
    def __init__(self, popsize: int, numits: int, noffspring: int, filename: str) -> None:
        '''SubClass EAMotifsReal, a real representation of EA algorithm for motif finding. 

        Parameters
        ----------
        popsize : int
            size of population 
        numits : int
            number of iterations to perform
        noffspring : int
            number of new individuals (descendants)
        filename : str
            name of the file to read the sequences, example: 'file.txt'
        '''
        self.motifs = MotifFinding()
        self.motifs.readFile(filename)
        indsize = len(self.motifs)*len(self.motifs.alphabet)
        EvolAlgorithm.__init__(self, popsize, numits, noffspring, indsize)

    def initPopul(self, indsize: int) -> None:
        '''Method that initializes the population with a given size for individuals (real representation).

        Parameters
        ----------
        indsize : int
            size of individuals (list of genes)
        '''
        maxvalue = self.motifs.seqSize(0) - self.motifs.motifSize
        self.popul = PopulReal(self.popsize, indsize, 0.0, maxvalue, [])

    def profile(self, pwm_list: list) -> list:
        '''Creates a PWM profile of a given list of probabilistic values, with normnalized values (sum of each column(nucleotide) is 1).

        Parameters
        ----------
        pwm_list : list
            List with values

        Returns
        -------
        mat : list
            It returns a list with dictionaries of the probabilistic profile
        '''
        mat = []
        for i in range(0, len(pwm_list), len(self.motifs.alphabet)):
            dic = dict(zip(self.motifs.alphabet, pwm_list[i:]))
            mat.append(dic)
        n_mat = mat
        for n in self.motifs.alphabet:
            soma = 0
            for dici in n_mat:
                soma += dici[n]
            for dici in n_mat:
                dici[n] = dici[n]/soma
        return n_mat, mat

    
    def prob_seq(self, seq: str, profile: list) -> float:
        '''Calculates and returns the probability of a given sequence by the associated profile
        
        Parameters
        ----------
        seq : str
            Sequence of DNA, RNA or AMINO ACIDS
        profile : list
            List with a dictionary of the probabilistic profile (PWM)

        Returns
        -------
        p : float
            A float number representative of the probability of the given sequence
        '''
        assert len(seq) == len(profile), 'Sequence size does not match associated profile size'
        p = 1
        for b, c in zip(seq, profile):
            p *= c[b]
        return round(p, 5)
    
    def seq_most_probable(self, seq: str, profile: list) -> tuple:
        '''Calculates the probability of each subsequence of a given sequence, and returns the subsequence with the highest probability according to the associated profile

        Parameters
        ----------
        seq : str
            Sequence of DNA, RNA or AMINO ACIDS
        profile : list
            List with a dictionary of the probabilistic profile (PWM or PSSM)

        Returns
        -------
        seq : str
            Subsequence with the highest probability
        ind : int
            Index of seq (Subsequence with the highest probability)
        '''
        list_seq = [seq[I:I+len(profile)] for I in range(len(seq) - len(profile) + 1)]
        probs = [self.prob_seq(s, profile) for s in list_seq]
        score_max = max(probs)
        ind = probs.index(score_max)
        seq = [list_seq[I] for I,p in enumerate(probs) if p == score_max][0]
        return seq, ind

    def consensus(self, profile) -> str:
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
        # p_max = 0 
        cons = ''
        # total = len(self.list_seq) + self.n*self.pseudo #só é aplicável ao pwm
        for dic in profile:
            m = max(*dic.values())
            # p_max += int(m*total)
            key = [k for k, v in dic.items() if v == m][0]
            cons += key
        return cons #, p_max

    def evaluate(self, indivs: list) -> None:
        '''Method that builds a pwm profile (normalized) for each individual, determines the most probable position of the motif and calculates the score of the final solution, setting the fitness for that individual. 

        Parameters
        ----------
        indivs : list
            list that represents the individuals (list of lists with float numbers - pwm profile)
        '''
        for i in range(len(indivs)):
            ind = indivs[i]
            sol = ind.getGenes()
            #criar função auxiliar pwm ou pssm (soma de cada coluna é igual a 1 - normalizar valores)
            pwm, o = self.profile(sol) # criar o perfil pwm com a lista anterior (por ordem de 4 em 4 (se abcedario é DNA))
            # para cada sequencia, calcular a posição mais provavel do motif no pwm criado (vetor search)
            search = []
            for x in range(len(self.motifs.seqs)):
                ## fazer função auxiliar para seq mais provavel (nao dá para utilizar função da Classe Motifs porque não é criada uma instancia, fazemos o perfil direto da lista)
                seq_prob, index = self.seq_most_probable(self.motifs.seqs[x], pwm)
                search.append(index)

            # calcular o score
            fit = self.motifs.score(search)
            ind.setFitness(fit)



