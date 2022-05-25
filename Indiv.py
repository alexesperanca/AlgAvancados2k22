#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2022 by Group 7 (MSc Bioinformatics - University of Minho).  All rights reserved.

"""
This module provides the :class:`Indiv` class, that implements individuals with binary (IndivInt subclass) or real representations (IndivReal subclass).
This class includes diverse methods, such as: mutation(), crossover(), and others.
"""

from random import randint, random, shuffle, uniform

class Indiv:

    def __init__(self, size: int, genes: list = [], lb: int = 0, ub: int = 1) -> None:
        '''Class to implement individuals with binary and real representations, with the atributes:
        - lb/ub (lower and upper limits of the range for representing genes)
        - fitness (stores fitness value for each individual)

        Parameters
        ----------
        size : int
            size of the list of genes
        genes : list, optional
            list of genes (representative of genome), by default []
        lb : int, optional
            lower limits of the range for representing genes, by default 0
        ub : int, optional
            upper limits of the range for representing genes, by default 1
        '''
        self.lb = lb
        self.ub = ub
        self.genes = genes
        self.fitness = None
        if not self.genes: # random generation of list of genes:
            self.initRandom(size)

    # comparadores.
    # Permitem usar sorted, max, min

    def __eq__(self, solution: list) -> bool:
        '''Method auxiliary that determines if the solution is an instance of Indiv class.

        Parameters
        ----------
        solution : list
            individual (instance of class if true)

        Returns
        -------
        bool
            True if solution is instance of this class
        '''
        if isinstance(solution, self.__class__):
            return self.genes.sort() == solution.genes.sort()
        return False

    def __gt__(self, solution: list) -> bool:
        '''Method auxiliary that determines if the current fitness is greater than the solutions' fitness.

        Parameters
        ----------
        solution : list
            individual (instance of class if true)

        Returns
        -------
        bool
            True if current fitness is greater than its solution
        '''
        if isinstance(solution, self.__class__):
            return self.fitness > solution.fitness
        return False

    def __ge__(self, solution: list) -> bool:
        '''Method auxiliary that determines if the current fitness is greater or equal than the solutions' fitness.

        Parameters
        ----------
        solution : list
            individual (instance of class if true)

        Returns
        -------
        bool
            True if current fitness is greater or equal than its solution
        '''
        if isinstance(solution, self.__class__):
            return self.fitness >= solution.fitness
        return False

    def __lt__(self, solution: list) -> bool:
        '''Method auxiliary that determines if the solution fitness is greater than the current.

        Parameters
        ----------
        solution : list
            individual (instance of class if true)

        Returns
        -------
        bool
            True if the solution fitness is greater than the curent fitness
        '''
        if isinstance(solution, self.__class__):
            return self.fitness < solution.fitness
        return False

    def __le__(self, solution: list) -> bool:
        '''Method auxiliary that determines if the solution fitness is greater or equal than the current.

        Parameters
        ----------
        solution : list
            individual (instance of class if true)

        Returns
        -------
        bool
            True if the solution fitness is greater or equal than the curent fitness
        '''
        if isinstance(solution, self.__class__):
            return self.fitness <= solution.fitness
        return False

    def __str__(self) -> str:
        ''' Method that writes the information about the object: the list of genes and its fitness
        
        Returns
        ------------
        str
            String with the information of the object
        '''
        return f"{str(self.genes)} {self.getFitness()}"

    def __repr__(self) -> str:
        ''' Method that shows the information about the object
        
        Returns
        ------------
        str
            String with the information of the object
        '''
        return self.__str__()

    def setFitness(self, fit: int | float) -> None:
        '''Method that sets the new fitness of the individual 

        Parameters
        ----------
        fit : int
            fitness of the individual
        '''
        self.fitness = fit

    def getFitness(self) -> int | float:
        '''Method that shows the fitness of the individual 

        Returns
        -------
        int
            fitness of the individual
        '''
        return self.fitness

    def getGenes(self) -> list:
        '''Method that shows the list of genes of the individual 

        Returns
        -------
        list
            list of genes of the individual
        '''
        return self.genes

    def initRandom(self, size: int) -> None:
        '''Method that generates a list of genes of the individual (random int numbers between upper and lower bounds) 

        Parameters
        ----------
        size : int
            number of genes to generate
        '''
        self.genes = []
        for _ in range(size):
            self.genes.append(randint(self.lb, self.ub))

    def mutation(self) -> None:
        '''Method for binary representations that alters a single gene (mutation)
        '''
        s = len(self.genes)
        pos = randint(0, s-1)
        if self.genes[pos] == 0:
            self.genes[pos] = 1
        else:
            self.genes[pos] = 0

    def crossover(self, indiv2) -> tuple:
        '''Method that makes a crossover between two individuals

        Parameters
        ----------
        indiv2 : type[Indiv]
            _description_

        Returns
        -------
        tuple
            _description_
        '''
        return self.one_pt_crossover(indiv2)

    def one_pt_crossover(self, indiv2: list) -> tuple:
        '''Auxiliary method that makes a crossover between two individuals

        Parameters
        ----------
        indiv2 : list
            individual (instance of Indiv class)

        Returns
        -------
        tuple
            the new individuals with list of genes crossed-over 
        '''
        offsp1 = []
        offsp2 = []
        s = len(self.genes)
        pos = randint(0, s-1)
        for i in range(pos):
            offsp1.append(self.genes[i])
            offsp2.append(indiv2.genes[i])
        for i in range(pos, s):
            offsp2.append(self.genes[i])
            offsp1.append(indiv2.genes[i])
        res1 = self.__class__(s, offsp1, self.lb, self.ub)
        res2 = self.__class__(s, offsp2, self.lb, self.ub)
        return res1, res2


class IndivInt (Indiv):

    def __init__(self, size: int, genes: list = [], lb: int = 0, ub: int = 1) -> None:
        '''Subclass to implement individuals with binary representation.

        Parameters
        ----------
        size : int
            size of the list of genes
        genes : list, optional
            list of genes (representative of genome), by default []
        lb : int, optional
            lower limits of the range for representing genes, by default 0
        ub : int, optional
            upper limits of the range for representing genes, by default 1
        '''
        self.lb = lb
        self.ub = ub
        self.genes = genes
        self.fitness = None
        if not self.genes:
            self.initRandom(size)

    def initRandom(self, size: int) -> None:
        '''Method that generates a list of genes of the individual (random int numbers between upper and lower bounds) 

        Parameters
        ----------
        size : int
            number of genes to generate
        '''
        self.genes = []
        for _ in range(size):
            self.genes.append(randint(self.lb, self.ub))


class IndivReal(Indiv):
    
    def __init__(self, size: int, genes: list = [], lb: float = 0.0, ub: float = 1.0) -> None:
        '''Subclass to implement individuals with real representation.

        Parameters
        ----------
        size : int
            size of the list of genes
        genes : list, optional
            list of genes (representative of genome), by default []
        lb : int, optional
            lower limits of the range for representing genes, by default 0.0
        ub : int, optional
            upper limits of the range for representing genes, by default 1.0
        '''
        self.lb = lb
        self.ub = ub
        self.genes = genes
        self.fitness = None
        if not self.genes:
            self.initRandom(size)

    def initRandom(self, size: int) -> None:
        '''Method that generates a list of genes of the individual (random float numbers between upper and lower bounds) 

        Parameters
        ----------
        size : int
            number of genes to generate
        '''
        self.genes = []
        for _ in range(size):
            self.genes.append(uniform(self.lb, self.ub))

    def mutation(self) -> None:
        '''Method for real representations that alters a single gene (mutation)
        '''
        s = len(self.genes)
        pos = randint(0, s-1)
        self.genes[pos] = uniform(self.lb, self.ub)
