#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2022 by Group 7 (MSc Bioinformatics - University of Minho).  All rights reserved.

"""
This module provides the `Popul` class and its subclasses `PopulInt` and `PopulReal`, for creating populations of individuals (needed for Evolutionary Algorithm: `EvolAlgorithm` class).
This class includes a binary and a real representation of the population.

"""

from Indiv import Indiv, IndivInt, IndivReal
from random import random


class Popul:

    def __init__(self, popsize: int, indsize: int, indivs: list = []) -> None:
        '''Class that implements the population with a given size.

        Parameters
        ----------
        popsize : int
            size of population 
        indsize : int
            size of individuals (list of genes)
        indivs : list, optional
            list of genes, by default []
        '''
        self.popsize = popsize
        self.indsize = indsize
        if indivs:
            self.indivs = indivs
        else:
            self.initRandomPop()

    def getIndiv(self, index: int) -> list:
        '''Method that returns a specific individual from the list, given its index.

        Parameters
        ----------
        index : int
            index of the individual on the list self.indivs

        Returns
        -------
        list
            list of information about individual: list of genes and fitness (Instance of the class Indiv)
        '''
        return self.indivs[index]

    def initRandomPop(self) -> None:
        '''Method that initializes the population (creates instances of Indiv class)
        '''
        self.indivs = []
        for _ in range(self.popsize):
            indiv_i = Indiv(self.indsize, [])
            self.indivs.append(indiv_i)

    def getFitnesses(self, indivs: list = None) -> list:
        '''Method that returns the fitness information about all the individuals

        Parameters
        ----------
        indivs : list, optional
            list of individuals, by default None

        Returns
        -------
        list
            list of fitness of all the individuals
        '''
        fitnesses = []
        if not indivs:
            indivs = self.indivs
        for ind in indivs:
            fitnesses.append(ind.getFitness())
        return fitnesses

    def bestSolution(self) -> list:
        '''Method that returns the best solution of all the individuals (that has the best fitness).

        Returns
        -------
        list
            the best solution
        '''
        return max(self.indivs)

    def bestFitness(self) -> int | float:
        '''Method that returns the best fitness of all the individuals

        Returns
        -------
        int | float
            the best fitness 
        '''
        indv = self.bestSolution()
        return indv.getFitness()


    def selection(self, n: int, indivs: list = None) -> list:
        '''_summary_

        Parameters
        ----------
        n : int
            size of selection (parents)
        indivs : _type_, optional
            list that represents the individuals, by default None

        Returns
        -------
        list
            selection list
        '''
        res = []
        fitnesses = list(self.linscaling(self.getFitnesses(indivs)))
        for _ in range(n):
            sel = self.roulette(fitnesses)
            fitnesses[sel] = 0.0
            res.append(sel)
        return res

    def roulette(self, f):
        tot = sum(f)
        val = random()
        acum = 0.0
        ind = 0
        while acum < val:
            acum += (f[ind] / tot)
            ind += 1
        return ind-1

    def linscaling(self, fitnesses):
        mx = max(fitnesses)
        mn = min(fitnesses)
        res = []
        for f in fitnesses:
            val = (f-mn)/(mx-mn)
            res.append(val)
        return res

    def recombination(self, parents, noffspring):
        offspring = []
        new_inds = 0
        while new_inds < noffspring:
            parent1 = self.indivs[parents[new_inds]]
            parent2 = self.indivs[parents[new_inds+1]]
            offsp1, offsp2 = parent1.crossover(parent2)
            offsp1.mutation()
            offsp2.mutation()
            offspring.append(offsp1)
            offspring.append(offsp2)
            new_inds += 2
        return offspring

    def reinsertion(self, offspring):
        tokeep = self.selection(self.popsize-len(offspring))
        ind_offsp = 0
        for i in range(self.popsize):
            if i not in tokeep:
                self.indivs[i] = offspring[ind_offsp]
                ind_offsp += 1


class PopulInt(Popul):

    def __init__(self, popsize, indsize, ub, indivs=[]):
        self.ub = ub
        Popul.__init__(self, popsize, indsize, indivs)

    def initRandomPop(self):
        self.indivs = []
        for _ in range(self.popsize):
            indiv_i = IndivInt(self.indsize, [], 0, self.ub)
            self.indivs.append(indiv_i)


class PopulReal(Popul):
  
    def __init__(self, popsize, indsize, lb=0.0, ub=1.0, indivs=[]):
        self.ub = ub
        self.lb = lb
        Popul.__init__(self, popsize, indsize, indivs)

    def initRandomPop(self):
        self.indivs = []
        for _ in range(self.popsize):
            indiv_r = IndivReal(self.indsize, [], self.lb, self.ub)
            self.indivs.append(indiv_r)
