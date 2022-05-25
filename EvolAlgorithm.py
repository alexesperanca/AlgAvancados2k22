#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2022 by Group 7 (MSc Bioinformatics - University of Minho).  All rights reserved.

"""
This module provides the `EvolAlgorithm` class, a Genetic and Evolutionary Algorithm. 
This work with populations, where each individual represents a solution to a given problem.
Solutions generally encoded in a sequence of symbols of a given alphabet (genome) that allow the application of operators.
New solutions created from mutation operators or sexual recombination.

"""

from Popul import Popul

class EvolAlgorithm:

    def __init__(self, popsize: int, numits: int, noffspring: int, indsize: int) -> None:
        '''Class that implements the Evolutionary Algorithm.

        Parameters
        ----------
        popsize : int
            size of population 
        numits : int
            number of iterations to perform
        noffspring : int
            number of new individuals (descendants)
        indsize : int
            size of individuals (list of genes)
        '''
        self.popsize = popsize
        self.numits = numits
        self.noffspring = noffspring
        self.indsize = indsize

    def initPopul(self, indsize: int) -> None:
        '''Method that initializes the population with a given size for individuals.

        Parameters
        ----------
        indsize : int
            size of individuals (list of genes)
        '''
        self.popul = Popul(self.popsize, indsize)

    def evaluate(self, indivs: list) -> None:
        '''Method that calculates the score for each individual, setting its fitness.

        Parameters
        ----------
        indivs : list
            list that represents the individuals solution. 
        '''
        for i in range(len(indivs)):
            ind = indivs[i]
            fit = 0.0
            for x in ind.getGenes():
                if x == 1:
                    fit += 1.0
            ind.setFitness(fit)
        return None

    def iteration(self) -> None:
        '''Method auxiliary of the Evolutionary Algorithm cycle (based on the number of iterations wanted).
        '''
        parents = self.popul.selection(self.noffspring)
        offspring = self.popul.recombination(parents, self.noffspring)
        self.evaluate(offspring)
        self.popul.reinsertion(offspring)

    def run(self) -> None:
        '''Method that runs the evolutionary algorithm cycle until finding the best solution or the number of iterations is reached.
        '''
        self.initPopul(self.indsize)
        self.evaluate(self.popul.indivs)
        self.bestsol = self.popul.bestSolution()
        for i in range(self.numits+1):
            self.iteration()
            bs = self.popul.bestSolution()
            if bs > self.bestsol:
                self.bestsol = bs
            # print("Iteration:", i, " ", "Best: ", self.bestsol)

    def printBestSolution(self):
        print("Best solution: ", self.bestsol.getGenes())
        print("Best fitness:", self.bestsol.getFitness())


def test():
    ea = EvolAlgorithm(100, 20, 50, 10)
    ea.run()
    ea.printBestSolution()


if __name__ == "__main__":
    test()
