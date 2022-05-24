from random import randint, random, shuffle, uniform
from typing import Type

class Indiv:

    def __init__(self, size: int, genes: list = [], lb: int = 0, ub: int = 1) -> None:
        '''Class to implement individuals with binary representations, with the atributes:
        - lb/ub (lower and upper limits of the range for representing genes)
        - fitness (stores fitness value for each individual

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
        if isinstance(solution, self.__class__):
            return self.genes.sort() == solution.genes.sort()
        return False

    def __gt__(self, solution):
        if isinstance(solution, self.__class__):
            return self.fitness > solution.fitness
        return False

    def __ge__(self, solution):
        if isinstance(solution, self.__class__):
            return self.fitness >= solution.fitness
        return False

    def __lt__(self, solution):
        if isinstance(solution, self.__class__):
            return self.fitness < solution.fitness
        return False

    def __le__(self, solution):
        if isinstance(solution, self.__class__):
            return self.fitness <= solution.fitness
        return False

    def __str__(self):
        return f"{str(self.genes)} {self.getFitness()}"

    def __repr__(self):
        return self.__str__()

    def setFitness(self, fit) -> None:
        self.fitness = fit

    def getFitness(self):
        return self.fitness

    def getGenes(self):
        return self.genes

    def initRandom(self, size) -> None:
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

    def one_pt_crossover(self, indiv2) -> tuple:
        '''Auxiliary method that makes a crossover between two individuals

        Parameters
        ----------
        indiv2 : type[Indiv]
            _description_

        Returns
        -------
        tuple
            _description_
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

    def __init__(self, size, genes=[], lb=0, ub=1):
        self.lb = lb
        self.ub = ub
        self.genes = genes
        self.fitness = None
        if not self.genes:
            self.initRandom(size)

    def initRandom(self, size):
        self.genes = []
        for _ in range(size):
            self.genes.append(randint(self.lb, self.ub))

    def mutation(self):
        s = len(self.genes)
        pos = randint(0, s-1)
        self.genes[pos] = randint(self.lb, self.ub)


class IndivReal(Indiv):
    
    def __init__(self, size, genes=[], lb=0.0, ub=1.0):
        self.lb = lb
        self.ub = ub
        self.genes = genes
        self.fitness = None
        if not self.genes:
            self.initRandom(size)

    def initRandom(self, size):
        self.genes = []
        for _ in range(size):
            self.genes.append(uniform(self.lb, self.ub))

    def mutation(self):
        s = len(self.genes)
        pos = randint(0, s-1)
        self.genes[pos] = uniform(self.lb, self.ub)
