from EvolAlgorithm import EvolAlgorithm
from Popul import PopulInt, PopulReal
from MotifFinding import MotifFinding
from Motifs import Motifs


def createMatZeros(nl, nc):
    res = []
    for _ in range(0, nl):
        res.append([0]*nc)
    return res


def printMat(mat):
    for i in range(0, len(mat)):
        for j in range(len(mat[i])):
            print(f"{mat[i][j]:.3f}", end=' ')
        print()


class EAMotifsInt (EvolAlgorithm):
    def __init__(self, popsize, numits, noffspring, filename):
        self.motifs = MotifFinding()
        self.motifs.readFile(filename)
        indsize = len(self.motifs)
        EvolAlgorithm.__init__(self, popsize, numits, noffspring, indsize)

    def initPopul(self, indsize):
        maxvalue = self.motifs.seqSize(0) - self.motifs.motifSize
        self.popul = PopulInt(self.popsize, indsize, maxvalue, [])

    def evaluate(self, indivs):
        for i in range(len(indivs)):
            ind = indivs[i]
            sol = ind.getGenes()
            print(sol)
            fit = self.motifs.score(sol)
            ind.setFitness(fit)


class EAMotifsReal (EvolAlgorithm):
    
    def __init__(self, popsize, numits, noffspring, filename):
        self.motifs = MotifFinding()
        self.motifs.readFile(filename)
        indsize = len(self.motifs)*len(self.motifs.alphabet)       
        EvolAlgorithm.__init__(self, popsize, numits, noffspring, indsize)

    def initPopul(self, indsize):
        maxvalue = self.motifs.seqSize(0) - self.motifs.motifSize
        self.popul = PopulReal(self.popsize, indsize, 0.0, maxvalue, [])

    def evaluate(self, indivs):
        for i in range(len(indivs)):
            ind = indivs[i]
            sol = ind.getGenes()
            print(sol)
            ### fazer um motif através da lista de nº reais (nao pode ser pelos index)
            search = [int(x) for x in sol] #aqui ver a solução certa depois do motif correto
            
            ''' Notas:
            No método evaluate deverá, para cada indivíduo:
            - Construir uma matriz de dimensão A x L (a PWM)
            - Normalizar a matriz para que a soma de cada coluna seja 1
            - Determinar a posição mais provável do motif representado por essa PWM em cada sequência, construindo o vetor s
            - Calcular o score da solução s construída no passo anterior (este será a fitness)
            '''
            print(search)
            fit = self.motifs.score(search) #fazer uma 
            ind.setFitness(fit)


def test1():
    ea = EAMotifsInt(100, 1000, 50, "Trabalho/AlgAvancados2k22/Aula_4/Code/exemploMotifs.txt")
    ea.run()
    ea.printBestSolution()


def test2():
    ea = EAMotifsReal(100, 2000, 50, "Trabalho/AlgAvancados2k22/Aula_4/Code/exemploMotifs.txt")
    ea.run()
    ea.printBestSolution()


test1()
test2()
