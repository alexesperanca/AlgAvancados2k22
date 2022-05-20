from EvolAlgorithm import EvolAlgorithm
from Popul import PopulInt, PopulReal
from MotifFinding import MotifFinding


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

    def profile(self, lista):
        '''Creates a PWM profile of a given list of probabilistic values, with normnalized values (sum of each column(nucleotide) is 1).

        Parameters
        ----------
        profile : list
            List with values

        Returns
        -------
        matrix : list
            It returns list with a dictionary of the probabilistic profile
        '''
        mat = []
        for i in range(0, len(lista), len(self.motifs.alphabet)):
            dic = dict(zip(self.motifs.alphabet, lista[i:]))
            mat.append(dic)
        for n in self.motifs.alphabet:
            soma = 0
            for dici in mat:
                soma += dici[n]
            for dici in mat:
                dici[n] = dici[n]/soma
        return mat
    
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
    
    def seq_most_probable(self, seq: str, profile: list) -> str:
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
        list_seq = [seq[I:I+len(profile)] for I in range(len(seq) - len(profile) + 1)]
        probs = [self.prob_seq(s, profile) for s in list_seq]
        score_max = max(probs)
        ind = probs.index(score_max)
        seq = [list_seq[I] for I,p in enumerate(probs) if p == score_max][0]
    
        return seq, ind

    def evaluate(self, indivs):
        for i in range(len(indivs)):
            ind = indivs[i]
            sol = ind.getGenes()         
            #criar função auxiliar pwm ou pssm (soma de cada coluna é igual a 1 - normalizar valores)
            pwm = self.profile(sol) # criar o perfil pwm com a lista anterior (por ordem de 4 em 4 (se abcedario é DNA))
            # para cada sequencia, calcular a posição mais provavel do motif no pwm criado (vetor search)
            search = []            
            for x in range(len(self.motifs.seqs)):
                ## fazer função auxiliar para seq mais provavel (nao dá para utilizar função da Classe Motifs porque não é criada uma instancia, fazemos o perfil direto da lista)
                seq_prob, index = self.seq_most_probable(self.motifs.seqs[x], pwm)
                search.append(index)
            
            # calcular o score
            fit = self.motifs.score(search)
            ind.setFitness(fit)
            

def test1():
    ea = EAMotifsInt(100, 1000, 50, "AlgAvancados2k22/exemploMotifs.txt")
    ea.run()
    ea.printBestSolution()


def test2():
    ea = EAMotifsReal(100, 2000, 50, "AlgAvancados2k22/exemploMotifs.txt")
    ea.run()
    ea.printBestSolution()


test1()
test2()
