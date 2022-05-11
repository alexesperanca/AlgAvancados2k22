# -*- coding: utf-8 -*-
import re

class Automata:
    
    def __init__(self, alphabet, pattern):
        self.numstates = len(pattern) + 1
        self.alphabet = alphabet
        self.transitionTable = {}
        self.pattern = pattern
        self.buildTransitionTable(pattern)        
    
    def buildTransitionTable(self, pattern):
        for q in range(self.numstates):
            for a in self.alphabet:
                if q % 2 == 0:
                    if pattern.find(a) == 0:
                        self.transitionTable[(q, a)] = self._obtain_transition(q)
                    
                    elif pattern.find(a) == 1:
                        self.transitionTable[(q, a)] = max(0, q - 2)

                else:
                    if pattern.find(a) == 0:
                        if q - 2 < 0:
                            self.transitionTable[(q, a)] = q
                        else:
                            self.transitionTable[(q, a)] = q - 2

                    elif pattern.find(a) == 1:
                        self.transitionTable[(q, a)] = self._obtain_transition(q)
       
    def _obtain_transition(self, q):
        if q + 1 == self.numstates:
            return q - 1
        else:
            return q + 1

    def printAutomata(self):
        print ("States: " , self.numstates)
        print ("Alphabet: " , self.alphabet)
        print ("Transition table:")
        for k in self.transitionTable.keys():
            print (k[0], ",", k[1], " -> ", self.transitionTable[k])

    def nextState(self, current, symbol):
        return self.transitionTable[[current, symbol]]
        
    def applySeq(self, seq):
        q = 0
        res = [q]
        for i in seq:
            res.append(self.transitionTable[(q, i)])
            q = self.transitionTable[(q, i)]
        return res
        
    def occurencesPattern(self, text):
        q = 0  
        m = re.search(self.pattern, text)
        res = m.span()
        return res

def overlap(s1, s2):
    maxov = min(len(s1), len(s2))
    for i in range(maxov,0,-1):
        if s1[-i:] == s2[:i]: return i
    return 0
               
def test():
    auto = Automata("AC", "ACA")
    auto.printAutomata()
    print (auto.applySeq("CACAACAA"))
    print (auto.occurencesPattern("CACAACAA"))

test()

#States:  4
#Alphabet:  AC
#Transition table:
#0 , A  ->  1
#0 , C  ->  0
#1 , A  ->  1
#1 , C  ->  2
#2 , A  ->  3
#2 , C  ->  0
#3 , A  ->  1
#3 , C  ->  2
#[0, 0, 1, 2, 3, 1, 2, 3, 1]
#[1, 4]


