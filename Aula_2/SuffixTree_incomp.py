# -*- coding: utf-8 -*-
from trie import Trie

class SuffixTree(Trie):
    
    def __init__(self, seq:str):
        self.list_seqs = [seq[i:] for i in range(len(seq))]
        Trie.__init__(self, self.list_seqs)
        
    def add_suffix(self, p):
        self.insert(p)

    def find_pattern(self, pat):
        return self.match(pat)

    def get_leafes_below(self, node):
        res = []
        for n in self.tree[node]:
            print("N tem?", self.tree)
            string = node + n
            cur = n
            tree = self.tree[node][cur]
            while cur != "#$#":
                for new in tree:
                    cur = new
                    if cur != "#$#": string += cur
                    tree = tree[new]
                    break
            res.append(string)
        return res

    def repeats(self, pat):
        return len(self._match(pat))

def test1():
    print("* Teste 1 *\n")
    seq = "TACTA"
    st = SuffixTree(seq)
    st.print_tree()
    print()
    st.find_pattern("TA")
    st.find_pattern("ACG")
    print()
    st.add_suffix("GGAT")
    st.print_tree()
    print(st.get_leafes_below("C"))

def test2():
    print("\n* Teste 2 *\n")
    seq = "TACTA"
    st = SuffixTree(seq)
    st.find_pattern("TA")
    print(st.repeats("TA"))

test1()
test2()
        
            
    
    
