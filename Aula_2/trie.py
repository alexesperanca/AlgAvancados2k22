# import pprint

# class Trie:
# 	def __init__(self, *seq_list):
# 		self.tree = {}
# 		self.seq_list = seq_list
# 		self.ord = 0
# 		for seq in seq_list:
# 			self.insert(seq)

# 	def insert(self, seq):
# 		dic = self.tree
# 		for x in seq[0:]:
# 			if x not in dic:
# 				dic[x] = {}
# 			dic = dic[x]
# 		dic["#$#"] = self.ord					# Incrementa no fim de cada sequência ter sido totalmente adicionada, o membro "#$#": 0
# 		self.ord += 1

# 	def match(self, pat):
# 		for seq in self.seq_list:
# 			if pat in seq:
# 				print(f"Present in {seq}\n")
# 				return True

# 	def print_tree(self):
# 		pprint.pprint(self.tree, width = 1)

# test = Trie("CTG", "CATA", "CAAGG")
# test.print_tree()
# test.match("TG")

class Trie:
    def __init__(self):
        self.nodes = { 0:{} } # root node
        self.num = 0
    
    def print_trie(self):
        for k in self.nodes.keys():
            print (k, "->" , self.nodes[k])
    
    def add_node(self, origin, symbol):
        self.num += 1
        self.nodes[origin][symbol] = self.num
        self.nodes[self.num] = {}

    def add_pattern(self, p):
        pos = 0
        node = 0
        while pos < len(p):
            if p[pos] not in self.nodes[node].keys() :
                self.add_node(node, p[pos])
            node = self.nodes[node][p[pos]]
            pos += 1
    
    def trie_from_patterns(self, pats):
        for p in pats:
            self.add_pattern(p)

    def prefix_trie_match(self, text):
        pos = 0
        match = ""
        node = 0
        while pos < len(text):
            if text[pos] in self.nodes[node].keys() :
                node = self.nodes[node][text[pos]]
                match += text[pos]
                if self.nodes[node] == {}: return match
                else: pos += 1
            else: return None
        return None

    def trie_matches(self, text):
        # res = []
        for i in range(len(text)):
            m = self.prefix_trie_match(text[i:])
            if m != None:
                # res.append((i,m))
                print(f'The pattern {m} is in the position {i}.')
        # return res
    

# def test():
#     patterns = ["GAT", "CCT", "GAG"]
#     t = Trie()
#     t.trie_from_patterns(patterns)
#     t.print_trie()

def test2():
    patterns = ["AGAGAT", "AGC", "AGTCC", "CAGAT", "CCTA",
                "GAGAT", "GAT", "TC"]
    t = Trie()
    t.trie_from_patterns(patterns)
    # t.print_trie()
    t.prefix_trie_match("GAGATCCTA")
    t.trie_matches("GAGATCCTA")

# test()

test2()
