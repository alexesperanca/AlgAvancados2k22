import pprint

class Trie:
	def __init__(self, *seq_list):
		self.tree = {}
		self.seq_list = seq_list
		self.ord = 0
		for seq in seq_list:
			self.insert(seq)

	def insert(self, seq):
		dic = self.tree
		for x in seq[0:]:
			if x not in dic:
				dic[x] = {}
			dic = dic[x]
		dic["#$#"] = self.ord					# Incrementa no fim de cada sequência ter sido totalmente adicionada, o membro "#$#": 0
		self.ord += 1

	def match(self, pat):
		for seq in self.seq_list:
			if pat in seq:
				print(f"Present in {seq}\n")
				return True

	def print_tree(self):
		pprint.pprint(self.tree, width = 1)

test = Trie("CTG", "CATA", "CAAGG")
test.print_tree()
test.match("TG")