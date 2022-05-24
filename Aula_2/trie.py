import pprint

class Trie:
	def __init__(self, seq_list):
		self.tree = {}
		self.seq_list = seq_list
		self.ord = 0
		for seq in seq_list:
			self.insert(seq)

	def insert(self, seq):
		if seq not in self.seq_list: self.seq_list.append(seq)
		dic = self.tree
		for x in seq[0:]:
			if x not in dic:
				dic[x] = {}
			dic = dic[x]
		dic["#$#"] = self.ord
		self.ord += 1

	def _match(self, pat:str):
		seqs_present = []
		for seq in self.seq_list:
			if pat in seq:
				seqs_present.append(seq)
		return seqs_present

	def match(self, pat:str):
		seqs_present = self._match(pat)
		seqs = ", ".join(seqs_present)
		if seqs_present == []:
			print(pat, "is absent in the sequence list")
			return None
		else:
			print(f"{pat} present in " + seqs)
			return seqs_present

	def print_tree(self):
		pprint.pprint(self.tree, width = 1)

def test():
	test = Trie(["CTG", "CATA", "CAAGG"])
	test.insert("GGGA")
	test.print_tree()
	test.match("GG")

#test()