# -*- coding: utf-8 -*-
# Copyright 2022 by Group 7 (MSc Bioinformatics - University of Minho).  All rights reserved.

"""
This module provides the :class:`Trie` class and :class: `SuffixTree` class, allows the user to construct a wide tree with several ramifications enabling improved analysis.
`Trie` class includes diverse strategies for obtain pattern repeats, pattern recognition, and sequence addition, such as:
    - Incremention of sequences to the tree
    - Obtain mathes between a given sequence and patterns in the tree

`SuffixTrie` class has `Trie` as its parent class, allows the user to construct a wide tree of a given initial sequence.
This class includes diverse strategies for suffix addition, and tree analysis, such as:
	- Addition of suffix to the tree constructed
	- Obtainment of sequences by a singles leaf given to the class
	- Encountering of patterns in the tree crossing with a sequence
	- Sequence repeats identification in the tree
"""

import pprint

class Trie:
	'''Class that implements an tree structured with sequences provided and enables the user to iterate, cross sequences and recognize patterns 
	'''
	def __init__(self, seq_list:list):
		'''Creation of the tree structured as a dictionary with all the sequences given

		Parameters
		----------
		seq_list : list
			Sequences to be incremented in the tree (dictionary)
		'''
		self.tree = {}
		self.seq_list = seq_list
		self.ord = 0
		for seq in seq_list:
			self.insert(seq)

	def insert(self, seq: str):
		'''Method that inserts a given sequence in the tree

		Parameters
		----------
		seq : str
			Sequence to be incremented in the tree
		'''
		if seq not in self.seq_list: self.seq_list.append(seq)
		dic = self.tree
		for x in seq[0:]:
			if x not in dic:
				dic[x] = {}
			dic = dic[x]
		dic["#$#"] = self.ord
		self.ord += 1

	def _get_match(self, pat:str) -> str:
		'''Auxiliary function that searches for matches of the given pattern in the tree 

		Parameters
		----------
		pat : str
			Pattern to identify in the tree

		Returns
		-------
		str
			Returns a string representation of the hit obtained. Otherwise, None is returned
		'''
		pos = 0
		cur = pat[pos]
		tree = self.tree
		hit = ""
		while cur != "#$#":
			if cur in tree:
				if pat[pos] == cur: hit += cur
				else: return None
				tree = tree[cur]
				pos += 1
				for new in tree:
					cur = new
					try:
						if pat[pos] == new: break
					except:
						if hit != "" and hit in self.seq_list: return hit
						else: return None
			else: return None
		return hit

	def match(self, seq:str) -> list:
		'''Method that divides the pattern given in several fragmented sequences to run "_get_match" auxillary function to recognize if matches are found in the tree

		Parameters
		----------
		seq : str
			Sequence to run against the tree and identify if matches are found

		Returns
		-------
		list
			Return of list with the matches and positions of the sequences that were found
		'''
		res = []
		for i in range(len(seq)):
			m = self._get_match(seq[i:])
			if m != None and len(m) > 1: res.append((m, i))
		return res

	def trie_matches(self, pat:str):
		'''Method to print the matches and positions of the matches obtained from the "match" function

		Parameters
		----------
		pat : str
			Sequence to run against the tree and identify if matches are found
		'''
		matches = self.match(pat)
		if matches == []:
			print(pat, "is absent in the sequence list")
		else:
			for seq, pos in matches:
				print(f"{pat} contains {seq} in position {pos}")

	def print_tree(self):
		'''Pretty prints the tree
		'''
		pprint.pprint(self.tree, width = 1)

class SuffixTree(Trie):
	'''Subclass of Trie that depends on this to construct the suffix tree and use several of its methods

	Parameters
	----------
	Trie : class
		 Class that implements an tree structured with sequences provided and enables the user to iterate, cross sequences and recognize patterns 
	'''
	def __init__(self, seq:str):
		'''Costruction of the suffix tree from a given sequence. Sequence is deconstructed in several smaller sequences to run Trie class innitialization 

		Parameters
		----------
		seq : str
			Sequence model for the several sequences to tree construction
		'''
		self.list_seqs = [seq[i:] for i in range(len(seq))]
		Trie.__init__(self, self.list_seqs)
        
	def add_suffix(self, p:str):
		'''Suffix addition (sequence) to the tree previously constructed

		Parameters
		----------
		p : str
			Suffix (sequence) to add to the tree
		'''
		self.insert(p)

	def find_pattern_in_seq(self, seq:str) -> list:
		'''Discovery of the tree patterns and corresponding positions present in the provided sequence

		Parameters
		----------
		pat : str
			Sequence reference to cross with the tree

		Returns
		-------
		list
			List of hit patterns and positions
		'''
		if len(self.match(seq)) == 0:
			return 'No match!'
		else:
			return self.match(seq)

	def get_leafes_below(self, node:str) -> list:
		'''Obtaining of all the sequences from an given node of the tree. It requires to the node be an unique character since the tree nodes are unique chars from the sequences

		Parameters
		----------
		node : str
			Character to iterate over the tree

		Returns
		-------
		list
			List of strings obtained from the iteration of the node given
		'''
		assert node in self.tree.keys(), "Node inputted is not present in the tree"
		node = node.upper()
		res = []
		for n in self.tree[node]:
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

	def repeats(self, pat: str) -> int:
		'''Sequence to count the number of pattern occurrences happen in the tree 

		Parameters
		----------
		pat : str
			Sequence to cross from the tree and count the number of occurrences

		Returns
		-------
		int
			Number of pattern occurrences
		'''
		return len(self.match(pat))

## Testing of the functions (uncomment in the end of the file to see the class execution)

def test():
	print("\n* Teste 1 *\n")
	test = Trie(["CTG", "CATA", "CAAGG"])
	test.insert("GGGA")
	test.print_tree()
	test.trie_matches("CTGCATACAAGG")

def test2():
	print("\n* Teste 2 *\n")
	patterns = ["AGAGAT", "AGC", "AGTCC", "CAGAT", "CCTA",
	"GAGAT", "GAT", "TC"]
	t = Trie(patterns)
	t.print_tree()
	t.trie_matches("GAGATCCTA")

# def test3():
#     print("\n* Teste 3 *\n")
#     seq = "TACTA"
#     st = SuffixTree(seq)
#     st.print_tree()
#     print()
#     print(st.find_pattern("TATA"))
#     print(st.find_pattern("ACG"))
#     print()
#     st.add_suffix("GGAT")
#     st.print_tree()
#     print(st.get_leafes_below("C"))

# def test3():
#     print("\n* Teste 3 *\n")
#     seq = "TACTA"
#     seq2 = 'ACGT'
#     st = SuffixTree(seq)
#     st2 = SuffixTree(seq2)
#     # st.print_tree()
#     # print()
#     print(st.find_pattern("TATA"))
#     print(st2.find_pattern("AATTTCGACGTCGATTGAT"))
#     # print()
#     # st.add_suffix("GGAT")
#     # st.print_tree()
#     # print(st.get_leafes_below("C"))

def test4():
	print("\n* Teste 4 *\n")
	seq = "TACTA"
	st = SuffixTree(seq)
	print(st.find_pattern_in_seq("ACTA"))
	print(st.repeats("TA"))

#if __name__ == "__main__":
	#test()
	#test2()
	#test3()
	#test4()

# test3()