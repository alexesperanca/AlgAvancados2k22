# -*- coding: utf-8 -*-
# Copyright 2022 by Group 7 (MSc Bioinformatics - University of Minho).  All rights reserved.

"""
This module provides the :class:`BWT` class that facilitates the analysis of big sequences and pattern discovery.
This class includes diverse strategies, such as:
	- Building of Burrows-Wheeler matrix that allows the user to a faster analysis of the initial sequence provided
	- Encountering of patterns
	- Original sequence faster retrieval
	- Building of a Suffix Array for match search with the BWT matrix. Faster search is performed with this method 
"""

import re
import pprint

class BWT:
	def __init__(self, seq: str):
		'''Initialization of the Burrows-Wheeler matrix construction 

		Parameters
		----------
		seq : str
			Sequence to be the model for BWT matrix
		'''
		self.seq = seq
		self._build_BWT()

	def _build_BWT(self):
		'''Burrows-Wheeler (BWT) matrix construction where the BWT line and Ordered line are defined and joined in a dictionary
		'''
		self.combinations = sorted([(self.seq[i:] + self.seq[:i], i) for i in range(len(self.seq))])	# Obter todas as sequências ordenadas
		print(self.combinations)
		b, pos = zip(*[(s[-1], p) for s, p in self.combinations])										# Obter apenas último membro das combinações (BWT) e nº de combinação
		print(b, pos)
		bwt = list(zip(b, pos))
		print(bwt)
		fun1 = self._nucl_table()
		self.bwt_line = [fun1(x) for x in [i[0] for i in bwt]]											# Obter linha BWT com a função _nucl_table que conta cada entrada dos caracteres
		print(self.bwt_line)
		fun2 = self._nucl_table()
		self.ord_line = [fun2(x) for x in sorted([i[0] for i in bwt])]
		print(self.ord_line)									# Obter linha Ordenada com a função _nucl_table tbm
		self.bwt_dic = {k: v for k, v in zip(self.bwt_line, self.ord_line)}								# Dicionário para aceder + facilmente a cada membro seguinte para recuperação da sequência
		print(self.bwt_dic)

	def _nucl_table(self) -> str:
		'''Auxiliary function that builds a dictionary with the occurences of characters provided to return the current occurrence of each character

		Returns
		-------
		str
			String with the current occurrence of the character given
		'''
		d = {}
		def _add(x: str) -> str:
			nonlocal d
			idx = d.get(x, 0)
			d[x] = idx + 1
			return x + str(idx)
		return _add

	def original_seq(self) -> str:
		'''Method that retrieves the original sequence from the BWT matrix built

		Returns
		-------
		str
			Original sequence
		'''
		l = "$0"
		s = ''
		while len(s) != len(self.seq) - 1:
			s += self.bwt_dic[l][0]
			l = self.bwt_dic[l]
		return s

	def find_pattern(self, pat: str) -> list:
		'''Method to find a pattern in the BWT matrix built. The algorithm walks through the matrix starting from the last character of the pattern. When reaches the end, the position corresponds to the local of the match of the pattern
		Mathod like "last_to_first"

		Parameters
		----------
		pat : str
			Pattern to match the BWT matrix

		Returns
		-------
		list
			Positions where the pattern matches
		'''
		dic2 = {v: k for k, v in self.bwt_dic.items()}
		l = [i for i in dic2.keys() if i[0] == pat[-1]]
		pos = len(pat) - 2
		while pos != -1:
			l = [dic2[i] for i in l if dic2[i][0] == pat[pos]]
			pos -= 1
		res = []
		for i, s in enumerate(self.ord_line):
			if s in l:
				res.append(i)
		return res

	def suffixarray(self) -> tuple:
		'''Method that constructs the suffix array. Retrieves a list with the initial positions of each ordered suffix

		Returns
		-------
		tuple
			Tuple of initial positions of each suffix
		'''
		d = {}
		for seq, pos in self.combinations:
			sufix = re.search(r'.*[^ACTG]', seq)		# Obtém a sequência até o símbolo "$" para adicionar no dicionário
			d[pos] = sufix.group()
		return tuple(d.keys())
