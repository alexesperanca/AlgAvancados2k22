# -*- coding: utf-8 -*-

'''
autor: Alexandre Esperança
data: 25/03/2022
'''

'''
Passos:
	1. Construção de todas as possibilidades da string inserida
	2. Obtenção de uma lista com todos os membros finais das possibilidades e a sua posição
	3. Obtenção de uma lista com esses membros e a sua contabilização
	4. Mesmo do 3. mas com a lista ordenada
	5. Criação de um dicionário que interliga cada membro da lista criada em 3 com o da 4
'''

import re
import pprint

class BWT:
	def __init__(self, seq):
		self.seq = seq
		self._build_BWT()

	def _build_BWT(self):
		self.combinations = sorted([(self.seq[i:] + self.seq[:i], i) for i in range(len(self.seq))])
		#print(self.combinations)
		b, pos = zip(*[(s[-1], p) for s, p in self.combinations])
		bwt = list(zip(b, pos))
		fun1 = self._nucl_table()
		self.bwt_line = [fun1(x) for x in [i[0] for i in bwt]]
		#print(self.bwt_line)
		fun2 = self._nucl_table()
		self.ord_line = [fun2(x) for x in sorted([i[0] for i in bwt])]
		#print(self.ord_line)
		self.bwt_dic = {k: v for k, v in zip(self.bwt_line, self.ord_line)}
		#print(self.bwt_dic)

	def _nucl_table(self):
		d = {}
		def _add(x):
			nonlocal d
			idx = d.get(x, 0)
			d[x] = idx + 1
			return x + str(idx)
		return _add

	def original_seq(self):
		'''
		A dar diferente do PP idk why, mas coincide com o que o prof fez...
		'''
		l = "$0"
		s = ''
		while len(s) != len(self.seq) - 1:
			s += self.bwt_dic[l][0]
			l = self.bwt_dic[l]
		return s

	def find_pattern(self, pat):
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

	def suffixarray(self):
		d = {}
		for seq, pos in self.combinations:
			sufix = re.search(r'.*[^ACTG]', seq)
			d[pos] = sufix.group()
		return d


seq = "TAGACAGAGA$"
test = BWT(seq)
print(test.original_seq())
print(test.find_pattern("AGA"))
print(test.suffixarray())